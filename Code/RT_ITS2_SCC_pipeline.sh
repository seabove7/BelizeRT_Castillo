## Castillo/Bove RT 16S and ITS2 
# Pipeline by Colleen B Bove
#------------------------------

### Move files from TUFC to SCC

## Log on to your scc account using ssh: ssh bovec@scc1.bu.edu
# log on to TUCF's ftp 
# adding -i will get rid of prompts later on (replace XXX with the provided IP address)
ftp -i XXXXX

## you will be prompted for username and password (220711-0481M_Sarah_Davies_7102)
# Username: YOUR.NAME.HERE
# Password: YOUR.PASSWORD.HERE
## you might have to type the following to get ftp to prompt you fo username
user

# NOTE: tab fills are deactivated in ftp

## to pull multiple files to remote server use
mget -r Demultiplex_Summary-220711-0481M_Sarah_Davies_7102.html bovec@scc1.bu.edu:/projectnb/davieslab/bove/raw_seqs/RT_seqs

mget -r *.html bovec@scc1.bu.edu:/projectnb/davieslab/bove/raw_seqs/RT_seqs/fastq_Lane1



# if you get an error like WARNING! 1433673 bare linefeeds received in ASCII mode 
# this means there was a problem entering binary mode and even if you successfully transfer files, it will screw up any zipped files. 
# to get around this before doing mget do 
binary 
# then do mget as described above

# if you get an error like 426 Failure writing network stream
# this is likely because the remote server took too long to download the file, you can try #when you have better internet, or try downloading one file at a time. 


## Be sure to do TUFT's recommended integrity check after unzipping your files!! 
# selecting just the.gz files to compare for now
grep ".fastq.gz" md5sum.txt > md5sum_gz.txt 

# compare md5sum outputs by first generating one
md5sum *fastq.gz > checklist.txt

# and then check that they match
diff fastq_Lane1/checklist.txt md5sum_gz.txt

# check the number of files (158 total; means 79 PE files with the 79th being the unknown)
ls | grep .fastq. | wc -l


### Copy the fastq files to working folder
cp /projectnb/davieslab/bove/raw_seqs/RT_seqs/fastq_Lane1/*.gz /projectnb/davieslab/bove/RT_seqs/fastq_files





# --------------------------------------------------

#### Sample pre-processing

## Installing BBtools
# see instructions: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/

# move to home directory (/projectnb/davieslab/bove) from local computer
scp -r BBMap_38.96.tar.gz bovec@scc1.bu.edu:/projectnb/davieslab/bove

# unpack the package (will create new directory called bbmap)
tar xopf BBMap_38.96.tar.gz

# To test the installation run stats.sh against the PhiX reference genome (included in the BBTools download):
/projectnb/davieslab/bove/bbmap/stats.sh in=/projectnb/davieslab/bove/bbmap/resources/phix174_ill.ref.fa.gz


## Adaptors for the barcodes (Hyb_F*_i5/Hyb_R*_i7) used in this run
# Forward adaptor:  AATGATACGGCGACCAC
# Reverse adaptor:  CAAGCAGAAGACGGCATAC

# These along with the reverse complements of the adaptors are saved as a fasta file:
nano adaptors.fasta

# ---
>forward
AATGATACGGCGACCAC
>forwardrc
GTGGTCGCCGTATCATT
>reverse
CAAGCAGAAGACGGCATAC
>reverserc
GTATGCCGTCTTCTGCTTG
# ---


### Cleaning and the separating samples ----*

## Note: Illumina should have cut these out already, normal if you don't get any

## Renaming the files
# make a list of sample names (will save MP_TV_A1_S1_L001_R1_001.fastq to A1 -- due to -f 3)
ls *R1_001.fastq | cut -d '_' -f 3 > samples.list

# remove unneeded text from file names (do this for .gz files)
nano rename.sh
# ---
#!#/bin/bash

for file in $(cat samples.list); do  
mv ${file}_*R1*.fastq ${file}_R1.fastq; 
mv ${file}_*R2*.fastq ${file}_R2.fastq; 
done
# ---

## Remove reads that still have the adaptor sequence (they should not be there)
nano adapt_rm.sh
# ---
#!#/bin/bash

for file in $(cat samples.list); do 
/projectnb/davieslab/bove/bbmap/bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; 
done &>bbduk_NoIll.log
# ---


# You can check how many were removed like this: (I had none so nothing was removed)
grep "Total Removed:" bbduk_NoIll.log 
# example of what output should look like (x each sample)
# Total Removed:          	0 reads (0.00%) 	0 bases (0.00%)


# moving these into 'adapt_rm_files' direction
mv fastq_files/*NoIll.fastq adapt_rm_files/


## Remove first 4 bases (degenerate primers created them)
nano rm_4base.sh
# ---
#!#/bin/bash

for file in $(cat /projectnb/davieslab/bove/RT/fastq_files/samples.list); do 
/projectnb/davieslab/bove/bbmap/bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq ftl=4 out1=${file}_R1_NoIll_No4N.fastq out2=${file}_R2_NoIll_No4N.fastq; 
done &>bbduk_No4N.log
# ---

# moving these into 'degen_rm' direction
mv adapt_rm_files/*No4N.fastq degen_rm/




### ITS2 samples ----*
## NOTE: If submitting to SymPortal, we want to keep the primers at request (but I am leaving that code here for reference)

## ITS2 primers (SYM_VAR_REV/SYM_VAR_5.8S2) used in this run

# Forward ITS2 primer:  GAATTGCAGAACTCCGTGAACC
# Reverse ITS2 primer:  CGGGTTCWCTTGTYTGACTTCATGC

# saving these as a file to input
nano ITS2_primers.fasta

# ---
>forward
GAATTGCAGAACTCCGTGAACC
>reverse
CGGGTTCWCTTGTYTGACTTCATGC
# ---

## only keeping reads that start with the ITS2 primer
# higher k = more reads removed, but cannot surpass k=20 or 21
nano its_pull.sh
# ---
#!#/bin/bash

for file in $(cat /projectnb/davieslab/bove/RT/fastq_files/samples.list); do 
/projectnb/davieslab/bove/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq k=15 restrictleft=21 ref=ITS2_primers.fasta outm1=${file}_R1_NoIll_No4N_ITS.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_ITS.fastq outu2=${file}_R2_check.fastq; 
done &>bbduk_ITS.log

#for file in $(cat /projectnb/davieslab/bove/RT/fastq_files/samples.list); do 
#/projectnb/davieslab/bove/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq k=15 restrictleft=21 literal=GAATTGCAGAACTCCGTGAACC,CGGGTTCWCTTGTYTGACTTCATGC outm1=${file}_R1_NoIll_No4N_ITS.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_ITS.fastq outu2=${file}_R2_check.fastq; 
#done &>bbduk_ITS.log
# ---


## Now you can check the log file (bbduk_ITS.log) to see # counts pulled with the primers ('Contaminants') -- these are saved as {file}_R1_NoIll_No4N_ITS.fastq
head -25 bbduk_ITS.log

# Input:                          55520 reads             13713440 bases.
# Contaminants:                   33506 reads (60.35%)    8275982 bases (60.35%)
# Total Removed:                  33506 reads (60.35%)    8275982 bases (60.35%) *** These are the counts with ITS2 primers
# Result:                         22014 reads (39.65%)    5437458 bases (39.65%)


# moving these into 'ITS2_out_files' direction
mv degen_rm/*ITS.fastq ITS2_out_files/
mv degen_rm/*check.fastq ITS2_out_files/


## To get the total counts per raw file:
grep -c '@' *.fastq >> raw_counts.txt


## ----- Stopping here with the ITS2 files since SymPortal requests primers remain in the reads

# copy the 'ITS.fastq' files to the ITS_SymPortal_29Jun22 directory (wihtin final) and then rename them
cp ITS2_out_files/*ITS.fastq final_samples/ITS_SymPortal_29Jun22

for file in *.fastq; do
mv -- "$file" "${file/_NoIll_No4N/}"
done



## Use cutadapt to remove primer:
module load cutadapt/4.1

nano cutprime.sh
# ---
#!#/bin/bash

module load cutadapt/4.1

for file in $(cat /projectnb/davieslab/bove/RT/fastq_files/samples.list); do
cutadapt -g GAATTGCAGAACTCCGTGAACC -a GCATGAAGTCAYACAAGWGAACCCG -G CGGGTTCWCTTGTYTGACTTCATGC -A GGTTCACGGAGTTCTGCAATTC -n 2 --discard-untrimmed -o ${file}_ITS2_R1.fastq -p ${file}_ITS2_R2.fastq ${file}_R1_NoIll_No4N_ITS.fastq ${file}_R2_NoIll_No4N_ITS.fastq;
done &> clip.log
# ---

# this will save the files as the original name (i.e., A1_ITS2_R1.fastq)


## Move the final files to the 'ITS2_sym_15July22' directory
mv ITS2_out_files/*R1.fastq final_samples/ITS2_sym_15July22/
mv ITS2_out_files/*R2.fastq final_samples/ITS2_sym_15July22/





### 16S samples ----*

## 16S primers (16s_806r/16s_515f) used in this run

# Forward 16S primer:  GTGYCAGCMGCCGCGGTAA
# Reverse 16S primer:  GGACTACNVGGGTWTCTAAT

# saving these as a file to input (the 16s_primers.fasta are the bad primers we double checked)
nano 16S_primers.fasta

# ---
>forward
GTGYCAGCMGCCGCGGTAA
>reverse
GGACTACNVGGGTWTCTAAT
# ---

## only keeping reads that start with the 16S primer
# higher k = more reads removed, but cannot surpass k=20 or 21
nano 16s_pull.sh
# ---
#!#/bin/bash

for file in $(cat /projectnb/davieslab/bove/RT/fastq_files/samples.list); do 
/projectnb/davieslab/bove/bbmap/bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq k=7 restrictleft=20 ref=16S_primers.fasta outm1=${file}_R1_NoIll_No4N_16S.fastq outu1=${file}_R1_16Scheck.fastq outm2=${file}_R2_NoIll_No4N_16S.fastq outu2=${file}_R2_16Scheck.fastq; 
done &>bbduk_16S.log
# ---

# moving these into '16S_out_files' direction
mv degen_rm/*16S.fastq 16S_out_files/
mv degen_rm/*16Scheck.fastq 16S_out_files/


## Use cutadapt to remove primer:
module load cutadapt/4.1

nano cutprime.sh
# ---
#!#/bin/bash

module load cutadapt/4.1

for file in $(cat /projectnb/davieslab/bove/RT/fastq_files/samples.list)
do
cutadapt -g GTGYCAGCMGCCGCGGTAA -a ATTAGAWACCCVNGTAGTCC -G GGACTACNVGGGTWTCTAAT -A TTACCGCGGCMGCTGYCAC -n 2 --discard-untrimmed -o ${file}_16S_R1.fastq -p ${file}_16S_R2.fastq ${file}_R1_NoIll_No4N_16S.fastq ${file}_R2_NoIll_No4N_16S.fastq
done &> clip.log
# ---

# this will save the files as the original name (i.e., A1_16S_R1.fastq)


## Move the final files to the '16S_15July22' directory
mv 16S_out_files/*R1.fastq final_samples/16S_15July22/
mv 16S_out_files/*R2.fastq final_samples/16S_15July22/





#------------------ Notes on what the bbduk.sh and cutadapt commands mean 

### Notes on bbduk.sh use:
# in1=            Main input file
# in2=            Input for 2nd read of pairs in a different file 
# ref=	          Comma-delimited list of reference files.
# out1=	          Write reads here that do not contain kmers matching the database.
# out2=	          Use this to write 2nd read of pairs to a different file.
# k=              Kmer (nucleotide sequence of a certain length) length used for finding contaminants. Contaminants shorter than k will not be found.  k must be at least 1.
# restrictleft=   If positive, only look for kmer matches in the leftmost X bases.
# literal=        Comma-delimited list of literal reference sequences.
# ftl=            Trim the leftmost n bases  
# outm1=          Write reads here that fail filters. In default kfilter mode, this means any read with a matching kmer. In any mode, it also includes reads that fail filters such as minlength, mingc, maxgc, entropy, etc.  In other words, it includes all reads that do not go to 'out'. (These are the ones with the primers to proceed with)
# outu1=          Reads that pass all filtering criteria.
# outm2=          (same as outm1 for the second file)
# outu2=          (same as outu1 for the second file)


### Notes on cutadapt use:
# -g Sequence of an adapter ligated to the 5' end (forward primer)
# -a Sequence of an adapter ligated to the 3' end (RC of the reverse primer)
# -G 5' adapter to be removed from R2 (reverse primer)
# -A 3' adapter to be removed from R2 (RC of the forward primer)
# -o Write trimmed reads to FILE.
# -p paired-output FILE
# -max-n 0 Discard reads with more than COUNT 'N' bases (means 0 Ns allowed)





### -------------------------------------------- ###
# Pipeline for uploading files to NCBI from BU SCC #
### -------------------------------------------- ###

## Moved the RT samples to load into the SRA here: /projectnb/davieslab/bove/RT/final_samples/ITS2_sym_15July22/SRA_RT_submission


# 1. Navigate to the source folder where the files for submission are (/projectnb/davieslab/bove/RT/final_samples/ITS2_sym_15July22/SRA_RT_submission);

# 2. Establish an FTP connection using the credentials below (replace XXXXXX with your provided address):
    ftp -i
    open XXXXXXX

    # Address: YOUR.FTP.ADDRESS
	# Username: YOUR.USER.NAME
	# Password: YOUR.PASSWORD

# 3. Navigate to your account folder:
    # From the command line use the 'cd' command (replace XXXXXX with your provided directory):
    cd uploads/XXXXX

# 4. Create a subfolder (required!) with a meaningful name:
    mkdir RT_ITS2
    #mkdir RT_TagSeq (this is the one with the TagSeq files)

# 5. Navigate to the target folder you just created:
    cd RT_ITS2

# 6. Copy ALL files into the target folder:
    mput *
