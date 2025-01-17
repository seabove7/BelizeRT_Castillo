# Gene expression plasticity facilitates acclimatization of a long-lived Caribbean coral across divergent reef environments

#### **Authors:** Karl D. Castillo<sup>#</sup>, Colleen B. Bove<sup>#</sup>, Annabel M. Hughes, Maya E. Powell, Justin B. Ries, Sarah W. Davies

*<sup>#</sup> denotes co-first authors* 

**Citation:**  Castillo KD, Bove CB *et al.*Gene expression plasticity facilitates acclimatization of a long-lived Caribbean coral across divergent reef environments. *Sci Rep* (2024). https://doi.org/10.1038/s41598-024-57319-0

#### **Zenodo DOI:** <a href="https://zenodo.org/badge/latestdoi/248840574"><img src="https://zenodo.org/badge/248840574.svg" alt="DOI"></a>


### **Abstract:**  
Local adaptation can increase fitness under stable environmental conditions. However, in rapidly changing environments, compensatory mechanisms enabled through plasticity may better promote fitness. Climate change is causing devastating impacts on coral reefs globally and understanding the potential for adaptive and plastic responses is critical for reef management. We conducted a four-year, three-way reciprocal transplant of the Caribbean coral *Siderastrea siderea* across forereef, backreef, and nearshore populations in Belize to investigate the potential for environmental specialization versus plasticity in this species. Corals maintained high survival within forereef and backreef environments, but transplantation to nearshore environments resulted in high mortality, suggesting that nearshore environments present strong environmental selection. Only forereef- sourced corals demonstrated evidence of environmental specialization, exhibiting the highest growth in the forereef. Gene expression profiling 3.5 years post-transplantation revealed that transplanted coral hosts exhibited profiles more similar to other corals in the same reef environment, regardless of their source location, suggesting that transcriptome plasticity facilitates acclimatization to environmental change in *S. siderea*. In contrast, algal symbiont (*Cladocopium goreaui*) gene expression showcased functional variation between source locations that was maintained post-transplantation. Our findings suggest limited acclimatory capacity of some *S. siderea* populations under strong environmental selection and highlight the potential limits of coral physiological plasticity in reef restoration.

<br/>

---

### Code Written By:
* The calcification, survival, and temperature scripts (*BelizeRT_CalcSurvData.Rmd*) were written by Colleen B. Bove (colleenbove@gmail.com)
* The gene expression scripts (*RT_Host_Sym_GE.R*) and SNP analyses (directory: *PopStructure_SNPs*) were written by Sarah W. Davies (daviessw@gmail.com)
* The ITS2 scripts (*reciprocal_transplants_ITS2.Rmd*) were written by Annabel M. Hughes (amhughes@bu.edu) and Colleen B. Bove (colleenbove@gmail.com)

**For questions regarding data or analyses of each section should be directed to the appropriate author.**

---

### Repository contains the following:
1. **[Code](https://github.com/seabove7/BelizeRT_Castillo/tree/main/Code)**
    * **CCAplast_function.R** -- Custom function to calculate gene expression plasticity based on CCA called in the *RT_Host_Sym_GE.R* script.
   * **Ssid_SNP_pipeline.txt** -- Pipeline run on Boston University's SCC to identify SNPs in host and symbiont samples (used in *PopStructure_SNPs* directory)
    * **uniHeatmap.R** -- Function for creating heatmaps in the *RT_Host_Sym_GE.R* script.
2. **[Data](https://github.com/seabove7/BelizeRT_Castillo/tree/main/Data)**
    * **219_20220803T065042_DBV_20220803T143825.profiles.absolute.abund_and_meta.txt** -- ITS2 *Symportal* raw output of symbiont type counts.
    * [GeneExpression_data](https://github.com/seabove7/BelizeRT_Castillo/tree/main/Data/GeneExpression_data) -- Directory containing data and output files necessary for the *RT_Host_Sym_GE.R* script for GE analyses.
    * **ITS2_symportal_RT.csv** -- ITS2 absolute counts cleaned up for use in *reciprocal_transplants_ITS2.Rmd* script.
    * **RT_calc_boot.csv** -- The bootstrapped calcification rates per sample over the entire experimental timeline calculated in the *BelizeRT_CalcSurvData.Rmd* script.
    * **RT_calcification.csv** -- The raw calcification rates per sample used in the *BelizeRT_CalcSurvData.Rmd* script.
    * **RT_SSID_microPrep_ITS2_meta.csv** -- Sample metadata for matching well ID to sample ID for ITS2 analyses (used in *reciprocal_transplants_ITS2.Rmd*).
    * **RT_survival_table.csv** -- Resulting summary table of coral samples across time (created in *BelizeRT_CalcSurvData.Rmd*).
    * **RT_survival.csv** -- Raw sample survival data per individual used in *BelizeRT_CalcSurvData.Rmd* script.
    * **RT_TP_calc_boot.csv** -- The bootstrapped calcification rates per sample per time point calculated in the *BelizeRT_CalcSurvData.Rmd* script.
    * **RTE_Temperature_Nov2014_Oct2015.csv** -- Temperature data collected *in situ* at each reef location used for temperature analyses in the *BelizeRT_CalcSurvData.Rmd* script.
3. **[Figures](https://github.com/seabove7/BelizeRT_Castillo/tree/main/Figures)**
    * **Figure1_map_temp** -- Figure 1 from the manuscript (with both PDF and PNG versions)
    * **Figure2_survival_calc** -- Figure 2 from the manuscript (with both PDF and PNG versions)
    * **Figure3_RTE_GE** -- Figure 3 from the manuscript (with Illustrator, PDF, and PNG versions)
    * [Supplemental_Figures](https://github.com/seabove7/BelizeRT_Castillo/tree/main/Figures/Supplemental_Figures) -- Directory containing all supplemental figures from the manuscript
4. **[GO_enrichment_Host](https://github.com/seabove7/BelizeRT_Castillo/tree/main/GO_enrichment_Host)** -- Directory containing all GO enrichment input and output files with the accompanying script for the coral host analyses
5. **[GO_enrichment_Symb](https://github.com/seabove7/BelizeRT_Castillo/tree/main/GO_enrichment_Symb)** -- Directory containing all GO enrichment input and output files with the accompanying script for the algal symbiont analyses
6. **[PopStructure_SNPs](https://github.com/seabove7/BelizeRT_Castillo/tree/main/PopStructure_SNPs)** -- Directory containing all input and output files with the accompanying script for performing the SNP analysis of population structure for the coral host and algal symbionts. 
7. **BelizeRT_CalcSurvData** -- R markdown script and html output with all code and analyses of temperature, survival, and calcification included in manuscript 
8. **RT_Host_Sym_GE.R** -- R script with all the code and analyses of gene expression of the host and symbionts included in manuscript
9. **reciprocal_transplants_ITS2** -- R markdown script and html output with all code and analyses of algal symbiont diversity included in manuscript
