## Mortality burdens attributable to temperature and seasonal variation
**R script** files to estimate cause-specific mortality burdens attributable to ambient temperature and seasonal variation in Spain.<br>
<br>
Tobías A, Iñíguez C, Royé D, Hashizume M, Madaniyazi L. Comparing mortality burdens attributable to ambient temperature and seasonal variation across major causes of death in Spain. <b>Environmental Epidemiology</b> 2025 (accepted for publication).

---

The **R script** files below reproducde step by step the analysis using an open access dataset for the city of Valencia (Spain). 
<br>
* **main_analysis.R** reproduces the examples from the published manuscript. 

* and the following **ancillary functions** also need to be uploaded: 

    * **valencia.csv** - example dataset for the city of Valencia (from the compilation of open access time-series datasets for studying temperature-mortality association by <a href="https://pubmed.ncbi.nlm.nih.gov/39071964/" target="_blank">Tobías et al. 2024</a>).  
    * **cyclic.R** - function to generate the clyclic spline to fit seasonality (adapted from the R packge DLNM by <a href="https://pubmed.ncbi.nlm.nih.gov/22003319/" target="_blank">Gasparrini 2017</a>).  
    * **findmax.R** - function to estimate the peak of seasonality (adapted from the R function to estimate the minimum the minimum mortality temperature by <a href="https://pubmed.ncbi.nlm.nih.gov/27748681/" target="_blank">Tobías et al. 2017</a>). 
    * **findmin.R** - function to estimate the trhough of seasonality (adapted from the R function to estimate the minimum the minimum mortality temperatureby <a href="https://pubmed.ncbi.nlm.nih.gov/27748681/" target="_blank">Tobías et al. 2017</a>). 
    * **attrdl.R** - function to estimate the attributable fraction of temperature-response (R function to estimate attributable measures from DLNM by <a href="https://pubmed.ncbi.nlm.nih.gov/24758509/" target="_blank">Gasparrini and Leone 2014</a>).
    * **attrs.R** - function to estimate the attributable fraction of seasonality (adapted from the R function to estimate attributable measures from DLNM by <a href="https://pubmed.ncbi.nlm.nih.gov/24758509/" target="_blank">Gasparrini and Leone 2014</a>). 
