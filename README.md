# Estimating genus-specific effects of non-native honey bees and urbanization on wild bee communities in the Eastern United States 

### Gabriela M. Quinlan, Jeffrey W. Doser, Melanie Kammerer, Christina M. Grozinger 

### In Review

### Please contact the first author for questions about the code or data used in the analysis: Gabriela M. Quinlan (gmq5021@psu.edu)

---------------------------------

## Abstract

Non-native species have the potential to detrimentally affect native species through resource competition, disease transmission, and other forms of antagonism. The western honey bee (Apis mellifera) is one such species that has been widely introduced beyond its native range for hundreds of years. There are strong concerns in the United States and other countries around the strain high-density, managed honey bee populations could pose to already imperiled wild bee communities. While there is some experimental evidence of honey bees competing with wild bees for resources, few studies have connected landscape-scale honey bee apiary density with down-stream consequences for wild bee communities. Here, using a unique dataset from Maryland, US and joint species distribution models, we provide the largest scale, most phylogenetically resolved assessment of non-native honey bee density effects on wild bee abundance to date. As Maryland primarily consists of urban beekeeping, we also assessed the relative impact of developed land on wild bee communities. We found poor evidence for negative effects of apiary density on the wild bee community overall, with the strongest evidence for negative effects of urban beekeeping on late-season, specialist genera (e.g., long-horned bees) and small, ground nesting, season-long foragers (e.g., green sweat bees). We also found positive effects of developed land for some genera including invasive Anthidium and other urban garden-associated genera. We additionally discuss methodological insights based on sampling efficiency of different methods (hand netting, pan trapping, vane trapping) across genera. Notably, the effect of apiary density and developed land were small relative to sampling method, so these effects should not be over interpreted. These findings offer some of the best evidence for the effects of managed honey bees on wild bee communities and we discuss several avenues to ameliorate potentially detrimental effects of urban beekeeping on the most imperiled wild bee groups. 

## Repository Directory

### [code/](./code/)

Contains all code to process data, fit models, and generate figures and results shown in the manuscript. note that the raw data processing files will not run successfully using the files on GitHub, as many of the raw data sets are too large for inclusion on GitHub, or are proprietary (apiary locations). All files will run successfully if run in the order specified by the file names starting from file "3a-main-null.R".  

+ `1a-reclassCDL.R`: script to reclass the Cropland Data Layer into the five categories shown in Figure 1 that were used in subsequent analyses (will not run successfully). 
+ `1b-mdClean.R`: script to clean the Maryland data from Sam Droege and colleagues at the USGS (will not run successfully).
+ `1c-get-apiary-surface.R`: script to extract the metric of apiary density across a 100x100m grid of Maryland to generate Figure 1 (will not run successfully).  
+ `1d-spAbundance-data-prep.R`: prepares data into format for fitting a JSDM in `spAbundance` (will not run successfully).
+ `2-get-inits.R`: script to extract initial values from a previous model fit (will not run successfully).
+ `3a-main-null.R`: script to run the null JSDM without either developed landcover or apiary density as covariates in the model. 
+ `3b-main-developed.R`: script to run a JSDM with developed landcover as a covariate in the model.
+ `3c-main-apiary.R`: script to run a JSDM using spAbundance using apiary density as a covariate in the model.
+ `3d-main-devel-apiary.R`: script to run a JSDM using apiary density and developed landcover as covariates in the model.
+ `4a-combineChains.R`: a function to process three chains run separately in `spAbundance` into a single object.
+ `4b-post-process-chains.R`: combines the individual chains that were run in parallel into a single `spAbundance` model object.
+ `4c-extract-samples.R`: script to extract the relevant information from the full `spAbundance` model objects and save in a smaller, more manageable file.
+ `5-convergence-gof-assessment.R`: assess model convergence and goodness of fit for all candidate models.
+ `6-hierarchical-partitioning.R`: script to perform hierarchical partitioning using the four candidate model results.
+ `7-summary.R`: script to summarize all analysis results and generate all figures provided in the manuscript.

### [data/](./data/)

Contains the data sets formatted for fitting models in `spAbundance`. Note that many of the raw data sets are too large for inclusion on GitHub, or are proprietary (apiary locations). 

+ `data-spAbundance-flat.csv`: a flat file of the data needed for fitting models in `spAbundance`. 
+ `reclassMDDC.tiff`: the reclassed CDL data used as covariates in the models.
+ `spAbundance-data.rda`: R data file containing the data necessary for fitting models in `spAbundance`.
+ `developed-for-violin-plot.rda`: some data used to generate the violin plot in Figure S3. 

### [figures/](.figures/)

Contains all figures included in the manuscript and supplemental information.

### [results/](.results/)

Contains some results files generated from the JSDMs fit in `spAbundance`. Note that most results files are too large for inclusion on GitHub, but can be generated by running the scripts in the `code` directory.

+ `beta-subset-prob-results.csv`: a subset of results from the JSDMs used to generate information in Tables 1 and 2 in the manuscript. 
+ `hp-results.csv`: the results from the hierarchical partitioning.
+ `inits-devel-apiary.rda`: initial values needed for running `3d-main-devel-apiary.R`. 
+ `inits-apiary.rda`: initial values needed for running `3c-main-apiary.R`. 
+ `inits-developed.rda`: initial values needed for running `3b-main-developed.R`.
+ `inits-null.rda`: initial values needed for running `3a-main-null.R`. 
+ `effect-sizes-apiary-pred.csv`: estimates of wild bee genera relative abundance at values of high and low apiary density (used to generate Figure S4). 
