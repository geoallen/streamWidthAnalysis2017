<b>streamWidthAnalysis2017</b>

<b> Description:</b><br>

This repository contains R code used to analyze the data and produce figures presented in:
Allen et al. (2018) "Similarity of stream width distributions across headwater systems" published in <i>Nature Communications</i> (DOI: <a href="https://www.nature.com/articles/s41467-018-02991-w">10.1038/s41467-018-02991-w</a>).

<b>Installation:</b>

- You will have to have a working `R` installation to use these routines. Download and install `R` according to the instructions on [the R website](https://www.r-project.org/)
- Download this repository into your working directory.
- Download raw field and GIS data into the same working directory (DOI: <a href="https://zenodo.org/record/1034385#.WpSVhxPwZE4">10.5281/zenodo.1034384</a>).
- To check and make sure your installation is correct, compare your installation folder with the file tree shown at the bottom of this README.

<b>Usage:</b>

- Open a Terminal and `cd` into the directory where you installed the data and scripts.
- Type `./smallStreamsAnalysis.R` to run the script, reproducing figures and analysis. The directories `/figures` and `/tables` will be created, and contain all the products of the analysis.

<b>Contributing:</b>

If you wish to modify or use this code for your own work, please contact the authors and contributors to this repository. We want to know what you're doing, and are willing to offer advice and some minimal assistance.

<b>Credits:</b>

Owner and research lead: @geoallen

Code contributors:
- @ericbarefoot

<b>File structure:</b>

To insure that the code runs correctly, and your file structure is correct, please make sure that when you download all the data, it matches the following structure:
```
.
├── EDfig2_EDtab3_GOF.R
├── EDtab1_EDtab2_catchment_attributes.R
├── fig1_widthMap.R
├── fig2_distributions.R
├── fig3_widthModel.R
├── figures
├── LICENSE
├── locationStreamSurveys
│   ├── blueduck.csv
│   ├── caribou.csv
│   ├── elder.csv
│   ├── kings.csv
│   ├── sagehen.csv
│   ├── stony.csv
│   └── v40.csv
├── README.md
├── repeatStreamSurveys
│   ├── stony_subcatchment_20151027.csv
│   ├── stony_subcatchment_20151209.csv
│   ├── stony_subcatchment_20160202.csv
│   ├── stony_subcatchment_20160214.csv
│   ├── stony_subcatchment_20160304a.csv
│   └── stony_subcatchment_20160304b.csv
├── smallStreamsAnalysis.R
├── streamWidthModelOutput
│   ├── modeledWidthTab4_3_Blue Duck.csv
│   ├── modeledWidthTab4_3_Elder.csv
│   ├── modeledWidthTab4_3_Kings.csv
│   ├── modeledWidthTab4_3_Sagehen.csv
│   └── modeledWidthTab4_3_V40.csv
└── tables
```
