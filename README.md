<b>streamWidthAnalysis2017</b>

<b> Description:</b><br>

This repository contains R code used to analyze the data and produce figures presented in:
Allen et al. (2018) "Similarity of stream width distributions across headwater systems" published in <i>Nature Communications</i> (DOI: <a href="https://www.nature.com/articles/s41467-018-02991-w">10.1038/s41467-018-02991-w</a>).

<b>Installation:</b>

- You will have to have a working `R` installation to use these routines. Download and install `R` according to the instructions on [the R website](https://www.r-project.org/)
- Either download this repository into your working directory, or clone it from git into the desired directory with:
```
git clone https://github.com/geoallen/streamWidthAnalysis2017.git
```
- Download raw field and GIS data into the same working directory from zenodo (DOI: <a href="https://zenodo.org/record/1034385#.WpSVhxPwZE4">10.5281/zenodo.1034384</a>).
- To check and make sure your installation is correct, compare your installation folder with the file tree shown at the bottom of this README.

<b>Usage:</b>

Using your preferred method, users can open up an R session either in RStudio or the terminal, and `source` the top level script. To do this:
- open R in your preferred way, and set the working directory to the installation folder.
- run the following:
```
source('smallStreamsAnalysis.R')
```

Optionally, if users want finer grained control, and are comfortable with the command line, the following also works:

- Open a Terminal and `cd` into the directory where you installed the data and scripts.
- Enter `./smallStreamsAnalysis.R` on UNIX or just `smallStreamsAnalysis` in the terminal to run the script, reproducing figures and analysis. The directories `/figures` and `/tables` will be created, and contain all the products of the analysis.
- Users on UNIX systems (MacOS and Linux) who have GNU `make` installed can simply run `make all` to generate all analyses. `make clean` will restore the repository to a pre-run configuration.

The command line interface for these scripts accepts the following options, which can be accessed in the command line with `./smallStreamsAnalysis.R --help`

```
usage:  smallStreamsAnalysis.r [options]

options:
  --fig1_2       Generate figures 1 and 2.
  --fig3         Generate figure 3.
  --fig4         Generate draft figure 4 and table 1.
  --EDtab1_2     Generate supplementary tables 1 and 2.
  --EDtab3_fig2  Generate supplementary table 3 and figure 2.
  --help         Usage.
```

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
