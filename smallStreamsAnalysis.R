#!/usr/bin/env Rscript

# allen-headwater-stream-widths-distributions-study

# contains R code used to analyze the data presented in Allen et al., "Similarity of stream width distributions in headwater catchments"

# to run:
# download all data a scripts into the same working directory
# open top level R script: "smallStreamsAnalysis.R"
# specify workingDir to your working directory path
#run code to reproduce figures and analysis

# fieldwork_analysis.R
# George Allen, Oct 2017

########################################################
# install and load  necessary libraries

packages = c('here', 'devtools', 'foreign', 'MASS', 'shapefiles')

ipak = function(pkg){
  new_pkg = pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new_pkg)) {
    install.packages(new_pkg, dependencies = TRUE, repos = 'http://cran.us.r-project.org')
  }
}

gitPackages = c('docopt')

gitPackagePaths = c('docopt/docopt.R')

gitPak = function(pkg, path) {
  require(devtools)
  new_pkg = pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new_pkg)) {
    install_github(path, dependencies = TRUE)
  }
}

ipak(packages)
gitPak(gitPackages, gitPackagePaths)

#############################################################################################
# Load required packages

require(here) # for easy relative path management.
require(docopt) # for passing arguments from command line.

# command line usage information

'
usage:  smallStreamsAnalysis.r [options]

options:
  --fig1_2       Generate figures 1 and 2.
  --fig3         Generate figure 3.
  --fig4         Generate figure 4 and table 1.
  --EDtab1_2     Generate supplementary tables 1 and 2.
  --EDtab3_fig2  Generate supplementary table 3 and figure 2.
  --EDtab4       Generate supplementary table 4.
  --help         Usage.
' -> doc

opts = docopt(doc)

#############################################################################################
# Define input parameters:

# data table file names:

data_folders = c('locationStreamSurveys', 'repeatStreamSurveys', 'streamWidthModelOutput')

dataNames = c('konza', 'sagehen', 'angelo', 'caribou', 'v40', 'blueduck', 'stony',
            'stony_subcatchment_20151027', 'stony_subcatchment_20151209',
            'stony_subcatchment_20160202', 'stony_subcatchment_20160214',
            'stony_subcatchment_20160304a', 'stony_subcatchment_20160304b')

inTabPaths = here('data', dataNames, 'field',  paste0(dataNames, '_field_dat.csv'))

if (F %in% file.exists(inTabPaths)){message("Field Data CSV files are missing")}

# figure labels:
tabNames = c('Konza', 'Sagehen', 'Angelo', "Caribou", "V40", "Blue Duck", "Stony",
             "2015-10-27", "2015-12-09", "2016-02-02",
             "2016-02-14", "2016-03-04a", "2016-03-04b")

# make sure folders exist for outputs

if(!file.exists(here('figures'))) {dir.create(here('figures')); message("'figures' directory did not exist. Creating...")} # checks for "figures" folder and creates it if it doesnt exist.
if(!file.exists(here('tables'))) {dir.create(here('tables')); message("'tables' directory did not exist. Creating...")} # checks for "tables" folder and creates it if it doesnt exist.

#############################################################################################
# Figure Scripts:

# Figures 1 and 2- stream width map generator:
# Makes a preliminary draft figure, which was later annotated manually in Adobe Illustrator and ArcGIS.
# Figure 2 is comprised of the bottom half of the draft figure.
if (opts$fig1_2 | all(!unlist(opts))) {
  message('compiling figures 1 and 2...')
  source(here('fig1_2_widthMap.R'))
  pdfOut = here('figures', 'fig1_2_widthMap_draft.pdf')
  fig1_2_widthMap(inTabPaths, tabNames, pdfOut)
}

# Figure 3 - stream width distributions:
if (opts$fig3 | all(!unlist(opts))) {
  message('compiling figure 3...')
  source(here('fig3_distributions.R'))
  pdfOut = here('figures', 'fig3_distributions.pdf')
  fig3_distributions(inTabPaths, tabNames, pdfOut)
}

# Figure 4 - modeled stream widths:
if (opts$fig4 | all(!unlist(opts))) {
  message('compiling figure 4 and table 1...')
  source(here('fig4_widthModel.R'))
  pdfOut = here('figures', 'fig4_widthModel4_3.pdf')
  csvOut = here('tables', 'modeledWidthTab4_3')
  fig4_widthModel(inTabPaths, tabNames, csvOut, pdfOut)
}

# Figure 5 was produced in Adobe Illustrator

#############################################################################################
# Extended Data Figures and Tables:

# ED Figure 1 was produced in Adobe Illustrator

# ED table 1 and table 2 - catchment attributes:
if (opts$EDtab1_2 | all(!unlist(opts))) {
  message('compiling supplementary tables 1 and 2...')
  source(here('EDtab1_EDtab2_catchment_attributes.R'))
  csvOut = here('tables', 'EDtable1.csv')
  EDtab1_EDtab2_catchment_attributes(inTabPaths, tabNames, csvOut, workingDir)
}

# ED Figure 2 and ED Table 3 - quantify GOF for distributions:
if (opts$EDtab3_fig2 | all(!unlist(opts))) {
  message('compiling supplementary table 3 and supplementary figure 2...')
  source(here('EDfig2_EDtab3_GOF.R'))
  modTabDir = here('tables')
  pdfOut = here('figures', 'EDfig2_GOF.pdf')
  csvOut = here('tables', 'EDtable3_GOF.csv')
  EDfig2_EDtab3_GOF(inTabPaths, modTabDir, tabNames, pdfOut, csvOut)
}
# is there a script for this?
# # ED table 4 - efflux calculation:
# source(here('EDtable3.R'))
# csvOut = here('tables', 'EDtable3.csv')
# EDtable4(inTabPaths, tabNames, csvOut, workingDir)
