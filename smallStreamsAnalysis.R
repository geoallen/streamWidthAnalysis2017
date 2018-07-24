# allen-headwater-stream-widths-distributions-study

# contains R code used to analyze the data presented in Allen et al., "Similarity of stream width distributions in headwater catchments"

# to run:
# download all data a scripts into the same working directory
# open top level R script: "smallStreamsAnalysis.R"
# specify workingDir to your working directory path
#run code to reproduce figures and analysis

# fieldwork_analysis.R
# George Allen, Oct 2017

#############################################################################################
# Load required packages

require(here) # for easy relative path management.

#############################################################################################
# Define input parameters:

# contains pointers to run:

workingDir = '/home/eric/Dropbox/programs/git_repos/streamWidthAnalysis2017'

# data table file names:

data_folders = c('locationStreamSurveys', 'repeatStreamSurveys', 'streamWidthModelOutput')

locNames = c('kings', 'sagehen', 'elder', 'caribou', 'v40', 'blueduck', 'stony') # these were not the same as the data files in the repo.
repNames = c('stony_subcatchment_20151027', 'stony_subcatchment_20151209',
           'stony_subcatchment_20160202', 'stony_subcatchment_20160214',
           'stony_subcatchment_20160304a', 'stony_subcatchment_20160304b')
loc_paths = here('locationStreamSurveys', paste0(locNames, '.csv'))
rep_paths = here('repeatStreamSurveys', paste0(repNames, '.csv'))
inTabPaths = c(loc_paths, rep_paths)
# inQRecordPaths = here('', fNames, 'discharge_records/') # is this essential? I can't immediately see where it gets used but these data are not in the repo.
if (F %in% file.exists(inTabPaths)){message("Field Data CSV files are missing")}

# figure labels:
tabNames = c('Kings', 'Sagehen', 'Elder', "Caribou", "V40", "Blue Duck", "Stony",
             "2015-10-27", "2015-12-09", "2016-02-02",
             "2016-02-14", "2016-03-04a", "2016-03-04b")

#############################################################################################
# Figure Scripts:

if(!file.exists(here('figures'))) {dir.create(here('figures')); message("'figures' directory did not exist. Creating...")} # checks for "figures" folder and creates it if it doesnt exist.

# Figure 1 - stream width map generator:
source(here('fig1_widthMap.R'))
pdfOut = here('figures', 'widthMap.pdf')
fig1_widthMap(inTabPaths, tabNames, pdfOut) # we have a memory issue with this one. It crashed R for me...

# This one ^ refers to GIS data that isn't provided in the repo...

# Figure 2 - stream width distributions:
# this one put functions last. they should be first, or put somewhere else. Otherwise the function doesn't know what they are and it throws an error.
# Was this working on your machine?
source(here('fig2_distributions.R'))
pdfOut = here('figures', 'fig2_distributions.pdf')
fig2_distributions(inTabPaths, tabNames, pdfOut)

# Figure 3 - modeled stream widths:
# this one has hardwired paths to some result.
# it also puts functions last. They can be last if they're outside the function call...
source(here('fig3_widthModel.R'))
pdfOut = here('figures', 'fig3_widthModel4_3.pdf')
csvOut = here('tables', 'modeledWidthTab4_3.csv')
fig3_widthModel(inTabPaths, tabNames, csvOut, pdfOut)

# Figure 4 was produced in Adobe Illustrator

#############################################################################################
# Extended Data Figures and Tables:

if(!file.exists(here('tables'))) {dir.create(here('tables')); message("'tables' directory did not exist. Creating...")} # checks for "tables" folder and creates it if it doesnt exist.

# ED Figure 1 was produced in Adobe Illustrator

# ED table 1 and table 2 - catchment attributes:
source(here('EDtable_catchment_attributes.R'))
csvOut = here('tables', 'EDtable1.csv')
EDtable_catchment_attributes(inTabPaths, tabNames, csvOut, workingDir)

# ED Figure 2 and ED Table 3 - quantify GOF for distributions:
source(here('EDfig2_EDtab3_GOF.R'))
modTabDir = here('tables')
pdfOut = here('figures', 'EDfig2_GOF.pdf')
csvOut = here('tables', 'EDtable3_GOF.csv')
EDfig2_EDtab3_GOF(inTabPaths, modTabDir, tabNames, pdfOut, csvOut)

# ED table 4 - efflux calculation:
source(here('EDtable3.R'))
csvOut = here('tables', 'EDtable3.csv')
EDtable4(inTabPaths, tabNames, csvOut, workingDir)
