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
# inQRecordPaths = paste0(workingDir, '/', fNames, 'discharge_records/') # is this essential? I can't immediately see where it gets used but these data are not in the repo.
if (F %in% file.exists(inTabPaths)){message("Field Data CSV files are missing")}

# figure labels:
tabNames = c('Kings', 'Sagehen', 'Elder', "Caribou", "V40", "Blue Duck", "Stony",
             "2015-10-27", "2015-12-09", "2016-02-02",
             "2016-02-14", "2016-03-04a", "2016-03-04b")

#############################################################################################
# Figure Scripts:

if(!file.exists(here('figures'))) {dir.create(here('figures')); message("'figures' directory did not exist. Creating...")} # checks for "figures" folder and creates it if it doesnt exist.

# Figure 1 - stream width map generator:
source(paste0(workingDir, '/fig1_widthMap.R'))
pdfOut = paste0(workingDir, '/figures/widthMap.pdf')
fig1_widthMap(inTabPaths, tabNames, pdfOut) # we have a memory issue with this one. It crashed R for me...

# Figure 2 - stream width distributions:
source(paste0(workingDir, '/fig2_distributions.R'))
pdfOut = paste0(workingDir, '/figures/fig2_distributions.pdf')
fig2_distributions(inTabPaths, tabNames, pdfOut)

# Figure 3 - modeled stream widths:
source(paste0(workingDir, '/fig3_widthModel.R'))
pdfOut = paste0(workingDir, '/figures/fig3_widthModel4_3.pdf')
csvOut = paste0(workingDir, '/tables/modeledWidthTab4_3')
fig3_widthModel(inTabPaths, tabNames, csvOut, pdfOut)

# Figure 4 was produced in Adobe Illustrator


#############################################################################################
# Extended Data Figures and Tables:

if(!file.exists(here('tables'))) {dir.create(here('tables')); message("'tables' directory did not exist. Creating...")} # checks for "tables" folder and creates it if it doesnt exist.

# ED Figure 1 was produced in Adobe Illustrator

# ED table 1 and table 2 - catchment attributes:
source(paste0(workingDir, '/EDtable_catchment_attributes.R'))
csvOut = paste0(workingDir, '/tables/EDtable1.csv')
EDtable_catchment_attributes(inTabPaths, tabNames, csvOut, workingDir)

# ED Figure 2 and ED Table 3 - quantify GOF for distributions:
source(paste0(workingDir, '/EDfig2_EDtab3_GOF.R'))
modTabDir = paste0(workingDir, '/tables')
pdfOut = paste0(workingDir, '/figures/EDfig2_GOF.pdf')
csvOut = paste0(workingDir, '/tables/EDtable3_GOF.csv')
EDfig2_EDtab3_GOF(inTabPaths, modTabDir, tabNames, pdfOut, csvOut)

# ED table 4 - efflux calculation:
source(paste0(workingDir, '/EDtable3.R'))
csvOut = paste0(workingDir, '/tables/EDtable3.csv')
EDtable4(inTabPaths, tabNames, csvOut, workingDir)
