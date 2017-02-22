# smallStreamsAnalysis.R
# George Allen, Feb 2017  

# contains R code used to analyze the data and create figures and tables
# presented in Allen et al., "Similarity of stream width distributions in headwater catchments" 

# to run:
# download all data and R scripts into the same working directory
# open top level R script: "smallStreamsAnalysis.R" 
# specify workingDir as your working directory path (below)
# select and run codes to reproduce figures and analysis

#############################################################################################
# Define input parameters. This section contains pointers 
# and constants to run scripts:


workingDir = 'E:/misc/2015_09_01_Small_Stream_Width_Fieldwork'

# data table file names:
fNames = c('konza', 'sagehen', 'angelo', 'caribou', 'v40', 'blueduck', 'stony',
           'stony_subcatchment_20151027', 'stony_subcatchment_20151209',
           'stony_subcatchment_20160202', 'stony_subcatchment_20160214',
           'stony_subcatchment_20160304a', 'stony_subcatchment_20160304b')
inTabPaths = paste0(workingDir, '/', fNames, '/field/', fNames, '_field_dat.csv')
inQRecordPaths = paste0(workingDir, '/', fNames, 'discharge_records/')
if (F %in% file.exists(inTabPaths)){message("Field Data CSV files are missing")}

# figure labels:
tabNames = c('Kings', 'Sagehen', 'Elder', "Caribou", "V40", "Blue Duck", "Stony",
             "2015-10-27", "2015-12-09", "2016-02-02",
             "2016-02-14", "2016-03-04a", "2016-03-04b")


#############################################################################################
# this section loads codes used to produce figures in the main manuscript:

# Figure Scripts:

# Figure 1 was created using ArcMap. 

# Alternative Fig. 1 - stream width map generator:
source(paste0(workingDir, '/R/Fig1_widthMap.R'))
pdfOut = paste0(workingDir, '/figures/Fig1_widthMap.pdf')
Fig1_widthMap(inTabPaths, tabNames, pdfOut)

# Figure 2 - stream width distributions:
source(paste0(workingDir, '/R/fig2_distributions.R'))
pdfOut = paste0(workingDir, '/figures/fig2_distributions.pdf')
fig2_distributions(inTabPaths, tabNames, pdfOut)

# Figure 3 - modeled stream widths:
source(paste0(workingDir, '/R/fig3_widthModel.R'))
pdfOut = paste0(workingDir, '/figures/fig3_widthModel4_3.pdf')
csvOut = paste0(workingDir, '/tables/modeledWidthTab4_3.csv')
fig3_widthModel(inTabPaths, csvOut, pdfOut)

# Figure 4 was produced in Adobe Illustrator

#############################################################################################
# Extended Data Figures and Tables:

# ED Figure 1 was produced in Adobe Illustrator 

# ED table 1 and table 2 - catchment attributes:
source(paste0(workingDir, '/R/EDtab1_EDtab2_catchment_attributes.R'))
csvOut = paste0(workingDir, '/tables/EDtab1_EDtab2.csv')
EDtab1_EDtab2_catchment_attributes(inTabPaths, tabNames, csvOut, workingDir)

# ED Figure 2 and ED Table 3 - quantify GOF for distributions:
source(paste0(workingDir, '/R/ED_fig2_tab3_GOF.R'))
pdfOut = paste0(workingDir, '/figures/ED_fig2_GOF.pdf')
csvOut = paste0(workingDir, '/tables/EDtable3_GOF.csv')
ED_fig2_tab3_GOF(inTabPaths, tabNames, pdfOut, csvOut)

# ED table 4 - CO2 efflux calculation:
source(paste0(workingDir, '/R/EDtab4_co2_flux.R'))
csvOut = paste0(workingDir, '/tables/EDtab4_co2_flux.csv')
EDtab4_co2_flux(inTabPaths, tabNames, csvOut, workingDir)