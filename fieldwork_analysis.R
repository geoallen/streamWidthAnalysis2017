# allen-headwater-stream-widths-distributions-study

# contains R code used to analyze the data presented in Allen et al., "Similarity of stream width distributions in headwater catchments" 

# to run:
# download all data a scripts into the same working directory
# open top level R script: "smallStreamsAnalysis.R" 
# specify workingDir as your working directory path
# select and run code to reproduce figures and analysis
  
# fieldwork_analysis.R
# George Allen, Feb 2017  

#############################################################################################
# Define input parameters:  

# contains pointers to run:

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
# Figure Scripts:

# Figure 1 - stream width map generator:
source(paste0(workingDir, '/R/widthMap.R'))
pdfOut = paste0(workingDir, '/figures/widthMap.pdf')
widthMap(inTabPaths, tabNames, pdfOut)

# Figure 2 - stream width distributions:
source(paste0(workingDir, '/R/fig2_distributions.R'))
pdfOut = paste0(workingDir, '/figures/fig2_distributions.pdf')
fig2_distributions(inTabPaths, tabNames, pdfOut)

# Figure 3 - modeled stream widths:
source(paste0(workingDir, '/R/fig3_widthModel.R'))
pdfOut = paste0(workingDir, '/figures/fig3_widthModel4_3.pdf')
csvOut = paste0(workingDir, '/tables/modeledWidthTab4_3')
fig3_widthModel(inTabPaths, tabNames, csvOut, pdfOut)

# Figure 4 was produced in Adobe Illustrator


#############################################################################################
# Extended Data Figures and Tables:

# ED Figure 1 was produced in Adobe Illustrator 

# ED table 1 and table 2 - catchment attributes:
source(paste0(workingDir, '/R/EDtable_catchment_attributes.R'))
csvOut = paste0(workingDir, '/tables/EDtable1.csv')
EDtable_catchment_attributes(inTabPaths, tabNames, csvOut, workingDir)

# ED Figure 2 and ED Table 3 - quantify GOF for distributions:
source(paste0(workingDir, '/R/EDfig2_EDtab3_GOF.R'))
modTabDir = paste0(workingDir, '/tables')
pdfOut = paste0(workingDir, '/figures/EDfig2_GOF.pdf')
csvOut = paste0(workingDir, '/tables/EDtable3_GOF.csv')
EDfig2_EDtab3_GOF(inTabPaths, modTabDir, tabNames, pdfOut, csvOut)

# ED table 4 - efflux calculation:
source(paste0(workingDir, '/R/EDtable3.R'))
csvOut = paste0(workingDir, '/tables/EDtable3.csv')
EDtable4(inTabPaths, tabNames, csvOut, workingDir)








io=0;if(io==1){


# join shapefile XYZ data to field tables:
source(paste0(workingDir, '/R/shapefile2fielddat.R'))
pdfOut = paste0(workingDir, '/figures/shapefile2fieldTabMap.pdf')
shapefile2fielddat(inTabPaths, tabNames, fNames, pdfOut, workingDir)




# join field data to shapefile tables:
source(paste0(workingDir, '/R/fieldDat2shapefile.R'))
fieldDat2shapefile(inTabPaths, tabNames, workingDir)


# normal_distributions_byOrder.R
source(paste0(workingDir, '/R/normal_distributions_byOrder.R'))
pdfOut = paste0(workingDir, '/figures/normal_distributions_byOrder.pdf')
normal_distributions_byOrder(inTabPaths, tabNames, pdfOut, 'run')



# mean_vs_SD_width_plot scatter:
source(paste0(workingDir, '/R/mean_vs_SD_width_plot.R'))
pdfOut = paste0(workingDir, '/figures/mean_vs_SD_width_plot.pdf')
mean_vs_SD_width_plot(inTabPaths, tabNames, pdfOut, 'run')


# width vs downstream distance plots:
source(paste0(workingDir, '/R/width_vs_Distance_bySegment.R'))
pdfOut = paste0(workingDir, '/figures/width_vs_Distance_bySegment.pdf')
width_vs_Distance_bySegment(inTabPaths, tabNames, pdfOut, 'run')


# width vs downstream distance with Histograms plots:
source(paste0(workingDir, '/R/width_vs_Distance_bySegment_withHistograms.R'))
pdfOut = paste0(workingDir, '/figures/width_vs_Distance_bySegment_withHistograms.pdf')
width_vs_Distance_bySegment_withHistograms(inTabPaths, tabNames, pdfOut, 'run')


# compare stream length of NHD flowlines to field survey:
source(paste0(workingDir, '/R/NHDPlusV21_2_fieldStreamLengthCompare.R'))
pdfOut = paste0(workingDir, '/figures/NHDPlusV21_2_fieldStreamLengthCompare.pdf')
csvOut = paste0(workingDir, '/tables/NHDPlusV21_2_fieldStreamLengthCompare.csv')
NHDPlusV21_2_fieldStreamLengthCompare(inTabPaths, tabNames, pdfOut, csvOut, workingDir, 'run')


# mean_vs_SD_width_plot scatter:
source(paste0(workingDir, '/R/mean_vs_SD_width_plot.R'))
pdfOut = paste0(workingDir, '/figures/mean_vs_SD_width_plot.pdf')
mean_vs_SD_width_plot(inTabPaths, tabNames, pdfOut, 'run')


# spline fits one plot:
source(paste0(workingDir, '/R/splines_onePlot.R'))
pdfOut = paste0(workingDir, '/figures/splines_onePlot.pdf')
splines_onePlot(inTabPaths, tabNames, pdfOut, 'run')


# density kernel one plot:
source(paste0(workingDir, '/R/kernel_onePlot.R'))
pdfOut = paste0(workingDir, '/figures/kernel_onePlot_20cm.pdf')
bwSeq = seq(3, 50, 1)
meanModeW = vector()
for (i in 1:length(bwSeq)){
  x = kernel_onePlot(inTabPaths, tabNames, pdfOut, bwSeq[i], 'run')
  meanModeW = append(meanModeW, mean(x))
}


# histogram by stream order:
source(paste0(workingDir, '/R/simple_hist_byOrder.R'))
pdfOut = paste0(workingDir, '/figures/stacked_hist.pdf')
simple_hist_byOrder(inTabPaths, tabNames, pdfOut, 'run')


# width by order box plot:
source(paste0(workingDir, '/R/boxplot_widthByOrder.R'))
pdfOut = paste0(workingDir, '/figures/boxplot_widthByOrder.pdf')
boxplot_widthByOrder(inTabPaths, tabNames, pdfOut, 'run')


# Width histograms by segment:
source(paste0(workingDir, '/R/widthHists_bySegment.R'))
pdfOut = paste0(workingDir, '/figures/widthHists_bySegment.pdf')
widthHists_bySegment(inTabPaths, tabNames, pdfOut, 'run')


# Width CDFs one plot:
source(paste0(workingDir, '/R/CDFs_onePlot.R'))
pdfOut = paste0(workingDir, '/figures/CDFs_onePlot.pdf')
CDFs_onePlot(inTabPaths, tabNames, pdfOut, 'run')

 
# segment length table generator:
source(paste0(workingDir, '/R/segLengthTabGen.R'))
csvOut = paste0(workingDir, '/tables/segmentLength_noTrim.csv')
segLengthTabGen(inTabPaths, tabNames, csvOut, 'run')


# test whether distribution is LogNorm or Power Law:
source(paste0(workingDir, '/R/logNorm_distrib_tester.R'))
pdfOut = paste0(workingDir, '/figures/logNormal_distrib_tester.pdf')
logNorm_distrib_tester(inTabPaths, tabNames, pdfOut, 'run')





# segment width and length variability: 



# segment length and width box plot:
#source(paste0(workingDir, '/R/boxplot_seg_widthAndLengthByOrder.R')
#pdfOut = paste0(workingDir, '/figures/boxplot_seg_widthAndLengthByOrder.pdf'
#boxplot_seg_widthAndLengthByOrder(inTabPaths, tabNames, pdfOut, 'norun')


}




