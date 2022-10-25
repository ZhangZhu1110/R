# We have a single cell gene expression profile with 113 samples and around 100 genes.
# They are derived from 10 patients and 4 vascular sites.
# Therefore, sach sample has 2 labels, patients and locations.
# I aim to demonstrate the data in different angle of views.
# The following codes exhibit the comparison of different visualization strategy.
# The corresponding report can be browsed in https://github.com/ZhangZhu1110/Courses-Reports-in-NTU/blob/main/Data%20visualization.pdf


library(ggplot2)
library(pheatmap)
setwd("C:/Users/zhang/Desktop/Data Visualization/Project/")
data = read.csv('CTC.csv', header = TRUE)
row.names(data) <- data[, 1]
data = data[, -1]


# single gene heterogeneity in location
# facet is not necessary
ggplot(data = data) + geom_boxplot(aes(y = RPL23AP7, group = location, fill = location))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background =
                                                                                                  element_blank(), axis.line = element_line(colour = "black"))

ggplot(data = data) + geom_boxplot(aes(y = RPL23AP7, group = location, fill = location))+ facet_wrap(~location, nrow = 2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background =
                                                                                                                                    element_blank(), axis.line = element_line(colour = "black"))
# linear correlation between 2 genes
ggplot(data = data, mapping = aes(x = MMRN1, y = JAM3)) + geom_point()
ggplot(data = data, mapping = aes(x = MMRN1, y = JAM3)) + geom_point() + geom_smooth()
#ggplot(data = data, mapping = aes(x = MMRN1, y = JAM3)) + geom_point(aes(color = location)) + geom_smooth() + facet_wrap(~location, nrow = 2)

#single gene heterogeneity in patient / location
# without facet_wrap
ggplot(data = data) + geom_boxplot(aes(y = LRMP, fill = patient))
ggplot(data = data) + geom_boxplot(aes(y = LRMP, fill = location))
# with facet_wrap
ggplot(data = data) + geom_boxplot(aes(y = FOS, fill = patient)) + facet_wrap(~patient, nrow = 2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot(data = data) + geom_boxplot(aes(y = FOS, fill = location)) + facet_wrap(~location, nrow = 2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# density plot and box plot
ggplot(data = data) + aes(x = LRMP, fill = location) + geom_density(alpha = 0.3)
ggplot(data = data) + aes(x = LRMP, fill = patient) + geom_density(alpha = 0.3)

# emphasize location
ggplot(data = data) + geom_boxplot(aes(y = LSM5, fill = location)) + facet_grid(location~patient)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# emphasize patient
ggplot(data = data) + geom_boxplot(aes(y = LSM5, fill = patient)) + facet_grid(location~patient)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


# investigate the difference between 2 genes among location / patient
# mapping continuous variable to color, size and shape
ggplot(data = data, aes(location, TCTA)) + geom_point(aes(colour = MYL4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggplot(data = data, aes(location, TCTA)) + geom_point(aes(size = MYL4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
# the x axis and y axis should be gene(continuous variable), rather than location(categorical variable)
ggplot(data = data, aes(MYL4, TCTA)) + geom_point(aes(shape = location)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggplot(data = data, aes(MYL4, TCTA)) + geom_point(aes(color = location)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggplot(data = data, aes(MYL4, TCTA)) + geom_point(aes(size = location)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
# better solution
ggplot(data = data, aes(MYL4, TCTA)) + geom_point(aes(color = location)) + facet_wrap(~location, nrow = 2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggplot(data = data, aes(MYL4, TCTA)) + geom_point(aes(color = patient)) + facet_wrap(~patient, nrow = 2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))

# the sample size of location and patients
# map categorical variable to size, color and shape
ggplot(data = data) + geom_count(aes(x = location, y = patient)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggplot(data = data) + geom_count(aes(x = location, y = patient), color = 'blue') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggplot(data = data) + geom_count(aes(x = location, y = patient, shape = location))
ggplot(data = data) + geom_count(aes(x = location, y = patient, size = location))
# mapping patient which has more unique value to aesthetic factor, worse result 
ggplot(data = data) + geom_count(aes(x = location, y = patient, shape = patient))
ggplot(data = data) + geom_count(aes(x = location, y = patient, size = patient))
#better solution
#emphasize patient
ggplot(data = data) + geom_bar(aes(x = location, fill = location)) + facet_wrap(~patient, nrow = 2)
#emphasize location
ggplot(data = data) + geom_bar(aes(x = patient, fill = patient)) + facet_wrap(~location, nrow = 2)


# location1
cols_remain<-c('TCERG1',
               'FOS',
               'ANGPT1',
               'ADORA2A',
               'LSM5',
               'RABL2A',
               'TRBV7-4',
               'MIR584',
               'C20orf96',
               'GSTM1',
               'NDRG3',
               'WDR61',
               'SUSD3',
               'XRN2',
               'CLK1',
               'MYL1',
               'TMSB4Y',
               'RAMP3',
               'RSPH9'
)

#location 2

cols_remain<-c('RPL23AP7',
               'TSPAN7',
               'AMBP',
               'SNORD3B-1',
               'SNORD3B-2',
               'C20orf96',
               'COMMD5',
               'RNASE6',
               'SNORD3D',
               'FMO5',
               'CHMP4C',
               'CPB2',
               'LUC7L3',
               'AKR1C1',
               'PLA2G16',
               'ANXA1',
               'S100A13',
               'FABP1',
               'MIR584',
               'FADD',
               'MT1E',
               'APOC4',
               'AKR1C2'
)

#location 3

cols_remain<-c('ALDOB',
               'MTCH2',
               'FCER1A',
               'ETFB',
               'HRG',
               'MRPL40',
               'MT1H',
               'PTPN6',
               'DCXR',
               'ANGPTL3',
               'MRPL12',
               'AKR1C4'
)

#overall
cols_remain<-c('TCERG1',
               'FOS',
               'ANGPT1',
               'ADORA2A',
               'LSM5',
               'RABL2A',
               'TRBV7-4',
               'MIR584',
               'C20orf96',
               'GSTM1',
               'NDRG3',
               'WDR61',
               'SUSD3',
               'XRN2',
               'CLK1',
               'MYL1',
               'TMSB4Y',
               'RAMP3',
               'RSPH9','RPL23AP7',
               'TSPAN7',
               'AMBP',
               'SNORD3B-1',
               'SNORD3B-2',
               'C20orf96',
               'COMMD5',
               'RNASE6',
               'SNORD3D',
               'FMO5',
               'CHMP4C',
               'CPB2',
               'LUC7L3',
               'AKR1C1',
               'PLA2G16',
               'ANXA1',
               'S100A13',
               'FABP1',
               'MIR584',
               'FADD',
               'MT1E',
               'APOC4',
               'AKR1C2',
               'ALDOB',
               'MTCH2',
               'FCER1A',
               'ETFB',
               'HRG',
               'MRPL40',
               'MT1H',
               'PTPN6',
               'DCXR',
               'ANGPTL3',
               'MRPL12',
               'AKR1C4'
)

newdata<-data[ ,colnames(data) %in% cols_remain]

ma <- as.matrix(newdata)
row.names(ma) <- row.names(data)
annot_colors=list(location=c(HV="pink",PA="yellow", PV='light blue', PoV = 'green'))

pheatmap(ma, scale = "column", clustering_method = "average",
         annotation_colors = annot_colors,
         annotation_row = data[, 86, drop =FALSE],
         show_rownames = TRUE)


pheatmap(ma, scale = "column", clustering_method = "average",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_colors = annot_colors,
         annotation_row = data[, 13006, drop =FALSE],
         show_rownames = FALSE,
         angle_col = '45', cellwidth = 12)

pheatmap(ma, scale = "column", clustering_method = "average",
         color = colorRampPalette(c("white", "firebrick3"))(100),
         annotation_colors = annot_colors,
         annotation_row = data[, 86, drop =FALSE],
         show_rownames = FALSE,
         angle_col = '45')

pheatmap(ma, scale = "column", clustering_method = "average",
         color = colorRampPalette(c("white", "navy"))(100),
         annotation_colors = annot_colors,
         annotation_row = data[, 86, drop =FALSE],
         show_rownames = FALSE,
         angle_col = '45')


