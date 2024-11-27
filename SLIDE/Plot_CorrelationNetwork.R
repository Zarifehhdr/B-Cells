################################################
##  This code plots correlation networks ## 
################################################
library(ggplot2) 
library(dplyr) 
library(tidyr) 
library(stringr)
library(qgraph)
library(scales)  

#load X which you used for SLDIE
xdf <- read.csv(".../Data/Concat/scRNA_Male_Donor/X_GC_PB.csv",row.names=1)  
# Load the LF features you want to plot the correlation networks of them: the txt file is in the out folder of SLIDE resluts
data <- read.table(".../feature_list_Z5.txt", header = TRUE, stringsAsFactors = FALSE)
lis = xdf[,c(data[1:dim(data)[1],1])]

# Create a vector to specify the shape for each gene
# "ellipse" for blue genes (Y=0) and "rectangle" for red genes (Y=1)
shape_vector <- ifelse(data$names %in% blue_genes, "ellipse", "rectangle")

# correlation matrix
corr_matrix <- cor(lis)

blue_genes <- data$names[data$color == "Blue"] # High in state corresponding with Y=0
red_genes <- data$names[data$color != "Blue"] # High in state corresponding with Y=0

# Separate colors for circles and squares
circle_color <- "#7BDE7B"  # Change this to the desired color for ellipse
square_color <- "#B83636"  # Change this to the desired color for rectangle

# Create a vector to specify the color for each gene
color_vector <- ifelse(shape_vector == "ellipse", circle_color, square_color)

# path to save the network plot 
pdf(".../Z5_corrnet.pdf", height = 8, width = 8)

qgraph(
  corr_matrix,
  layout = "circular",
  colors = color_vector,
  threshold = 0.4, # You can change this threshold value
  repulsion = 0.4, 
  labels = colnames(corr_matrix),
  label.font = 2,
  label.scale.equal = TRUE,
  label.prop = 0.95,
  shape = shape_vector,  
  posCol = "springgreen4",  # Positive correlations
  negCol = "orange",        # Negative correlations
  height = 5,
  width = 5,
  label.cex = 1.1
)

legend("topright", legend = c("Positive Correlation", "Negative Correlation"), title = "Edge Colors", col = c("springgreen4", "orange"), lty = 1, lwd = 3.5, cex = 1,bty = "n")
legend(
  "topleft",
  legend = c("High in GC", "High in PB"), # change the legend for other states
  title = "Genes",
  pch = c(21, 22),
  col = "black",  # Border color of shapes
  pt.bg = c(circle_color, square_color),  # Specify fill colors for shapes in the legend
  pt.cex = 2,
  cex = 1,
  bty = "n"  # No box around the legend
)
dev.off()




