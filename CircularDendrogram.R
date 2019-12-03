# Install these packages if not installed already.

list.of.packages <- c("ape","ggtree","phangorn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list=ls())


load("treeScriptWorkspace.Rdata")
myfile <- file.choose(new=FALSE)  # It will prompt you to choose a csv file.
csv_data <- read.csv(myfile,header=TRUE,sep=",")

# You may need to adjust the following lines depending on the format of your csv file.
rownames(csv_data) <- csv_data[,1] # set first column to row names 

# Setting up to color the nodes based on a certain feature.
Aff <- csv_data$Aff
rownames(Aff) <- rownames(csv_data)
Affdata <- cbind(as.data.frame(row.names(csv_data)),as.data.frame(Aff))

# Subset the descriptors/fingerprints and assign it to csv_data.
csv_data <- csv_data[,-ncol(csv_data)]        # drop last column. subset your data (don't include endpoints like pIC50, etc.)
head(csv_data)                                # To see what the csv head looks like


#make sure you normalize columns. (for continuous data, not binary)
matTrans <- data.matrix(csv_data)
matTrans <- scale(matTrans)
head(matTrans)

d <- dist(matTrans)  # Default is Euclidean
# d <- distance(matTrans, method = 'tanimoto')

tupgma <- upgma(d, method="ward.D2")  # apply clustering to your data with ward linkage
# pKidata <- cbind(as.data.frame(row.names(csv_data)),as.data.frame(pki))


# Adjust the parameters accordingly.
t4 <- ggtree(tupgma, layout="circular", size=1)
t4 <- t4 %<+% Affdata+ geom_text(aes(color = Aff, label=label, angle=angle, fontface="bold"), hjust=-0.05, size=1) +
  geom_tippoint(aes(color=Aff),alpha=0.5, size=2.5)+
  scale_color_continuous(low='red', high='darkgreen') +
  theme(legend.position="right",text = element_text(size = 15)) +
  geom_treescale(x=15,y=15,width = 100, offset = NULL,
                 color = "white", linesize = 1E-100, fontsize = 1E-100)

print(t4)
# dev.off()

