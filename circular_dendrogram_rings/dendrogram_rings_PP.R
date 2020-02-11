list.of.packages <- c("ape","ggtree","phangorn","ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list=ls())


load("treeScriptWorkspace.Rdata")
myfile <- file.choose(new=FALSE)
csv_data <- read.csv(myfile,header=TRUE,sep=",")
rownames(csv_data) = csv_data[,1]
head(csv_data)

# Convert csv_data frame to numeric. Some columns are factors.
csv_data[] <- lapply(csv_data, function(x) as.numeric(as.character(x)))

# pKidata <- cbind(as.data.frame(row.names(csv_data)),as.data.frame(pki))
pAffinity <- csv_data$pAffinity
pAffinitydata <- cbind(as.data.frame(row.names(csv_data)),as.data.frame(pAffinity))

# select variables ""pRandom1",	"pRandom2"
to_color <- c("pRandom1",	"pRandom2")
outer_variables <- csv_data[to_color]

# exclude variables "ppAffinityinity",	"pRandom1",	"pRandom2" from data to cluster
ecfp6 <- names(csv_data) %in% c("pAffinity", "pRandom1",	"pRandom2")
ecfp6_data <- csv_data[!ecfp6]

#make sure you normalize columns! (for continuous data, not binary)
# matTrans <- scale(ecfp6_data)
# head(sortMat)

d <- dist(ecfp6_data)
tupgma <- upgma(d, method="ward.D2")


t4 <- ggtree(tupgma, layout="circular", size=1)
t4 <- t4 %<+% pAffinitydata+ geom_text(aes(color = pAffinity, label=label, angle=angle, fontface="bold"), hjust=-0.05, size=1) +
  geom_tippoint(aes(color=pAffinity),alpha=0.5, size=2.5)+
  scale_colour_gradient(low='red', high='darkgreen', na.value = "grey50") +
  theme(legend.position="right",text = element_text(size = 15)) +
  geom_treescale(x=15,y=15,width = 15, offset = NULL,
                 color = "white", linesize = 1E-100, fontsize = 1E-100)

print(t4)

t3<- gheatmap(t4,outer_variables,width=0.25, colnames_angle=-80,hjust='left', 
              low = "red",
              high = "darkgreen",
              color="grey50",
              font.size=3.2,offset = 11.5)+ 
  theme(legend.position="left")
print(t3)

open_tree(t3,11.8) %>% rotate_tree(11.8)

