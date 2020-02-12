list.of.packages <- c("ape","ggtree","phangorn","ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list=ls())


load("treeScriptWorkspace.Rdata")
myfile <- file.choose(new=FALSE)
original_csv_data <- read.csv(myfile,header=TRUE,sep=",")
csv_data <- read.csv(myfile,header=TRUE,sep=",")
rownames(csv_data) <- csv_data[,1]
head(csv_data)

Dock_score <- original_csv_data$Dock_score
rownames(original_csv_data) <- original_csv_data[,1]
rownames(Dock_score) <- rownames(original_csv_data)
Dock_scoredata <- cbind(as.data.frame(row.names(original_csv_data)),as.data.frame(Dock_score))

Aff_score <- csv_data$Aff
pred <- csv_data$DNN
# scatter plot

sp <- ggplot(csv_data, aes(x=Aff_score, y=pred, color=Dock_score, size = -Dock_score)) + geom_point() # , size = 5
sp + scale_color_continuous(low='darkgreen', high='red') +
  theme(legend.position="right")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_jitter(alpha = 0.1) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) + 
  theme(text = element_text(color="black", size=15)) +
  scale_size(range = c(1,6))

