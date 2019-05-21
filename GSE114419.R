#Install all the packages first
#affy, affyPLM, limma, ArrayExpress, Biobase
#Biostrings,genefilter

#expression for installation

BiocManager::install("Oligo")
library(affy)
###############################################
###############################################
library(oligo)
#Loading data -- Affymatrix
celpath = "F:/May 19/GSE114419/"
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)

#Retrive the annotation of the data
ph = data@phenoData


#Retrive probe annotation 
feat = data@featureData

#Retrive experiment annotation
exp = data@experimentData

###Pause and read the experiment to understand the 
##samples ... to name them appropriately
###FOR GSE114419 data
ph@data[,1] = c("Control1","Control2","Control3","PCOS1","PCOS2","PCOS3")

# Quality COntrol Check
op = par(mfrow = c(2,3))
for(i in 1:6){hist(data[,i],lwd =2,which ='all',ylab='Density',xlab='Log2 Intensities',main=ph@data$index[i] )}


#Quality COntrol Histogram
color = c('green','green','green','red','red','red')
hist(data[1:6],lwd = 2,which='all',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data')


#QC Box plots
name = "boxplot.jpg"
jpeg(name)
boxplot(data,which='all',col='red',names=ph@data$index)
dev.off()

#QC MA plot
for (i in 1:6)
{
  name = paste("MAplot",i,".jpg",sep="")
  jpeg(name)
  MAplot(data,which=i)
  dev.off()
}

#Calculating QC of array
#calculates the quality measure
#data_qc = qc(data)
#avbg(data.qc)
#The above value must be commparable to all the arrays

## NORMALIZATION
#using RMA method
data.rma = rma(data)
data.matrix = exprs(data.rma)

#Or using GCRMA 
library("gcrma")
data.gcrma = gcrma(data)
data.matrix = exprs(data.gcrma)


##Rechecking the Quality
#Boxplot
name = "boxplot_after_gcrma.jpg"
jpeg(name)
boxplot(data.matrix,which='pm',col='red',names=ph@data$sample)
dev.off()

#MA 
for (i in 1:6)
{
  name = paste("MAplot_after_gcrma",i,".jpg",sep="")
  jpeg(name)
  MAplot(data,which=i)
  dev.off()
}

#PCA
#color=c('green','red','red','red','green','green','green','red')
#data.PC = prcomp(t(data.matrix),scale.=TRUE)
#plot(data.PC$x[1:2],col=color)


########### DE ANALYSIS ######## FOR Control Vs Test ###
#                 VIA LIMMA  
ph@data[ ,2] = c("control","control","control","PCOS","PCOS","PCOS")
colnames(ph@data)[2] = "source"
groups = ph@data$source

f = factor(groups,levels = c("control","PCOS"))
#Most general model for two sample analysis
design = model.matrix( ~ 0 + f)
colnames(design) = c("control","PCOS")
# make sure limma library is loaded
data.fit = lmFit(data.matrix,design)
#publishes the top 15 genes
data.fit$coefficients[1:15,]

# allocating the division for comparision
contrast.matrix = makeContrasts(PCOS-control,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)

#Performing t-test
data.fit.eb = eBayes(data.fit.con)

#the log2fold change value
data.fit.eb$coefficients[1:15,]

#the statistical values from t-test
data.fit.eb$t[1:15,]
data.fit.eb$p.value[1:15,]

#Creating the volcano plot
name = "Volcano.jpg"
jpeg(name)
#coef = 1,when its linear , 2 or 3, when comparision between 2 seperate instances, 4 when comparision
#between before and after treatment setups
volcanoplot(data.fit.eb,coef=1L,highlight=10)
dev.off()
#
#
#######DE Genes using TopTable Method
options(digits=2)
tab = topTable(data.fit.eb,coef=1,number=200,adjust.method="BH")
#
topgenes = tab[tab[, "adj.P.Val"] < 0.7, ]
dim(topgenes)
topgenes
