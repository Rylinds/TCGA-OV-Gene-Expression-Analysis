## Setup the packages
# I have suppressed messages here for clarity, but
# consider removing it in your own script.
# If viewing the messages doesn't elucidate the issue, contact Nolan.
suppressMessages(source(
  "https://utexas.box.com/shared/static/ftih4imqaj4d7tj7mssrhufckwexn63m.r"
))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(DESeq2))
#### Download and prepare demo data ####
# This file was made using save()
download.file(
  url = "https://utexas.box.com/shared/static/orm7anj6iewp6g6saa5h89qv4gro92k1.zip",
  destfile = "demoObjects.zip"
)
unzip("demoObjects.zip")


## Load the demo data
# See setup demo objects section to get this file
load("demoObjects.Rdata")
## Source the function file
source("wk11Hw_functions.R")
## Show the environment after loading data
print(ls())




## Test nicePca on pca1Out_demo
outA <- nicePCA(pca1Out_demo,dds1_demo$group2)
outA + geom_vline(xintercept = 120 ,color="red")

outB <- nicePCA(pca2Out_demo,dds2_demo$group2)
outB



## Test the niceVolcano function
# Using 0.01 & 2.5/-2.5 as the padj log2fc cutoffs
volc1Out <- niceVolcano(resOutput_demo,0.01,2.5)
#Ensure that the data returned is ggplot class
"ggplot"%in%class(volc1Out)
volc1Out

volc2Out <- niceVolcano(resOutput_demo_subset,0.05)
"ggplot"%in%class(volc2Out)
volc2Out



# Test resSummary function
resSummary(resOutput_demo)

resSummary(resOutput_demo_subset)


