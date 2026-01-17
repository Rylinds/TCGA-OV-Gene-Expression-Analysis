#### Section 1: Set options and environment ####

# A reminder about comments in this pipeline
library(ggplot2)
library(ggrepel)

# Create an empty list object to store options inside of. Keeps things organized
opt <- list()

# Describe the path used as the working directory
# opt$wd <- "~//project/TCGA-OV"

# Create a version of the R version that is more path-friendly by removing special characters
opt$rString <- gsub(" ","_",gsub(" \\(.*","",(R.version.string)))

# Use the R version to create a path to store the packages (this helps avoid permission issues on various computers)
# opt$yourPckDir   <- file.path("./Rpackages",opt$rString)

# Describes a path to a package directory on the Edupod to speed up setup when using the Edupod
# opt$sharedPckDir <- "hidden Edupod path"

# Specify packages to be installed from the CRAN server (the typical repository of R packages)
opt$cranPackages <- c("BiocManager","ggplot2","ggrepel","remotes")

# Specify packages to be installed from the Bioconductor server (a biology specialized repository of packages)
opt$biocPackages <- c("DESeq2","Biostrings","SummarizedExperiment","BioinformaticsFMRP/TCGAbiolinksGUI.data","BioinformaticsFMRP/TCGAbiolinks")

# Watch for pop-up box(es) and select "no" if it asks to install from sources needing compilation.
# It saves a copy of itself to your opt$wd directory if you want to read it.

# Downloads and runs a shared script file to configure your computer.
source("public redaction")

# If this doesn't work, start over from the beginning and watch carefully for error codes.

# Check to see if the trickiest to install package is able to be loaded.
library(TCGAbiolinks)

# Due to an error in the release version of TCGAbiolinks,
# you cannot use the typical method of installing TCGAbiolinks.

# Ensures an installed package version avoids a known release error
if(packageVersion("TCGAbiolinks")<"2.37.1"){ #This comparison actually works due to class of packageVersion() output.
    stop("You need to load a more recent version of TCGAbiolinks. This is implemented in my configuration file.\n\nTry deleting the existing TCGAbiolinks installations, fully restarting RStudio (don't save namespace) and repeating the configuration code.")
}



#### Section 2: Load sample sheets and check quality ####
# The group data.frames result from searching the GDC and exporting "Sample sheets"

# Structure sample sheets files into data.frame objects
group1 <- read.csv("TCGA-OV_Group1.csv", fill=T, header = 1)
group2 <- read.csv("TCGA-OV_Group2.csv", fill=T, header = 1)

# The purpose of this section is to check for certain issues with sample files.

# Check that group1 and group2 were structured correctly
if( !"data.frame"%in%class(group1) ){ stop("Oh, no! Your group1 object is not a data.frame. Go fix this!") }
if( !"data.frame"%in%class(group2) ){ stop("Oh, no! Your group2 object is not a data.frame. Go fix this!") }

# Check column names in the sample sheet needed in the code later on.
opt$checkedColumns <- c("Project.ID","Sample.ID","Data.Type")
if( !all(opt$checkedColumns%in%colnames(group1)) ){ stop("Oh, no! You are missing important column names. Check your input file / code for group 1.") }
if( !all(opt$checkedColumns%in%colnames(group2)) ){ stop("Oh, no! You are missing important column names. Check your input file / code for group 2.") }

# This is needed when there are multiple files describing the same sample
#    e.g. if a sequencing run was bad, they might have sequenced the same sample twice.
# Ideally fix this by deleting entries from sample sheets.

# Identify repeated entries that could distort sample comparisons
if( length(unique(group1$Sample.ID)) < nrow(group1) ){stop("Oh, no! You have some duplicate values in the Sample.ID column of group 1; investigate and correct this.")}
if( length(unique(group2$Sample.ID)) < nrow(group2) ){stop("Oh, no! You have some duplicate values in the Sample.ID column of group 2; investigate and correct this.")}
if( length(unique(group1$Case.ID  )) < nrow(group1) ){stop("Oh, no! You have some duplicate values in the Case.ID column of group 1. This is nearly always in error; investigate and correct this if needed.")}
if( length(unique(group2$Case.ID  )) < nrow(group2) ){stop("Oh, no! You have some duplicate values in the Case.ID column of group 2. This is nearly always in error; investigate and correct this if needed.")}

# Prevent excessive memory use by limiting total sample count
if((nrow(group1)+nrow(group2))>250 ){
    stop("You have too many samples and might run into computer issues. ",
         "Strongly consider subsetting your sample sheet!!!")
}else if(max(c(nrow(group1),nrow(group2)))>100){
    stop("You have a large number of samples in one of your groups. ",
         "Consider randomly subsetting your sample sheets to to reduce their size ",
         "and avoid calculation issues.")
}

# Produce an error if the number of cases is too small and likely of limited statistical reliability
if(min(c(nrow(group1),nrow(group2)))<6){
    stop("You have too few samples in one of your groups; check your inputs. ")
}

# Check that the data.type needed for differential expression analysis is described
if(!all(c(group1$Data.Type,group2$Data.Type)=="Gene Expression Quantification")){
    stop("This pipeline is only built to use Gene Expression Quantification data. ")
}

#### Section 3: Search for data to download ####
##### Search GDC for all files associated with target project #####

# If this doesn't work, start over from the beginning and watch carefully for error codes.

# Load the primary library used to search the GDC dataset
library(TCGAbiolinks)

# Create a character vector that determines what projects to search through on the GDC based on the grouping files
opt$project <- unique(c(group1$Project.ID,group2$Project.ID))

# Do a broad search of the GDC that returns all relevant files
query1 <- GDCquery(
    project = opt$project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    access = "open"
)

##### Investigate query results to identify targeted samples ####
# Extract the data.frame from inside the initial
samDf <- query1$results[[1]]

# Determine which search results correspond to the samples in the group files
samDf$inGrp1      <- samDf$sample.submitter_id%in%group1$Sample.ID
samDf$inGrp2      <- samDf$sample.submitter_id%in%group2$Sample.ID
samDf$sampleLogic <- samDf$inGrp1 | samDf$inGrp2

# Edit your sample files to remove problematic samples if needed.

# Visualize and explore your search results to double check it's working.
View(samDf[order(-samDf$sampleLogic),])

# Save the desired sample barcodes (a.k.a. cases) to use to filter another query of the database
opt$initBarcodes <- samDf$sample.submitter_id[samDf$sampleLogic]

##### Search GDC again but limit results to targeted samples ####
# Search the GDC and filter the results by the file barcodes
query2 <- GDCquery(
    project = opt$project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    access = "open",
    barcode = opt$initBarcodes #Note that this is the change from the first search.
)


##### Check query results #####

# Check to make sure your grouping files have data
length(group1$Sample.ID) #Should be >0 ---- TCGA-OV has 17 - pass
length(group2$Sample.ID) #Should be >0 ---- TCGA-OV has 10 - pass

# Check that the grouping files do not overlap:
sum(group1$Sample.ID%in%group2$Sample.ID) #Should be 0 ---- pass
sum(group2$Sample.ID%in%group1$Sample.ID) #Should be 0 ---- pass

# Check that the second search had some results
samDf2 <- query2$results[[1]]
nrow(samDf2) #Should be >0 ---- 27, pass

# These numbers are nebulous, but I usually suggest no larger than 150 in a group to avoid computation issues.
# Try randomly deleting rows from a sample sheet to make it smaller if needed.

# Check that your second search had results that overlapped with BOTH group files
samDf2$inGroup1 <- samDf2$sample.submitter_id%in%group1$Sample.ID
samDf2$inGroup2 <- samDf2$sample.submitter_id%in%group2$Sample.ID
sum(samDf2$inGroup1) #Should be >=6 and not >>>50 ---- 17, pass
sum(samDf2$inGroup2) #Should be >=6 and not >>>50 ---- 10, pass

## Q.32
# Check to make sure your search results did not have extra rows.
all(samDf2$inGroup1|samDf2$inGroup2) # ---- TRUE


#### Section 4: Obtain and structure the sequencing data into an R object ####
##### Download and prepare data into an object ####
# This may require retrying if it fails.
# Don't download >3 Gb on Edupod server.

# # Download the data files from the GDC
GDCdownload(query = query2, method = "api", files.per.chunk = 10)

# This is the step effected by the error in TCGAbiolinks
# Therefore I have added logic to help avoid issues.
if(packageVersion("TCGAbiolinks")>="2.37.1"){
    dds1 <- GDCprepare(query = query2)
}

##### Check the quality of the prepared data ####
# This information is stored inside of a complex object that requires special functions to pull out.
# Once the data is pulled out, I have converted it to a data.frame so that it could be viewed.

# Visualize the clinical data describing the samples to double check the preparation worked.
library(SummarizedExperiment) #Needed to access the function colData()
View(as.data.frame(colData(dds1)))

# See which sample IDs are duplicated (if any)
which(duplicated(dds1$sample_submitter_id)) #Should be none; remove with: dds1 <- dds1[,!duplicated(dds1$sample_submitter_id)]
# got integer(0) ---- pass

# Load packages for modifying the sequencing data object
library(SummarizedExperiment)
library(DESeq2)

# These are multi-dimensional data set with multiple sections:
#   - colData   #This section houses a table describing each sample, accessed via $ or colData(dds)
#   - rowRanges #This section houses a table describing each locus (typically genes), accessed via rowRanges(dds)
#   - assays    #This is a list of data.frames describing the counts per locus (rows) and sample (columns) combination, the raw count data can be accessed via assays(dds)[[1]].
#
# The extraction tools change behavior for this data:
#   - $ and [[]] are used to access / modify / and create sample descriptors. I recommend $.
#   - [,] is used to subset the loci (using the rows argument) or the samples (using the column argument)
#
# Documentation: https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html

# Add new columns to the clinical data for use in downstream data cleaning
dds1$group1 <- dds1$sample_submitter_id%in%group1$Sample.ID
dds1$group2 <- dds1$sample_submitter_id%in%group2$Sample.ID

# Add a factor column to the clinical data to control how DESeq2 will group the samples
dds1$comp   <- factor(x = dds1$group2,levels = c(FALSE,TRUE))

# This is also where you define what data DESeq2 should use to control the grouping.

# Convert the sequencing data into DESeqDataSet class (a subclass of RangedSummarizedExperiment)
dds1 <- DESeqDataSet(dds1,design = ~comp)

# Visualize your samples' clinical data and ensure they are what you expected
View(as.data.frame(colData(dds1)))




#### Section 5: Analyze count data #####
##### Detect bad samples based on PCA of normalized counts #####

# Create a function to visualize overall clustering of the samples via PCA
pcaFun <- function(ddsFun,plotToggle=T){
    ### Check arguments of pcaFun prior to analysis
    if(!"DESeq2"%in%installed.packages()){
        stop("DESeq2 needs to be available. Run the configuration file.")
    }
    if(class(ddsFun)=="DESeq2"){
      stop("ddsFun needs to be a DESeq2 class object")
    }
    if(!"comp"%in%colnames(colData(ddsFun))){
      stop("The column 'comp' needs to be the name of a factor column in ddsFun's sample descriptor section")
    }

    ### Calculate normalized values
    ddsFun <- DESeq2::estimateSizeFactors(ddsFun)
    normCounts <- as.data.frame(counts(ddsFun, normalized = T))

    ### Calculate the variance of counts across each gene in order to subset data
    # This is needed because only rows with > 0 variance can be used in PCA
    # The apply loops function var() (variance) across rows of normCounts
    perLocusVar   <- apply(normCounts,MARGIN = 1,var)

    ### Determine which genes to use based on variance
    pcaLocusLogic <- perLocusVar>median(perLocusVar[perLocusVar>0])
    if(sum(pcaLocusLogic)==0){stop("Your data did not contain loci with variable counts")}

    ### Subsets and transposes the dataset for analysis with prcomp
    pcaInput <- t(normCounts[which(pcaLocusLogic),])

    ### Calculate principal components with scaling (to make each locus equal in weight)
    pca <- prcomp(pcaInput,scale. = T)

    ### Calculate the percent variance explained by each principal component
    pca_var <- pca$sdev^2
    pve <- pca_var / sum(pca_var)

    ### Plots a scatterplot of the samples visualized via principal component analysis
    if(plotToggle){
        plot(
            x = pca$x[,1],y = pca$x[,2],asp = 1,
            xlab = paste0("PC1: ",round(pve[1],3)*100,"%"),
            ylab = paste0("PC2: ",round(pve[2],3)*100,"%"),
            col  = colData(ddsFun)$comp,
            pch  = as.numeric(colData(ddsFun)$comp)+1
        )
    }

    ### Return a list of the data
    outLs <- list(x=pca$x,pve=pve)
    return(outLs)
}

# Calculate and visualize how the samples cluster via PCA
pca1Out <- pcaFun(ddsFun = dds1)

# Analyze the PCA results to determine undesirable samples / barcodes
opt$pca1Cutoffs <- 100

abline(v = opt$pca1Cutoffs, col = "blue", lty = 2)   # right cutoff
#abline(v = -opt$pca1Cutoffs, col = "blue", lty = 2)  --- optional left cutoff, in case of left outliers

#abline(v = opt$pca1Cutoffs,col="blue")
dds1$pca1Logic  <- pca1Out$x[,1]>opt$pca1Cutoffs
badSamples      <- dds1$barcode[dds1$pca1Logic]

# Subset the data to certain samples based on the SummarizedExperiment documentation
dds2 <- dds1[,!dds1$barcode%in%badSamples]
dim(dds2)

# After initial sample filtering, repeat PCA to see if additional filtering is needed
pca2Out <- pcaFun(dds2)

##### Detect bad loci based having high frequency of 0 in either grouping#####

# Calculate proportion with zero raw counts per grouping
locusCounts <- as.data.frame(counts(dds2, normalized = F))
locusCounts_grp1At0 <- rowMeans(locusCounts[,dds2$group1]==0)
locusCounts_grp2At0 <- rowMeans(locusCounts[,dds2$group2]==0)


# The cutoff of 0.9 is somewhat arbitrary.
# The main problem is that loci with 0 in one group become Infinitely significant
# However, loci that are simply extremely low in one group also cause a problem!
# Such loci require modified methods to deal with.

# Visualize the proportion of the samples with 0 read depth per grouping
hist(x = c(locusCounts_grp1At0,locusCounts_grp2At0),100)
opt$maxPct0Cutoff  <- 0.9
abline(v = opt$maxPct0Cutoff,col="red")

# Determine which loci meet the read depth based filter across both groups
locusCounts_logic <-
    locusCounts_grp1At0 < opt$maxPct0Cutoff &
    locusCounts_grp2At0 < opt$maxPct0Cutoff

# Subset the data to certain loci based on the SummarizedExperiment documentation
dds <- dds2[locusCounts_logic,]
dim(dds)



##### Check subset count data #####

# Calculate your total gene and sample count
nrow(dds) #The number of genes should be >>> 1,000. TCGA-OV was above 38 thousand ---- 46,295, pass
ncol(dds) #The number of samples should be at least 12, probably less than ~100, and definitely less than 200. ---- 26, pass

# Calculate the number of samples per grouping
class(dds$comp) #This needs to be a factor type object
table(dds$comp) #You should have two values > 6. ---- 17, 9, pass

# Calculate the minimum sample raw depth
min(colMeans(assays(dds)[[1]])) #Needs to be above 0 ---- 472.0039, pass

# Calculate the maximum proportion of samples with 0 counts per group
max(rowMeans(assays(dds)[[1]][,dds$group1]==0)) #Needs to be less than 100% ---- 0.8823529, pass
max(rowMeans(assays(dds)[[1]][,dds$group2]==0)) #Needs to be less than 100% ---- 0.8888889, pass

#### Section 6: Analyze the data with DESeq2 #####
##### Do the differential expression analysis ####
# Be aware that we are only doing the most basic of analyses with DESeq2.
# It's a big topic: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Repeat the conversion to DESeqDataSet class object to make sure it's still well formatted
dds <- DESeqDataSet(dds,design = ~comp)

# Recalculate the size factors used in depth normalization after all the previous filtering
dds <- estimateSizeFactors(dds)

# Run the actual differential expression calculations.
dds <- DESeq(dds)

# Organize the results of the differential expression calculations into a table
res <- results(dds)

# Combine the results of differential expression with descriptions of the loci
resOutput <- cbind(as.data.frame(res),as.data.frame(rowRanges(dds)))

# Remove a duplicated column that causes problems downstream otherwise
resOutput <- resOutput[,!duplicated(colnames(resOutput))]
#write.csv(resOutput, "~/Desktop/Loci.csv", row.names=TRUE)

##### Analyze the results ####

# Plot a simple scatterplot / volcano plot of the analyzed genes
plot(resOutput$log2FoldChange,-log10(resOutput$padj))

###### End of Pipeline ########

source("PipelineHelpers.R")

resSummary(resOutput)

nicePCA(pca1Out,dds1$group2)+   # My first PCA plot before outlier removal
    geom_vline(xintercept = opt$pca1Cutoffs,color="red") #Showing the cuttoff I used for filtering samples based on pca1
nicePCA(pca2Out,dds2$group2)    # My second PCA plot after outlier removal
niceVolcano(resOutput,0.05,1.5) # My volcano plot

# Getting the supplementary file Samples.csv
#sample_info <- as.data.frame(colData(dds1))[
#  , c("barcode",
#      "group1", "group2",  # grouping variables
#      "patient", "sample", "shortLetterCode", "definition",
#      "tumor_descriptor", "sample_type", "figo_stage", "tissue_or_organ_of_origin", 
#      "primary_diagnosis", "age_at_diagnosis",
#      "race", "vital_status")
#]

# Create Samples.csv from the FINAL dataset (dds), after outlier removal
#sample_info <- as.data.frame(colData(dds))[ ,
#                                            c("barcode",
#                                              "group1", "group2",      # grouping variables
#                                              "patient", "sample",     # useful clinical columns
#                                              "shortLetterCode",
#                                              "definition",
#                                              "tumor_descriptor"
#                                            )
#]
#nrow(sample_info)
# Write to CSV
#write.csv(sample_info, file = "Samples.csv", row.names = FALSE)

