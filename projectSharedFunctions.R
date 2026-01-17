#Link to this file: https://utexas.box.com/shared/static/2q9p7yc3ixieypl8r7rbzex998nisy86.r
#Note to future Nolan: since this file downloads itself, be sure to update the link as needed.
message("Running version 20251027.06 of project setup script...")
#Run this code every time the file is sourced
if(exists("opt")){
    opt$config <- NULL
    opt$config$projectConfig <- function(opt){
        ## Check input
        if(class(opt)!="list"){
            stop("Object opt need to a list class object")
        }

        ## Check opt object column names
        reqNames <- c("wd","yourPckDir","cranPackages","biocPackages")
        reqNamesDesc<-c(
            "A length 1 character vector describing path to the working directory. Typically this is '~/inb321g/project'.",
            "A length >=1 character vector describing path to the directory it should install packages. Typically this is '~/inb321g/project/Rpackages'.",
            "A length >=1 character vector describing the packages to be installed from the CRAN repository. In otherwords, these are the packages installed using `install.packages()`.",
            "A length >=1 character vector describing the packages to be installed from the Bioconductor suite of packages. This is a collection of packages commonly used by biologists. It uses an install tool named `BiocManager::install()` to install. It is also capable of installing from github if you describe the package using the 'repoName/packageName' format."
        )
        if(!all(reqNames%in%names(opt))||
           !all(unlist(lapply(opt[reqNames],class))=="character")||
           !all(unlist(lapply(opt[reqNames],length))>0)){
            message("Object opt must have the following container names:\n",
                    paste0(reqNames,": ",reqNamesDesc,collapse = "\n\n"))
            stop("Exiting the setup unsuccessfully due to a malformed `opt` object. \nSee the above output for feedback and try again!")
        }
        message("\nChecked your opt object container names.")

        ##Set working directories
        #Setup Working directory first because packages might go in it
        dir.create(opt$wd,recursive = T,showWarnings = F)
        if(!dir.exists(opt$wd)){stop("Could not create working directory. Tested for: ",opt$wd)}
        setwd(opt$wd)
        message("\nWorking directory set to: ",getwd())

        ##Store a copy of this file in the working directory
        download.file(
            destfile = "projectSharedFunctions.R",
            url = "https://utexas.box.com/shared/static/2q9p7yc3ixieypl8r7rbzex998nisy86.r", #Nolan: Don't forget to update this if it changes!
            quiet = T,verbose=F
        )
        message("\nDownloaded a copy of this script to the wd for your records.")

        ##Make the shared directory nothing if it doesn't exist
        message("\nDetermining if the shared package directory path should be used:")
        if(!all(dir.exists(opt$sharedPckDir))){
            sharedPckDir <- NULL
            message("The shared package directory path\n '",opt$sharedPckDir,"'\n did not refer to an existing directory.\n This will therefore be ignored. \nThis is expected on your personal computer, but\n it is unusual on the EduPod.")
        }else{
            sharedPckDir <- opt$sharedPckDir
            message("R will attempt to use the shared packages directory \n(",opt$sharedPckDir,")\n when loading packages.")
        }

        ## Make the package directory container
        opt$pckDir <- unique(c(sharedPckDir,opt$yourPckDir))

        ## Change the package directory if you are Nolan (this code is here to help me keep the shared directory updated)
        if(path.expand("~/")=="/stor/home/nbb624/"){opt$pckDir <- opt$sharedPckDir}

        ##Setup package directory
        dir.create(opt$pckDir[1],recursive = T,showWarnings = F)
        if(!dir.exists(opt$pckDir[1])){
            stop("Could not create package directory. Tested for: ",opt$pckDir[1])
        }
        .libPaths(file.path(opt$pckDir)) #This function lets you both set and return the library path
        message("\nPath to where R should look for packages (the 'library path') set to: \n",
                paste0(.libPaths(),collapse = "\n"))

        #### Install and/or load packages ####
        ##To allow for github sourced packages, produce a repo stripped version of package vectors
        cleanedCran<-gsub(".*\\/","",opt$cranPackages)
        cleanedBioc<-gsub(".*\\/","",opt$biocPackages)

        ##Install CRAN derived packages
        if(!all(cleanedCran%in%installed.packages())){
            message("\n\nAbout to start installing CRAN packages.\nSelect 'no' to popup messages about whether or not you should install from sources requiring compilation.")
            install.packages(opt$cranPackages,lib = opt$pckDir[1],quiet = T)
            message("\nCRAN packages installed using `install.packages()`.")
        }

        ##Install Bioconductor derived packages
        if(!all(cleanedBioc%in%installed.packages())){
            message("\n\nAbout to start installing Bioconductor packages.\nSelect 'no' to popup messages about whether or not you should install from sources requiring compilation.")
            BiocManager::install(opt$biocPackages,update = T,force = T,ask = F,lib = opt$pckDir[1],quiet = T)
            message("\nBioconductor packages installed using `BiocManager::install()`.")
        }

        #### If using personal and shared class directory, don't allow key packages
        # To be implemented

        ##Load packages
        message("\nChecking for specific package installations...")
        if(!requireNamespace("BiocManager" ,quietly = T)){stop("BiocManager" ," could not be loaded! Try again and/or contact an instructor.")}
        if(!requireNamespace("DESeq2"      ,quietly = T)){stop("DESeq2"      ," could not be loaded! Try again and/or contact an instructor.")}
        if(!requireNamespace("TCGAbiolinks",quietly = T)){stop("TCGAbiolinks"," could not be loaded! Try again and/or contact an instructor.")}
        message("BiocManager, DESeq2, and TCGAbiolinks installations seem valid.")

        ##Add messages describing outcome
        opt$configMessage <- c(opt$configMessage,
                               paste0("Set the working directory to:\n  ",getwd())
        )
        opt$configMessage <- c(opt$configMessage,paste0(
            "\nSet R to look for packages in the following directories:\n  ",
            paste0(.libPaths(),collapse = "\n  ")
        ))

        ## Message success if successful
        message("\nConfiguration successful!!!!")
        return(opt$configMessage)
    }
    opt$config$configResult <- try(opt$config$projectConfig(opt = opt))
    if(class(opt$configResult)[1]=="try-error"){
        stop("Something went wrong! Here are some general troubleshooting tips:\n",
             " - You should see an error message above. Try to follow it first if you can.\n",
             " - If the issue is loading packages, quit the R session without saving the workspace.\n",
             " - You may need to do this a few times.\n",
             " - Each time you retry, you will need to recreate opt if you are doing it right.\n",
             " - Contact Nolan if this message keeps appearing, as this script may require an update for your system.\n\n",
             "Here is a print out of the error message this script 'caught':\n",
             opt$configResult
        )
    }
}else{
    stop("\nYou need to create a list class object named `opt`.\n\n",
         "`opt` should have the following containers and values:\n",
         "    wd = The path to the working directory\n",
         "    yourPckDir = The path to the package directory\n",
         "    cranPackages = The packages from the CRAN to be installed\n",
         "    biocPackages = The packages from bioconductor to be installed\n",
         "\nOnce you have made `opt`, rerun this script file via the source function.\n")
}
