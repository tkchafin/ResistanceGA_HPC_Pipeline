
# Optimizing Resistance Surface using ResistanceGA
# Need packages: ResistanceGA , rgdal, otuSummary
# Setting up the workspace
library(ResistanceGA)
library(otuSummary)
library(raster)
library(tidyverse)
library(GA)
library(MuMIn)
library(lme4)
library(gdistance)
library(ggplot2)
library(ggExtra)
library(Matrix)
library(utils)
library(akima)
library(plyr)
library(dplyr)
library(magrittr)
library(spatstat)
library(grDevices)

args <- commandArgs(trailingOnly=TRUE)

#if (length(args) < 6){
#  print("Wrong number of inputs")
#  quit()
#}

#WD<-args[1] #working directory
#all.comb<-args[1] #output directory name (e.g. all_comb/Run1/rep_1)
Plotdir<-args[1]
CORES<-as.numeric(args[2]) #number of cores per process
RSRC<-args[3] #path to R scripts to source
GENDIST<-args[4]
SAMPLES<-args[5]
MAXITER<-as.numeric(args[6])
FILES<-args[7:length(args)] #rasters to operate on

#WD<-"/Users/tkchafin/Downloads/WTD_RGA"
#all.comb<-"all_comb/rep_1"
#CORES<-4 
#RSRC<- "R"
#FILES<-c("rasters/MZ_1_MajorHighways.reduce.proj.tif", "rasters/MZ_2_RiverStreamOrder.reduce.proj.tif")

print(Plotdir)
print(FILES)

#setwd(WD)
#print(WD)

if(!file.exists(Plotdir)){
  dir.create(Plotdir, recursive=T)
}

#source the R files
file.sources = list.files(path=paste0(RSRC), pattern="*.R$")
sapply(paste0(RSRC, "/", file.sources),source)

# loading the data
#files <- list.files(path=".", pattern="*reduce.proj.tif$", full.names=TRUE, recursive=FALSE)

surfaces <- raster::stack(FILES)
samples <- read.table(file = SAMPLES, sep = '\t', header = TRUE)
sample.locales <- SpatialPoints(samples[,c(14,15)])

gen.dist <- read.table(file = GENDIST, sep = '\t', header = TRUE)
gen.dist <- otuSummary::matrixConvert(gen.dist, colname = c("ind1", "ind2", "dist"))
gen.dist <- gen.dist[,3]

print("ready to prep")

# Running the Resistance Surface Optimization
# Format inputs
GA.inputs <- GA.prep(ASCII.dir = surfaces,
                     Results.dir = ".",
                     method = 'LL',
                     max.cat = 500,
                     max.cont = 500,
                     maxiter=MAXITER,
                     seed = NULL,
                     parallel = CORES,
                     quiet = FALSE,
                     # this value should be justified
                     k.value = 3, 
                     Plots.dir=Plotdir)
saveRDS(GA.inputs, file=paste0(Plotdir, "/", "GA-inputs.rds"))

gdist.inputs <- gdist.prep(length(sample.locales),
                           samples = sample.locales,
                           response = gen.dist,
                           method = 'commuteDistance')
saveRDS(gdist.inputs, file=paste0(Plotdir, "/", "gdist-inputs.rds"))

# Runs single surface opitimization for each surface independently followed by optimization and
# model selection of combination surfaces
#need to tell process which daughter process it is (index) and how many total processes there are
#1-based index on daughter processes, BTW
comb.analysis <- all_comb_pipeline_SS(gdist.inputs,
                          GA.inputs,
                          ".",
                          max.combination = length(surfaces),
                          iters = 100,
                          replicate = 1,
                          sample.prop = 0.75,
                          nlm = FALSE,
                          dist_mod = TRUE,
                          null_mod = TRUE, 
                          Plots.dir = Plotdir)
