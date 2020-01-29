
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

WD<-args[1] #working directory
RSRC<-args[2] #path to R scripts to source
MAX_COMB<-args[3]
MAXITER<- args[4]
CORES<-args[5]
GENDIST<-args[6]
SAMPLES<-args[7]
FILES<-args[8:length(args)]

#source the R files
file.sources = list.files(path=paste0(RSRC), pattern="*.R$")
sapply(paste0(RSRC, "/", file.sources),source)

print(FILES)
surfaces <- raster::stack(FILES)
samples <- read.table(file = SAMPLES, sep = '\t', header = TRUE)
sample.locales <- SpatialPoints(samples[,c(14,15)])

gen.dist <- read.table(file = GENDIST, sep = '\t', header=T)
gen.dist <- otuSummary::matrixConvert(gen.dist, colname = c("ind1", "ind2", "dist"))
gen.dist <- gen.dist[,3]

#print("ready to prep")

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
                     Plots.dir=".")
saveRDS(GA.inputs, file="GA-inputs.MS.rds")

gdist.inputs <- gdist.prep(length(sample.locales),
                           samples = sample.locales,
                           response = gen.dist,
                           method = 'commuteDistance')
saveRDS(gdist.inputs, file="gdist-inputs.MS.rds")

all_comb_pipeline_make_batches(gdist.inputs,
                          GA.inputs,
                          max.combination = MAX_COMB)
