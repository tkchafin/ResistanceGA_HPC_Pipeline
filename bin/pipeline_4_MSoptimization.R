
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
GD<-args[4]
GA<-args[5]
MY_JOB<-args[6]

#source the R files
file.sources = list.files(path=paste0(RSRC), pattern="*.R$")
sapply(paste0(RSRC, "/", file.sources),source)

print(GD)
print(GA)
exit()
gdist.inputs <- readRDS(GD)
GA.inputs <- readRDS(GA)

print(paste0("MAX_COMB:", MAX_COMB))
print(paste0("MY_JOB:", MY_JOB))

all_comb_pipeline_MS(gdist.inputs,
                          GA.inputs,
                          max.combination = MAX_COMB,
                          my_job_index=MY_JOB)
