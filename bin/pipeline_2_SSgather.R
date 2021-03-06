
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
#all.comb<-args[2] #output directory name (e.g. all_comb/Run1/rep_1)
WD<-args[1]
RSRC<-args[2] #path to R scripts to source
#FILES<-args[3:length(args)]

#print(FILES)

#source the R files
file.sources = list.files(path=paste0(RSRC), pattern="*.R$")
sapply(paste0(RSRC, "/", file.sources),source)

gdist.inputs <- readRDS(paste0(WD, "/", "gdist-inputs.rds"))
GA.inputs <- readRDS(paste0(WD, "/", "GA-inputs.rds"))

#gather SS results
ss_full_results <- SS_optim_gather(gdist.inputs = gdist.inputs,
                                   GA.inputs = GA.inputs,
                                   nlm = FALSE,
                                   dist_mod = TRUE,
                                   null_mod = TRUE,
                                   max.combination=12,
                                   results_ss=WD)
saveRDS(ss_full_results, file="single_surface.RESULTS_ALL.rds")
