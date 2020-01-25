
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
WD<-"~/ResistanceGA_HPC_Pipeline/work/all_comb"
RSRC<-args[2] #path to R scripts to source
SS <- args[3]
GD <- args[4]
GA <- args[5]
iters<-args[6]

#source the R files
file.sources = list.files(path=paste0(RSRC), pattern="*.R$")
sapply(paste0(RSRC, "/", file.sources),source)

gdist.inputs <- readRDS(paste0(WD, "/", GD))
GA.inputs <- readRDS(paste0(WD, "/", GA))

#load all ms results
ms.results <- list.files(full.names = TRUE, path=WD, 
                       pattern = "*.MS.rds$") %>%
                   map_dfr(readRDS) %>% 
                   bind_rows()

ss.results <- readRDW(paste0(WD, "/", SS))

ms_results_gather(gdist.inputs,
                          GA.inputs,
                          ms.results=ms.results,
                          ss/results=ss.results)
