
# Optimizing Resistance Surface using ResistanceGA
# Need packages: ResistanceGA , rgdal, otuSummary
# Setting up the workspace
suppressWarnings(require(raster))
suppressWarnings(require(ResistanceGA))
#setwd("C:/Users/zbind/Desktop/ResistanceGA_Examples/WTD")
setwd(%BASE)
all.comb <- "C:/Users/zdzbinde/Desktop/ResistanceGA_Examples/WTD/all.comb/"
#CS.program <- paste('"C:/Program Files/Circuitscape/cs_run.exe"')


# loading the data
files <- list.files(path=".", pattern="*.proj.tif$", full.names=TRUE, recursive=FALSE)
surfaces <- raster::stack(paste0("C:/Users/zdzbinde/Desktop/ResistanceGA_Examples/WTD", files))
samples <- read.table(file = 'wtd_mz_list.tsv', sep = '\t', header = TRUE)
sample.locales <- SpatialPoints(samples[,c(14,15)])

gen.dist <- read.table(file = 'wtd_mz_prov_dist.tsv', sep = '\t', header = TRUE)
gen.dist <- otuSummary::matrixConvert(gen.dist, colname = c("ind1", "ind2", "dist"))
gen.dist <- gen.dist[,3]


# Running the Resistance Surface Optimization
# Format inputs
GA.inputs <- GA.prep(ASCII.dir = surfaces,
                     Results.dir = "all_comb",
                     method = 'LL',
                     max.cat = 500,
                     max.cont = 500,
                     seed = NULL,
                     parallel = TRUE,
                     quiet = FALSE,
                     # this value should be justified
                     k.value = 3)

gdist.inputs <- gdist.prep(length(sample.locales),
                           samples = sample.locales,
                           response = gen.dist,
                           method = 'commuteDistance')


# Runs single surface opitimization for each surface independently followed by optimization and
# model selection of combination surfaces
comb.analysis <- all_comb(gdist.inputs,
                          GA.inputs,
                          all.comb,
                          max.combination = 4,
                          iters = 1000,
                          replicate = 1,
                          sample.prop = 0.75,
                          nlm = FALSE,
                          dist_mod = TRUE,
                          null_mod = TRUE)
