library(tidyverse)

#RGA modified functions

SS_optim_worker <- function(CS.inputs = NULL,
                            gdist.inputs = NULL,
                            GA.inputs,
                            nlm = FALSE,
                            dist_mod = TRUE,
                            null_mod = TRUE) {
  if (!is.null(GA.inputs$scale)) {
    stop(
      "This function should NOT be used if you intend to apply kernel smoothing to your resistance surfaces"
    )
  }
  results_ss <- "."
  Plots.dir<-GA.inputs$Plots.dir

  if(!file.exists(results_ss)){
    dir.create(results_ss, recursive=T)
  }
  
  #RESULTS.cat <- list() # List to store categorical results within
  #RESULTS.cont <- list() # List to store continuous results within
  cnt1 <- 0
  cnt2 <- 0
  k.value <- GA.inputs$k.value
  #MLPE.list <- list()
  #cd.list <- list()
  #k.list <- list()
  
  # Optimize each surface in turn
  for (i in 1:GA.inputs$n.layers) {
    #if i not in list of job indices, skip
      r <- GA.inputs$Resistance.stack[[i]]
      names(r) <- GA.inputs$layer.names[i]
      print(paste0("Starting SS optimization on layer ", GA.inputs$layer.names[i]))
      
      # CIRCUITSCAPE ------------------------------------------------------------
      
      
      # * Categorical -----------------------------------------------------------
      
      # Processing of categorical surfaces
      if (!is.null(CS.inputs)) {
        if (GA.inputs$parallel != FALSE) {
          warning(
            "\n CIRCUITSCAPE cannot be optimized in parallel. \n Ignoring parallel arguement. \n If you want to optimize in parallel, use least cost paths and gdistance.",
            immediate. = TRUE
          )
        }
        if (GA.inputs$surface.type[i] == 'cat') {
          cnt1 <- cnt1 + 1
          names(r) <- GA.inputs$layer.names[i]
          
          single.GA <- ga(
            type = "real-valued", fitness = Resistance.Opt_single, Resistance = r,
            population = GA.inputs$population, selection = GA.inputs$selection,  pcrossover = GA.inputs$pcrossover,
            pmutation = GA.inputs$pmutation, crossover = GA.inputs$crossover, Min.Max = GA.inputs$Min.Max, GA.inputs = GA.inputs,
            CS.inputs = CS.inputs, lower = GA.inputs$min.list[[i]], upper = GA.inputs$max.list[[i]],
            popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]), maxiter = GA.inputs$maxiter,
            run = GA.inputs$run, keepBest = GA.inputs$keepBest, elitism = GA.inputs$percent.elite,
            mutation = GA.inputs$mutation, seed = GA.inputs$seed, iter = i, quiet = GA.inputs$quiet )
          
          if(dim(single.GA@solution)[1] > 1) {
            single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
          }
          
          single.GA@solution <-
            single.GA@solution / min(single.GA@solution)
          df <- data.frame(id = unique(r), t(single.GA@solution))
          r <- subs(r, df)
          names(r) <- GA.inputs$layer.names[i]
          
          Run_CS(CS.inputs, GA.inputs, r, EXPORT.dir = results_ss)
          
          Diagnostic.Plots(
            resistance.mat = paste0( results_ss, GA.inputs$layer.names[i], "_resistances.out"),
            genetic.dist = CS.inputs$response, plot.dir = GA.inputs$Plots.dir, type = "categorical", ID = CS.inputs$ID, ZZ = CS.inputs$ZZ)
          
          fit.stats <-
            r.squaredGLMM(
              MLPE.lmm(
                resistance = paste0(
                  results_ss,
                  GA.inputs$layer.names[i],
                  "_resistances.out"
                ),
                pairwise.genetic = CS.inputs$response,
                REML = F,
                ID = CS.inputs$ID,
                ZZ = CS.inputs$ZZ
              )
            )
          
          aic <-
            AIC(
              MLPE.lmm(
                resistance = paste0(
                  results_ss,
                  GA.inputs$layer.names[i],
                  "_resistances.out"
                ),
                pairwise.genetic = CS.inputs$response,
                REML = F,
                ID = CS.inputs$ID,
                ZZ = CS.inputs$ZZ
              )
            )
          
          LL <-
            logLik(
              MLPE.lmm(
                resistance = paste0(
                  results_ss,
                  GA.inputs$layer.names[i],
                  "_resistances.out"
                ),
                pairwise.genetic = CS.inputs$response,
                REML = F,
                ID = CS.inputs$ID,
                ZZ = CS.inputs$ZZ
              )
            )
          
          if (k.value == 1) {
            k <- 2
          } else if (k.value == 2) {
            k <- GA.inputs$parm.type$n.parm[i] + 1
          } else if (k.value == 3) {
            k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
          } else {
            k <- length(GA.inputs$layer.names[i]) + 1
          }
          
          n <- CS.inputs$n.Pops
          AICc <-
            (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
          # AICc <- (aic)+(((2*k)*(k+1))/(CS.inputs$n.Pops-k-1))
          
          RS <- data.frame(
            GA.inputs$layer.names[i],
            single.GA@fitnessValue,
            k,
            aic,
            AICc,
            fit.stats[[1]],
            fit.stats[[2]],
            LL[[1]],
            single.GA@solution
          )
          
          k <- GA.inputs$parm.type$n.parm[i]
          
          Features <- matrix()
          for (z in 1:(k)) {
            feature <- paste0("Feature", z)
            Features[z] <- feature
          }
          
          colnames(RS) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              "k",
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL",
              Features
            )
          
          #RESULTS.cat[[cnt1]] <- RS
          saveRDS(RS, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_RESULTS_CAT.SS.rds"), compress=TRUE)
          
          MLPE <- MLPE.lmm(
            resistance = paste0(
              results_ss,
              GA.inputs$layer.names[i],
              "_resistances.out"
            ),
            pairwise.genetic = CS.inputs$response,
            REML = TRUE,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
          saveRDS(MLPE, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_MLPE.SS.rds"), compress=TRUE)
          
          CD <- (read.table(paste0(
            results_ss,
            GA.inputs$layer.names[i],
            "_resistances.out"))[-1, -1])
          saveRDS(CD, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_CD.SS.rds"), compress=TRUE)
          
          #names(MLPE.list)[i] <- GA.inputs$layer.names[i]
          #names(cd.list)[i] <- GA.inputs$layer.names[i]
          
          # * Continuous -----------------------------------------------------------
        } else {
          # Processing of continuous surfaces
          cnt2 <- cnt2 + 1
          r <- SCALE(r, 0, 10)
          names(r) <- GA.inputs$layer.names[i]
          
          single.GA <- ga(
            type = "real-valued",
            fitness = Resistance.Opt_single,
            Resistance = r,
            population = GA.inputs$population,
            selection = GA.inputs$selection,
            pcrossover = GA.inputs$pcrossover,
            pmutation = GA.inputs$pmutation,
            crossover = GA.inputs$crossover,
            Min.Max = GA.inputs$Min.Max,
            GA.inputs = GA.inputs,
            CS.inputs = CS.inputs,
            lower = GA.inputs$min.list[[i]],
            upper = GA.inputs$max.list[[i]],
            popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]),
            maxiter = GA.inputs$maxiter,
            run = GA.inputs$run,
            keepBest = GA.inputs$keepBest,
            elitism = GA.inputs$percent.elite,
            mutation = GA.inputs$mutation,
            seed = GA.inputs$seed,
            iter = i,
            quiet = GA.inputs$quiet
          )
          
          if(dim(single.GA@solution)[1] > 1) {
            single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
          }
          
          # Using GA results, optimize with nlm
          start.vals <- single.GA@solution[-1]
          
          if (nlm == TRUE) {
            # Informed start values; these are the optimized values from the single parameter optimization
            EQ <- get.EQ(single.GA@solution[1])
            Optim.nlm <-
              nlm(
                Resistance.Optimization_cont.nlm,
                log(start.vals),
                Resistance = r,
                equation = single.GA@solution[1],
                get.best = FALSE,
                CS.inputs = CS.inputs,
                Min.Max = 'min',
                write.dir = results_ss
              )
            
            OPTIM <-
              Resistance.Optimization_cont.nlm(
                PARM = (Optim.nlm$estimate),
                Resistance = r,
                equation = single.GA@solution[1],
                get.best = TRUE,
                CS.inputs = CS.inputs,
                Min.Max = 'min',
                write.dir = results_ss
              )
            
            Diagnostic.Plots(
              resistance.mat = paste0(
                results_ss,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              genetic.dist = CS.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
            
            Plot.trans(
              PARM = exp(Optim.nlm$estimate),
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir,
              Name = GA.inputs$layer.names[i]
            )
            
            RS <-
              data.frame(GA.inputs$layer.names[i],
                         Optim.nlm$minimum,
                         EQ,
                         Cont.Param(exp(Optim.nlm$estimate)))
            colnames(RS) <-
              c("Surface",
                "obj.func",
                "R2m",
                "R2c",
                "Equation",
                "shape",
                "max")
            #RESULTS.cont[[cnt2]] <- RS
            saveRDS(RS, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_RESULTS_CONT.SS.rds"), compress=TRUE)
            
            MLPE <- MLPE.lmm(
              resistance = paste0(
                results_ss,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              pairwise.genetic = CS.inputs$response,
              REML = TRUE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
            saveRDS(MLPE, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_MLPE.SS.rds"), compress=TRUE)
            
            CD <- (read.table(paste0(
              results_ss,
              GA.inputs$layer.names[i],
              "_resistances.out"))[-1, -1])
            saveRDS(CD, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_CD.SS.rds"), compress=TRUE)
            
            #names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            #names(cd.list)[i] <- GA.inputs$layer.names[i]
            
          } else {
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r.tran <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r.tran) <- GA.inputs$layer.names[i]
            
            Run_CS(CS.inputs,
                   GA.inputs,
                   r.tran,
                   EXPORT.dir = results_ss)
            
            Diagnostic.Plots(
              resistance.mat = paste0(
                results_ss,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              genetic.dist = CS.inputs$response,
              plot.dir = Plots.dir,
              type = "continuous",
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = Plots.dir
            )
            
            fit.stats <- r.squaredGLMM(
              MLPE.lmm(
                REML = F,
                resistance = paste0(
                  results_ss,
                  GA.inputs$layer.names[i],
                  "_resistances.out"
                ),
                pairwise.genetic = CS.inputs$response,
                ID = CS.inputs$ID,
                ZZ = CS.inputs$ZZ
              )
            )
            aic <- AIC(
              MLPE.lmm(
                REML = F,
                resistance = paste0(
                  results_ss,
                  GA.inputs$layer.names[i],
                  "_resistances.out"
                ),
                pairwise.genetic = CS.inputs$response,
                ID = CS.inputs$ID,
                ZZ = CS.inputs$ZZ
              )
            )
            
            LL <-
              logLik(
                MLPE.lmm(
                  resistance = paste0(
                    results_ss,
                    GA.inputs$layer.names[i],
                    "_resistances.out"
                  ),
                  pairwise.genetic = CS.inputs$response,
                  REML = F,
                  ID = CS.inputs$ID,
                  ZZ = CS.inputs$ZZ
                )
              )
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            #k.list[[i]] <- k
            #names(k.list)[i] <- GA.inputs$layer.names[i]
            saveRDS(k, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_K.SS.rds"), compress=TRUE)
            
            n <- CS.inputs$n.Pops
            AICc <-
              (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
            # AICc <- (aic)+(((2*k)*(k+1))/(CS.inputs$n.Pops-k-1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            #RESULTS.cont[[cnt2]] <- RS
            saveRDS(RS, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_RESULTS_CONT.SS.rds"), compress=TRUE)
            
            MLPE <- MLPE.lmm(
              resistance = paste0(
                results_ss,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              pairwise.genetic = CS.inputs$response,
              REML = TRUE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
            saveRDS(MLPE, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_MLPE.SS.rds"), compress=TRUE)
            
            CD <- (read.table(paste0(
              results_ss,
              GA.inputs$layer.names[i],
              "_resistances.out"))[-1, -1])
            saveRDS(CD, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_CD.SS.rds"), compress=TRUE)
            
            #names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            #names(cd.list)[i] <- GA.inputs$layer.names[i]
            
          }
        } # Close if-else
      }
      
      # Optimize with gdistance -------------------------------------------------
      
      if (!is.null(gdist.inputs)) {
        # * Categorical -----------------------------------------------------------
        
        if (GA.inputs$surface.type[i] == 'cat') {
          cnt1 <- cnt1 + 1
          names(r) <- GA.inputs$layer.names[i]
          
          single.GA <- ga(
            type = "real-valued",
            fitness = Resistance.Opt_single,
            Resistance = r,
            population = GA.inputs$population,
            selection = GA.inputs$selection,
            pcrossover = GA.inputs$pcrossover,
            pmutation = GA.inputs$pmutation,
            crossover = GA.inputs$crossover,
            Min.Max = GA.inputs$Min.Max,
            GA.inputs = GA.inputs,
            gdist.inputs = gdist.inputs,
            lower = GA.inputs$min.list[[i]],
            upper = GA.inputs$max.list[[i]],
            parallel = GA.inputs$parallel,
            popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]),
            maxiter = GA.inputs$maxiter,
            run = GA.inputs$run,
            keepBest = GA.inputs$keepBest,
            elitism = GA.inputs$percent.elite,
            mutation = GA.inputs$mutation,
            seed = GA.inputs$seed,
            iter = i,
            quiet = GA.inputs$quiet
          )
          
          if(dim(single.GA@solution)[1] > 1) {
            single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
          }
          
          single.GA@solution <-
            single.GA@solution / min(single.GA@solution)
          df <- data.frame(id = unique(r), t(single.GA@solution))
          r <- subs(r, df)
          NAME <- GA.inputs$layer.names[i]
          names(r) <- NAME
          
          cd <- Run_gdistance(gdist.inputs, r)
          # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
          write.table(
            as.matrix(cd),
            file = paste0(results_ss, "/", NAME, "_", gdist.inputs$method,  "_distMat.csv"),
            
            sep = ",",
            row.names = F,
            col.names = F
          )
          writeRaster(r,
                      paste0(results_ss, "/", NAME, ".asc"),
                      overwrite = TRUE)
          
          # save(single.GA,
          #      file = paste0(results_ss, NAME, ".rda"))
          
          saveRDS(single.GA,
                  file = paste0(results_ss, "/", NAME, ".rds"))
          
          Diagnostic.Plots(
            resistance.mat = cd,
            genetic.dist = gdist.inputs$response,
            plot.dir = Plots.dir,
            type = "categorical",
            name = NAME,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ
          )
          
          fit.stats <- r.squaredGLMM(
            MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = F,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
          )
          
          aic <- AIC(
            MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = F,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
          )
          
          LL <- logLik(
            MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = F,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
          )
          
          if (k.value == 1) {
            k <- 2
          } else if (k.value == 2) {
            k <- GA.inputs$parm.type$n.parm[i] + 1
          } else if (k.value == 3) {
            k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
          } else {
            k <- length(GA.inputs$layer.names[i]) + 1
          }
          
          #k.list[[i]] <- k
          #names(k.list)[i] <- GA.inputs$layer.names[i]
          saveRDS(k, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_K.SS.rds"), compress=TRUE)
          
          n <- gdist.inputs$n.Pops
          AICc <-
            (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
          # AICc <- (aic)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
          
          
          RS <- data.frame(
            GA.inputs$layer.names[i],
            single.GA@fitnessValue,
            k,
            aic,
            AICc,
            fit.stats[[1]],
            fit.stats[[2]],
            LL[[1]],
            single.GA@solution
          )
          
          k <- GA.inputs$parm.type$n.parm[i]
          
          Features <- matrix()
          for (z in 1:(k)) {
            feature <- paste0("Feature", z)
            Features[z] <- feature
          }
          
          colnames(RS) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              "k",
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL",
              Features
            )
          
          #RESULTS.cat[[cnt1]] <- RS
          saveRDS(RS, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_RESULTS_CAT.SS.rds"), compress=TRUE)
          
          MLP <-  MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            REML = TRUE,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ
          )
          saveRDS(MLPE, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_MLPE.SS.rds"), compress=TRUE)
          
          CD <- as.matrix(cd)
          saveRDS(CD, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_CD.SS.rds"), compress=TRUE)
          #names(cd.list)[i] <- GA.inputs$layer.names[i]
          
          #names(MLPE.list)[i] <- GA.inputs$layer.names[i]
          
        } else {
          # Processing of continuous surfaces
          cnt2 <- cnt2 + 1
          r <- SCALE(r, 0, 10)
          names(r) <- GA.inputs$layer.names[i]
          
          single.GA <- ga(
            type = "real-valued",
            fitness = Resistance.Opt_single,
            Resistance = r,
            population = GA.inputs$population,
            selection = GA.inputs$selection,
            pcrossover = GA.inputs$pcrossover,
            pmutation = GA.inputs$pmutation,
            crossover = GA.inputs$crossover,
            Min.Max = GA.inputs$Min.Max,
            GA.inputs = GA.inputs,
            gdist.inputs = gdist.inputs,
            lower = GA.inputs$min.list[[i]],
            upper = GA.inputs$max.list[[i]],
            parallel = GA.inputs$parallel,
            popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]),
            maxiter = GA.inputs$maxiter,
            run = GA.inputs$run,
            keepBest = GA.inputs$keepBest,
            elitism = GA.inputs$percent.elite,
            mutation = GA.inputs$mutation,
            seed = GA.inputs$seed,
            iter = i,
            quiet = GA.inputs$quiet
          )
          
          # ** Second optim ---------------------------------------------------------
          
          
          # Using GA results, optimize with nlm
          start.vals <- single.GA@solution[-1]
          
          if (nlm == TRUE) {
            # Informed start values; these are the optimized values from the single parameter optimization
            EQ <- get.EQ(single.GA@solution[1])
            Optim.nlm <-
              nlm(
                Resistance.Optimization_cont.nlm,
                log(start.vals),
                Resistance = r,
                equation = single.GA@solution[1],
                get.best = FALSE,
                gdist.inputs = gdist.inputs,
                Min.Max = 'min',
                write.dir = results_ss
              )
            
            r <-
              Resistance.tran(
                transformation = EQ,
                shape = Optim.nlm$estimate[1],
                max = Optim.nlm$estimate[2],
                r = r,
                out = results_ss
              )
            
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(results_ss, "/", NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(results_ss, "/", NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA,
            #      file = paste0(results_ss, NAME, ".rda"))
            
            saveRDS(single.GA,
                    file = paste0(results_ss, "/", NAME, ".rds"))
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = exp(Optim.nlm$estimate),
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = Plots.dir,
              Name = GA.inputs$layer.names[i]
            )
            
            fit.stats <-
              r.squaredGLMM(
                MLPE.lmm(
                  resistance = cd,
                  pairwise.genetic = gdist.inputs$response,
                  REML = F,
                  ID = gdist.inputs$ID,
                  ZZ = gdist.inputs$ZZ
                )
              )
            
            
            aic <-
              AIC(
                MLPE.lmm(
                  resistance = cd,
                  pairwise.genetic = gdist.inputs$response,
                  REML = F,
                  ID = gdist.inputs$ID,
                  ZZ = gdist.inputs$ZZ
                )
              )
            
            RS <- data.frame(
              GA.inputs$layer.names[i],-single.GA@fitnessValue,
              fit.stats[[1]],
              fit.stats[[2]],
              single.GA@solution
            )
            
            # RS<-data.frame(GA.inputs$layer.names[i],Optim.nlm$minimum,EQ,Cont.Param(exp(Optim.nlm$estimate)))
            colnames(RS) <-
              c("Surface",
                "AICc",
                "R2m",
                "R2c",
                "Equation",
                "shape",
                "max")
            #RESULTS.cont[[cnt2]] <- RS
            saveRDS(RS, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_RESULTS_CONT.SS.rds"), compress=TRUE)
            
            # * Continuous -----------------------------------------------------------
            
          } else {
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(results_ss, "/", NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(results_ss, "/", NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA,
            #      file = paste0(results_ss, NAME, ".rda"))
            
            saveRDS(single.GA,
                    file = paste0(results_ss, "/", NAME, ".rds"))
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = Plots.dir
            )
            
            fit.stats <-
              r.squaredGLMM(
                MLPE.lmm(
                  resistance = cd,
                  pairwise.genetic = gdist.inputs$response,
                  REML = F,
                  ID = gdist.inputs$ID,
                  ZZ = gdist.inputs$ZZ
                )
              )
            
            
            aic <-
              AIC(
                MLPE.lmm(
                  resistance = cd,
                  pairwise.genetic = gdist.inputs$response,
                  REML = F,
                  ID = gdist.inputs$ID,
                  ZZ = gdist.inputs$ZZ
                )
              )
            
            LL <-
              logLik(
                MLPE.lmm(
                  resistance = cd,
                  pairwise.genetic = gdist.inputs$response,
                  REML = F,
                  ID = gdist.inputs$ID,
                  ZZ = gdist.inputs$ZZ
                )
              )
            
            MLPE <-  MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = TRUE,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            saveRDS(MLPE, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_MLPE.SS.rds"), compress=TRUE)
            
            CD <- as.matrix(cd)
            saveRDS(CD, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_CD.SS.rds"), compress=TRUE)
            #names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            #names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            #k.list[[i]] <- k
            #names(k.list)[i] <- GA.inputs$layer.names[i]
            saveRDS(k, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_K.SS.rds"), compress=TRUE)
            
            n <- gdist.inputs$n.Pops
            AICc <-
              (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
            # AICc <- (aic)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            #RESULTS.cont[[cnt2]] <- RS
            saveRDS(RS, file=paste0(results_ss, "/", GA.inputs$layer.names[i], "_RESULTS_CONT.SS.rds"), compress=TRUE)
          }
        } # Close if-else
    }
  } # Close ascii loop
}





####################################################
#
# Above is the "worker" code... Needs to write some Rdata
#
# for the below "gather" code to then load.
#
#
#
####################################################



SS_optim_gather <- function(CS.inputs = NULL,
                            gdist.inputs = NULL,
                            GA.inputs,
                            nlm = FALSE,
                            dist_mod = TRUE,
                            null_mod = TRUE,
                            max.combination,
                            results_ss = NULL) {
  
  
  if(!exists('gdist.inputs'))
    return(cat("ERROR: Please specify gdist.inputs"))
  
  if(!exists('GA.inputs'))
    return(cat("ERROR: Please specify GA.inputs"))
  
  GA.input_orig <- GA.inputs
  
  Results <- vector(mode = 'list', length = 1)
  
  # Make results data frame
  Results.cat <- data.frame()
  Results.cont <- data.frame()
  # cnt1<-0
  # cnt2<-0

  Results.cat<- list.files(full.names = TRUE, path=results_ss, 
                           pattern = "*_RESULTS_CAT.SS.rds$") %>%
                map_dfr(readRDS) %>% 
                bind_rows()
  
  Results.cont<- list.files(full.names = TRUE, path=results_ss, 
                           pattern = "*_RESULTS_CONT.SS.rds$") %>%
                map_dfr(readRDS) %>% 
                bind_rows()
  

  ##################################
  # Compile results into tables
  cat("\n")
  cat("\n")
  if (nrow(Results.cat) > 0) {
    Features <- array()
    for (i in 1:ncol(Results.cat) - 8) {
      feature <- paste0("Feature", i)
      Features[i] <- feature
    }
    colnames(Results.cat) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        'k',
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL",
        Features
      )
    Results.cat <-  Results.cat[order(Results.cat$AICc), ]
    write.table(
      Results.cat,
      paste0(results_ss, "/", "CategoricalResults.csv"),
      sep = ",",
      col.names = T,
      row.names = F
    )
  }
  
  if (ncol(Results.cont) > 0) {
    colnames(Results.cont) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        'k',
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL",
        "Equation",
        "shape",
        "max"
      )
    Results.cont <- Results.cont[order(Results.cont$AICc), ]
    write.table(
      Results.cont,
      paste0(results_ss, "/", "ContinuousResults.csv"),
      sep = ",",
      col.names = T,
      row.names = F
    )
  }
  
  # Full Results
  if (nrow(Results.cat) > 0 & nrow(Results.cont) > 0) {
    Results.All <- rbind(Results.cat[, c(1:8)], Results.cont[, c(1:8)])
  } else if (nrow(Results.cat) < 1 & nrow(Results.cont) > 0) {
    Results.All <- (Results.cont[, c(1:8)])
  } else {
    Results.All <- (Results.cat[, c(1:8)])
  }

  if (dist_mod == TRUE) {
    r <- GA.inputs$Resistance.stack[[1]]
    r <- reclassify(r, c(-Inf, Inf, 1))
    names(r) <- "dist"
    cd <- Run_gdistance(gdist.inputs, r)
    write.table( as.matrix(cd), file = paste0(results_ss, "/", 'Distance', "_", gdist.inputs$method, "_distMat.csv"),  sep = ",", row.names = F, col.names = F )
    Dist.AIC <- suppressWarnings(AIC(MLPE.lmm2( resistance = cd, response = gdist.inputs$response, ID = gdist.inputs$ID, ZZ = gdist.inputs$ZZ, REML = FALSE )))
    fit.stats <- r.squaredGLMM(MLPE.lmm2(resistance = cd, response = gdist.inputs$response, ID = gdist.inputs$ID, ZZ = gdist.inputs$ZZ, REML = FALSE))
    LL <- logLik(MLPE.lmm2(resistance = cd, response = gdist.inputs$response, ID = gdist.inputs$ID, ZZ = gdist.inputs$ZZ, REML = FALSE))
    dist.MLPE <-  MLPE.lmm2(resistance = cd, response = gdist.inputs$response, REML = TRUE, ID = gdist.inputs$ID, ZZ = gdist.inputs$ZZ)
    dist.CD<- as.matrix(cd)
    ROW <- nrow(gdist.inputs$ID)
    dist.k <- 2
    if (GA.inputs$method == "AIC") {
      dist.obj <- -Dist.AIC
    } else if (GA.inputs$method == "R2") {
      dist.obj <- fit.stats[[1]]
    } else {
      dist.obj <- LL[[1]]
    }
    k<-dist.k
    n <- gdist.inputs$n.Pops
    AICc <-
      (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
    # AICc <- (Dist.AIC)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
    Dist.AICc <- data.frame("Distance", dist.obj, dist.k, Dist.AIC, AICc, fit.stats[[1]], fit.stats[[2]], LL[[1]])
    colnames(Dist.AICc) <- c("Surface", paste0("obj.func_", GA.inputs$method), 'k', "AIC", "AICc", "R2m", "R2c", "LL")
    Results.All <- rbind(Results.All, Dist.AICc)
  }
  if (null_mod == TRUE){
    dat <- data.frame(gdist.inputs$ID, response = gdist.inputs$response)
    colnames(dat) <- c("pop1", "pop2", "response")
    
    # Fit model
    mod <-lFormula(response ~ 1 + (1 | pop1), data = dat, REML = FALSE)
    mod$reTrms$Zt <- gdist.inputs$ZZ
    dfun <- do.call(mkLmerDevfun, mod)
    opt <- optimizeLmer(dfun)
    Null.AIC <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    fit.stats <- r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    LL <- logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    ROW <- nrow(gdist.inputs$ID)
    k <- 1
    if (GA.inputs$method == "AIC") {
      null.obj <- -Null.AIC
    } else if (GA.inputs$method == "R2") {
      null.obj <- fit.stats[[1]]
    } else {
      null.obj <- LL[[1]]
    }
    n <- gdist.inputs$n.Pops
    AICc <- (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
    # AICc <- (Null.AIC)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
    Null.AICc <- data.frame("Null", null.obj, k, Null.AIC, AICc, fit.stats[[1]], fit.stats[[2]], LL[[1]])
    colnames(Null.AICc) <- c( "Surface", paste0("obj.func_", GA.inputs$method), 'k', "AIC", "AICc", "R2m", "R2c", "LL")
    Results.All <- rbind(Results.All, Null.AICc)
  }
  
  Results.All <- Results.All[order(Results.All$AICc), ]
  cat("\n")
  cat("\n")
  write.table( Results.All, paste0(results_ss, "/", "All_Results_Table_", gdist.inputs$method,".csv"), sep = ",", col.names = T, row.names = F)
  
  # Get parameter estimates
  if (!is.null(CS.inputs)) {
    MLPE.results <- MLPE.lmm_coef(
      resistance = results_ss,
      genetic.dist = CS.inputs$response,
      out.dir = results_ss,
      method = "cs",
      ID = CS.inputs$ID,
      ZZ = CS.inputs$ZZ
    )
    
  } else {
    MLPE.results <- MLPE.lmm_coef(
      resistance = results_ss,
      genetic.dist = gdist.inputs$response,
      out.dir = results_ss,
      method = "gd",
      ID = gdist.inputs$ID,
      ZZ = gdist.inputs$ZZ
    )
  }
  
  read_stuff <- function(f, p){
    dat<- readRDS(f)
    name<-str_replace(basename(f), p, "")
    return(data.frame(name, dat))
  }
  read_stuff_name <- function(f, p){
    name<-str_replace(basename(f), p, "")
    return(name)
  }
  read_stuff_dat <- function(f, p){
    dat<- readRDS(f)
    name<-str_replace(basename(f), p, "")
    return(dat)
  }
  
  k.list <- list.files(full.names = TRUE, path=results_ss, 
                            pattern = "*_K.SS.rds$")
  k.list <- lapply(k.list, read_stuff, "_K.SS.rds")
  k.list <- plyr::ldply(k.list)
  colnames(k.list) <- c("surface", "k")
  
  cd.list <- list.files(full.names = TRUE, path=results_ss, 
                       pattern = "*_CD.SS.rds$")
  vals <- lapply(cd.list, read_stuff_dat, "_CD.SS.rds")
  names <- lapply(cd.list, read_stuff_name, "_CD.SS.rds")
  cd.list<-vals
  names(cd.list)<-names
  
  MLPE.list <- list.files(full.names = TRUE, path=results_ss, 
                        pattern = "*_MLPE.SS.rds$")
  vals <- lapply(MLPE.list, read_stuff_dat, "_MLPE.SS.rds")
  names <- lapply(MLPE.list, read_stuff_name, "_MLPE.SS.rds")
  MLPE.list<-vals
  names(MLPE.list)<-names
  
  
  if (dist_mod == TRUE) {
    newdf <- data.frame("Distance", dist.k)
    names(newdf)<-c("surface", "k")
    k.list<-rbind(k.list, newdf)
    
    cd.list <- c(cd.list, dist.CD) 
    names(cd.list) <- c(names, "Distance")
    
    MLPE.list <- c(MLPE.list, dist.MLPE) 
    names(MLPE.list) <- c(names, "Distance")
  }
  
  saveRDS(k.list, file=paste0(results_ss, "/single_surface_K.ALL.RDS"))
  saveRDS(MLPE.list, file=paste0(results_ss, "/single_surface_MLPE.ALL.RDS"))
  saveRDS(cd.list, file=paste0(results_ss, "/single_surface_CD.ALL.RDS"))
  
  # Full Results
  if (nrow(Results.cat) > 0 & nrow(Results.cont) > 0) {
    RESULTS <-
      list(
        ContinuousResults = Results.cont,
        CategoricalResults = Results.cat,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = NULL,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list
      )
    
  } else if (nrow(Results.cat) < 1 & nrow(Results.cont) > 0) {
    RESULTS <-
      list(
        ContinuousResults = Results.cont,
        CategoricalResults = NULL,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list
      )
    
  } else if (nrow(Results.cat) > 0 & nrow(Results.cont) < 1) {
    RESULTS <-
      list(
        ContinuousResults = NULL,
        CategoricalResults = Results.cat,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list
      )
  } else {
    RESULTS <-
      list(
        ContinuousResults = NULL,
        CategoricalResults = NULL,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list
      )
  }
  
  unlink(GA.inputs$Write.dir, recursive = T, force = T)
  return(RESULTS)
  ###############################################################################################################
}











#first step, parallel single-surface optimizations
all_comb_pipeline_SS <- function(gdist.inputs,
                                 GA.inputs,
                                 results.dir,
                                 procnum,
                                 procs,
                                 max.combination = 4,
                                 iters = 1000,
                                 replicate = 1,
                                 sample.prop = 0.75,
                                 nlm = FALSE,
                                 dist_mod = TRUE,
                                 null_mod = TRUE,
                                 Plots.dir) {
  
  
  if(!exists('gdist.inputs')) 
    return(cat("ERROR: Please specify gdist.inputs"))
  
  if(!exists('GA.inputs')) 
    return(cat("ERROR: Please specify GA.inputs"))
  
  
  if(length(max.combination) > 2) {
    return(cat("ERROR: Please specify either a single value or a vector with the minimum and maximum value"))
  }
  
  # if(length(max.combination) == 2) {
  #   if(max.combination[2] > GA.inputs$n.layers) {
  #     return(cat("ERROR: Please specify a maximum combination that is less than or equal to the number of raster layers in the analysis"))
  #   }
  # }
  
  
  # Create combination list -------------------------------------------------
  mc <- max.combination
  
  if(length(max.combination) == 2) {
    if(mc[1] == 1) {
      min.combination <- 2
      max.combination <- mc[2]
      ss <- TRUE
    } else {
      min.combination <- mc[1]
      max.combination <- mc[2]
      ss <- FALSE
    } 
  } else {
    min.combination <- 2
    ss <- TRUE
  }
  
  if(max.combination > GA.inputs$n.layers) {
    max.combination <- GA.inputs$n.layers
  }
  comb.list <- vector(mode = "list", length = (max.combination - 1))
  
  
  GA.input_orig <- GA.inputs
  
  Results <- vector(mode = 'list', length = replicate)
  
  # Begin Replicate Loop --------------------------------------------------

    # Skip if min combination > 1
    if(ss == FALSE) {
      ss.results <- NULL
      AICc.tab <- NULL
      #dir.create(paste0(results.dir,'rep_',i))
    } else {  # Do single surface optimization
      
      # Single Surface optimization -----------------------------------------------------
      if(!is.null(GA.inputs$scale)) {
        stop(
          "all_comb_pipeline_SS  NOT WORKING YET FOR SCALED SURFACE OPTIMIZATION. EMAIL tkchafin@uark.edu FOR MORE INFORMATION"
        )
        # * Single Surface: scaled --------------------------------------------------------
        
        ss.results <- SS_optim.scale(gdist.inputs = gdist.inputs,
                                     GA.inputs = GA.inputs,
                                     nlm = nlm,
                                     dist_mod = dist_mod,
                                     null_mod = null_mod)
        
        AICc.tab <- ss.results$AICc
      } else {
        # * Single Surface --------------------------------------------------------
        
        ss.results <- SS_optim_worker(gdist.inputs = gdist.inputs,
                                      GA.inputs = GA.inputs,
                                      nlm = nlm,
                                      dist_mod = dist_mod,
                                      null_mod = null_mod)
        #NEXT: Call SS_optim_gather
        #THEN: Continue pipeline for MS optimizations
        AICc.tab <- ss.results$AICc
      }
    }
}



#########################################
#
# all_comb_pipeline_MS
#
#
#########################################
all_comb_pipeline_make_batches <- function(gdist.inputs,
                                 GA.inputs,
                                 max.combination) {
  
  if(length(max.combination) > 2) {
    return(cat("ERROR: Please specify either a single value or a vector with the minimum and maximum value"))
  }
  
  # Create combination list -------------------------------------------------
  max.combination<-as.numeric(max.combination)
  mc <- max.combination
  
  if(length(max.combination) == 2) {
    if(mc[1] == 1) {
      min.combination <- 2
      max.combination <- mc[2]
      ss <- TRUE
    } else {
      min.combination <- mc[1]
      max.combination <- mc[2]
      ss <- FALSE
    } 
  } else {
    min.combination <- 2
    ss <- TRUE
  }
  
  if(max.combination > GA.inputs$n.layers) {
    max.combination <- GA.inputs$n.layers
  }
  comb.list <- vector(mode = "list", length = (max.combination - 1))
  
  list.count <- 0
  surface.count <- 0
  for(i in min.combination:max.combination) {
    list.count <- list.count + 1
    print(max.combination)
    comb.list[[list.count]] <- t(combn(1:GA.inputs$n.layers, i))
    if(is.null(nrow(comb.list[[list.count]]))) {
      n.comb <- 1
    } else {
      n.comb <- nrow(comb.list[[list.count]])
    }
    surface.count <- surface.count + n.comb
  }
  
  all.combs <- list()
  comb.names <- list()
  row.index <- 0
  for(i in 1:length(comb.list)){
    combs <- comb.list[[i]]
    
    if(is.null(nrow(comb.list[[i]]))) {
      t.combs <- 1
    } else {
      t.combs <- nrow(comb.list[[i]])
    }
    
    for(j in 1:t.combs) {
      row.index <- row.index + 1
      all.combs[[row.index]] <- combs[j,]
      c.names <- GA.inputs$layer.names[combs[j,]]
      comb.names[[row.index]] <- paste(c.names, collapse = ".")
    }
  }
  
  #all.combs is a list of combinations (models) to test
  v<-c(seq(1:length(all.combs)))
  write_lines(v, "ms.batch")
}


all_comb_pipeline_MS <- function(gdist.inputs,
                     GA.inputs,
                     results.dir,
                     max.combination = NULL,
                     iters = 1000,
                     replicate = 1,
                     sample.prop = 0.75,
                     nlm = FALSE,
                     dist_mod = TRUE,
                     null_mod = TRUE, 
                     cores = 1,
                     my_job_index = NULL) {
  
  if(length(max.combination) > 2) {
    return(cat("ERROR: Please specify either a single value or a vector with the minimum and maximum value"))
  }
  
  # Create combination list -------------------------------------------------
  max.combination<-as.numeric(max.combination)
  my_job_index <- as.numeric(my_job_index)
  mc <- max.combination
  
  if(length(max.combination) == 2) {
    if(mc[1] == 1) {
      min.combination <- 2
      max.combination <- mc[2]
      ss <- TRUE
    } else {
      min.combination <- mc[1]
      max.combination <- mc[2]
      ss <- FALSE
    } 
  } else {
    min.combination <- 2
    ss <- TRUE
  }
  if(max.combination > GA.inputs$n.layers) {
    max.combination <- GA.inputs$n.layers
  }
  comb.list <- vector(mode = "list", length = (max.combination - 1))
  
  list.count <- 0
  surface.count <- 0
  for(i in min.combination:max.combination) {
    list.count <- list.count + 1
    comb.list[[list.count]] <- t(combn(1:GA.inputs$n.layers, i))
    if(is.null(nrow(comb.list[[list.count]]))) {
      n.comb <- 1
    } else {
      n.comb <- nrow(comb.list[[list.count]])
    }
    surface.count <- surface.count + n.comb
  }
  
  all.combs <- list()
  comb.names <- list()
  row.index <- 0
  for(i in 1:length(comb.list)){
    combs <- comb.list[[i]]
    
    if(is.null(nrow(comb.list[[i]]))) {
      t.combs <- 1
    } else {
      t.combs <- nrow(comb.list[[i]])
    }
    
    for(j in 1:t.combs) {
      row.index <- row.index + 1
      all.combs[[row.index]] <- combs[j,]
      c.names <- GA.inputs$layer.names[combs[j,]]
      comb.names[[row.index]] <- paste(c.names, collapse = ".")
    }
  }
  
  GA.input_orig <- GA.inputs
  
  Results <- vector(mode = 'list', length = replicate)
  
  #subset all.combs and comb.names for the chunk provided to this process

  # Optimization without scaling --------------------------------------------
  
  # * Multisurface ----------------------------------------------------------
  
  ms.cd <- vector(mode = 'list',
                  length = length(all.combs))
  
  ms.k <- vector(mode = 'list',
                 length = length(all.combs))
  
  AICc.tab_list <- vector(mode = 'list',
                          length = length(all.combs))
  
  ms.results <- vector(mode = "list", length = length(all.combs))
  
  #n_ss.cd <- length(ss.results$cd)
  #all.cd <- c(ss.results$cd,
   #           ms.cd)
  
  #for(j in 1:length(all.combs)) {
  j <- my_job_index
  print(paste0("My job is: ", j))
    results<-paste0("./MS_RESULTS_", comb.names[[j]])
    plots<-(paste0("./MS_RESULTS_", comb.names[[j]],"/Plots"))
    dir.create(results)
    dir.create(plots)
    
    # Select raster surfaces
    r.vec <- 1:GA.input_orig$n.layers
    drop.vec <- r.vec[!(r.vec %in% all.combs[[j]])]
    asc.comb <- dropLayer(GA.input_orig$Resistance.stack, drop.vec)
    
    if(!is.null(GA.input_orig$inputs$select.trans)) {
      s.trans <- GA.input_orig$inputs$select.trans[all.combs[[j]]]
    } else {
      s.trans <- GA.input_orig$inputs$select.trans
    }
    
    if(is.null(GA.input_orig$inputs$scale)) {
      ms.scale <- FALSE
    } else {
      ms.scale <- TRUE
    }
    
    if(sum(GA.input_orig$inputs$scale.surfaces) == 0) {
      sc.surf <- NULL
    } else {
      sc.surf <- GA.input_orig$inputs$scale.surfaces[all.combs[[j]]]
    }
    
    #sc.surf <- GA.input_orig$inputs$scale.surfaces[all.combs[[j]]]
    
    # Update GA.input 
    GA.inputs <- GA.prep(ASCII.dir = asc.comb,
                         Results.dir = results,
                         min.cat = GA.input_orig$inputs$min.cat,
                         max.cat = GA.input_orig$inputs$max.cat,
                         max.cont = GA.input_orig$inputs$max.cont,
                         min.scale = GA.input_orig$inputs$min.scale,
                         max.scale = GA.input_orig$inputs$max.scale,
                         cont.shape = NULL,
                         select.trans = s.trans,
                         method = GA.input_orig$inputs$method,
                         scale = ms.scale,
                         scale.surfaces = sc.surf,
                         k.value = GA.input_orig$inputs$k.value,
                         pop.mult = GA.input_orig$inputs$pop.mult,
                         percent.elite = GA.input_orig$inputs$percent.elite,
                         type = GA.input_orig$inputs$type,
                         pcrossover = GA.input_orig$inputs$pcrossover,
                         pmutation = GA.input_orig$inputs$pmutation,
                         maxiter = as.numeric(GA.input_orig$inputs$maxiter),
                         run = GA.input_orig$inputs$run,
                         keepBest = GA.input_orig$inputs$keepBest,
                         population = GA.input_orig$inputs$population,
                         selection = GA.input_orig$inputs$selection,
                         crossover = GA.input_orig$inputs$crossover,
                         mutation = GA.input_orig$inputs$mutation,
                         pop.size = GA.input_orig$inputs$pop.size,
                         parallel = as.numeric(cores),
                         seed = GA.input_orig$inputs$seed,
                         quiet = GA.input_orig$inputs$quiet,
                         Plots.dir=plots
    )
    
    # Update GA.input directories
    GA.inputs$Plots.dir <- paste0(".")
    
    GA.inputs$Results.dir <- paste0(".")
    
    ms.results <- MS_pipeline_optim(gdist.inputs = gdist.inputs,
                                GA.inputs = GA.inputs, 
                                NAME=comb.names[[j]])
    print(ms.results)
    
    #all.cd[[(j + n_ss.cd)]] <- ms.results[[j]]$cd[[1]]
    
    #ms.k[[j]] <- ms.results[[j]]$k
    
    #AICc.tab_list[[j]] <- ms.results[[j]]$AICc.tab
    saveRDS(ms.results, paste0(comb.names[[j]], ".MS.rds"))

}



ms_results_gather <- function(gdist.inputs,
                                 GA.inputs,
                                 results.dir,
                                 max.combination = NULL,
                                 iters = 1000,
                                 replicate = 1,
                                 sample.prop = 0.75,
                                 nlm = FALSE,
                                 dist_mod = TRUE,
                                 null_mod = TRUE,
                                 ms.results = NULL,
                                 ss.results = NULL) {



  all.cd <- c(ss.results$cd, ms.results$cd[[1]])

  ms.k <- ms.results$k

  AICc.tab_list[[j]] <- ms.results[[j]]$AICc.tab

  # Convert combination lists to data frames
  all.k <- rbind(ss.results$k,
               plyr::ldply(ms.k))

  names(all.cd) <- all.k$surface


  if(is.null(ss.results)) {
    all.AICc <- plyr::ldply(AICc.tab_list)
  } else {
    all.AICc <- rbind(ss.results$AICc,
                      plyr::ldply(AICc.tab_list))
  }

# all.AICc <- all.AICc %>% 
#   dplyr::arrange(., AICc) %>%
#   dplyr::mutate(., delta.AICc = AICc - min(AICc)) %>%
#   dplyr::mutate(., weight = (exp(-0.5 * delta.AICc)) / sum(exp(-0.5 * delta.AICc))) %>%
#   as.data.frame()
# 
# 
# # * Bootstrap -----------------------------------------------------
# 
# obs <- gdist.inputs$n.Pops
# genetic.mat <- matrix(0, obs, obs)
# genetic.mat[lower.tri(genetic.mat)] <- gdist.inputs$response
# 
# boot.results <- Resist.boot(mod.names = all.k[,1],
#                             dist.mat = all.cd,
#                             n.parameters = all.k[,2],
#                             sample.prop = sample.prop,
#                             iters = iters,
#                             obs = obs,
#                             genetic.mat = genetic.mat)
# 
# } # End scaled if-else
# 
# # Write AICc table and Boot Results to replicate results directory
# write.table(x = all.AICc,
#             paste0(results.dir,'rep_',i,"/",
#                    "All_Combinations_Summary.csv"),
#             row.names = F,
#             col.names = T,
#             sep = ",")
# 
# write.table(x = as.data.frame(boot.results),
#             paste0(results.dir,'rep_',i,"/",
#                    "Bootstrap_Results.csv"),
#             row.names = F,
#             col.names = T,
#             sep = ",")
# 
# if(replicate > 1) {
#   Results[[i]] <- list(summary.table = all.AICc,
#                        boot.results = boot.results,
#                        all.k = all.k,
#                        all.cd = all.cd,
#                        genetic.dist_mat = genetic.mat,
#                        ss.results = ss.results,
#                        ms.results = ms.results
#   )
#   names(Results)[i] <- paste0('rep_',i)
#   
# } else {
#   Results <- list(summary.table = all.AICc,
#                   boot.results = boot.results,
#                   all.k = all.k,
#                   all.cd = all.cd,
#                   genetic.dist_mat = genetic.mat,
#                   ss.results = ss.results,
#                   ms.results = ms.results
#   )
# 
# return(Results)
# 
# } # End function

