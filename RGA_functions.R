
## FOR ASSESSING AICc of fitted models, not exported
Resistance.Opt_AICc <-
  function(PARM,
           Resistance,
           CS.inputs = NULL,
           gdist.inputs = NULL,
           GA.inputs,
           Min.Max = 'max',
           iter = NULL,
           quiet = FALSE) {
    t1 <- proc.time()[3]
    
    
    EXPORT.dir <- GA.inputs$Write.dir
    
    r <- Resistance
    if (!is.null(iter)) {
      if (GA.inputs$surface.type[iter] == "cat") {
        PARM <- PARM / min(PARM)
        parm <- PARM
        df <-
          data.frame(id = unique(r), PARM) # Data frame with original raster values and replacement values
        r <- subs(r, df)
        
      } else {
        r <- SCALE(r, 0, 10)
        
        # Set equation for continuous surface
        equation <- floor(PARM[1]) # Parameter can range from 1-9.99
        
        # Read in resistance surface to be optimized
        SHAPE <- (PARM[2])
        Max.SCALE <- (PARM[3])
        
        # Apply specified transformation
        rick.eq <- (equation == 2 ||
                      equation == 4 ||
                      equation == 6 || equation == 8)
        if (rick.eq == TRUE & SHAPE > 5) {
          equation <- 9
        }
        
        if (equation == 1) {
          r <- Inv.Rev.Monomolecular(r, parm = PARM)
          EQ <- "Inverse-Reverse Monomolecular"
          
        } else if (equation == 5) {
          r <- Rev.Monomolecular(r, parm = PARM)
          EQ <- "Reverse Monomolecular"
          
        } else if (equation == 3) {
          r <- Monomolecular(r, parm = PARM)
          EQ <- "Monomolecular"
          
        } else if (equation == 7) {
          r <- Inv.Monomolecular(r, parm = PARM)
          EQ <- "Inverse Monomolecular"
          
        } else if (equation == 8) {
          r <- Inv.Ricker(r, parm = PARM)
          EQ <- "Inverse Ricker"
          
        } else if (equation == 4) {
          r <- Ricker(r, parm = PARM)
          EQ <- "Ricker"
          
        } else if (equation == 6) {
          r <- Rev.Ricker(r, parm = PARM)
          EQ <- "Reverse Ricker"
          
        } else if (equation == 2) {
          r <- Inv.Rev.Ricker(r, parm = PARM)
          EQ <- "Inverse-Reverse Ricker"
          
        } else {
          r <- (r * 0) + 1 #  Distance
          EQ <- "Distance"
        } # End if-else
      } # Close parameter type if-else
    } else {
      r <- SCALE(r, 0, 10)
      
      # Set equation for continuous surface
      equation <- floor(PARM[1]) # Parameter can range from 1-9.99
      
      # Read in resistance surface to be optimized
      SHAPE <- (PARM[2])
      Max.SCALE <- (PARM[3])
      
      # Apply specified transformation
      rick.eq <- (equation == 2 ||
                    equation == 4 || equation == 6 || equation == 8)
      if (rick.eq == TRUE & SHAPE > 5) {
        equation <- 9
      }
      
      if (equation == 1) {
        r <- Inv.Rev.Monomolecular(r, parm = PARM)
        EQ <- "Inverse-Reverse Monomolecular"
        
      } else if (equation == 5) {
        r <- Rev.Monomolecular(r, parm = PARM)
        EQ <- "Reverse Monomolecular"
        
      } else if (equation == 3) {
        r <- Monomolecular(r, parm = PARM)
        EQ <- "Monomolecular"
        
      } else if (equation == 7) {
        r <- Inv.Monomolecular(r, parm = PARM)
        EQ <- "Inverse Monomolecular"
        
      } else if (equation == 8) {
        r <- Inv.Ricker(r, parm = PARM)
        EQ <- "Inverse Ricker"
        
      } else if (equation == 4) {
        r <- Ricker(r, parm = PARM)
        EQ <- "Ricker"
        
      } else if (equation == 6) {
        r <- Rev.Ricker(r, parm = PARM)
        EQ <- "Reverse Ricker"
        
      } else if (equation == 2) {
        r <- Inv.Rev.Ricker(r, parm = PARM)
        EQ <- "Inverse-Reverse Ricker"
        
      } else {
        r <- (r * 0) + 1 #  Distance
        EQ <- "Distance"
      } # End if-else
    }
    File.name <- "resist_surface"
    if (cellStats(r, "max") > 1e6)
      r <-
      SCALE(r, 1, 1e6) # Rescale surface in case resistance are too high
    r <- reclassify(r, c(-Inf, 1e-06, 1e-06, 1e6, Inf, 1e6))
    
    
    if (!is.null(CS.inputs)) {
      writeRaster(
        x = r,
        filename = paste0(EXPORT.dir, File.name, ".asc"),
        overwrite = TRUE
      )
      CS.resist <-
        Run_CS2(
          CS.inputs,
          GA.inputs,
          r = r,
          EXPORT.dir = GA.inputs$Write.dir,
          File.name = File.name
        )
      
      # Replace NA with 0...a workaround for errors when two points fall within the same cell.
      # CS.resist[is.na(CS.resist)] <- 0
      
      # Run mixed effect model on each Circuitscape effective resistance
      AIC.stat <- suppressWarnings(AIC(
        MLPE.lmm2(
          resistance = CS.resist,
          response = CS.inputs$response,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ,
          REML = FALSE
        )
      ))
      ROW <- nrow(CS.inputs$ID)
      
    }
    
    if (!is.null(gdist.inputs)) {
      cd <- Run_gdistance(gdist.inputs, r)
      
      AIC.stat <- suppressWarnings(AIC(
        MLPE.lmm2(
          resistance = cd,
          response = gdist.inputs$response,
          ID = gdist.inputs$ID,
          ZZ = gdist.inputs$ZZ,
          REML = FALSE
        )
      ))
      ROW <- nrow(gdist.inputs$ID)
    }
    
    k <- length(PARM) + 1
    AICc <- (AIC.stat) + (((2 * k) * (k + 1)) / (ROW - k - 1))
    
    rt <- proc.time()[3] - t1
    if (quiet == FALSE) {
      cat(paste0("\t", "Iteration took ", round(rt, digits = 2), " seconds"), "\n")
      #     cat(paste0("\t", EQ,"; ",round(SHAPE,digits=2),"; ", round(Max.SCALE,digits=2)),"\n")
      cat(paste0("\t", "AICc = ", round(AICc, 4)), "\n")
      if (!is.null(iter)) {
        if (GA.inputs$surface.type[iter] != "cat") {
          cat(paste0("\t", EQ, " | Shape = ", PARM[2], " | Max = ", PARM[3]),
              "\n",
              "\n")
        }
      }
    }
    OPTIM.DIRECTION(Min.Max) * (AICc) # Function to be minimized/maximized
  }

# FUNCTIONS
OPTIM.DIRECTION <- function(x) {
  OPTIM <- ifelse(x == 'max', -1, 1)
  return(OPTIM)
}


Cont.Param <- function(PARM) {
  df <- data.frame(PARM[1], PARM[2])
  colnames(df) <- c("shape_opt", "max")
  row.names(df) <- NULL
  return(df)
}


read.matrix <-
  function(cs.matrix) {
    lower(read.table(cs.matrix)[-1, -1])
  }


read.matrix2 <- function(cs.matrix) {
  m <- read.table(cs.matrix)[-1, -1]
}


# Create ZZ matrix for mixed effects model
ZZ.mat <- function(ID) {
  Zl <-
    lapply(c("pop1", "pop2"), function(nm)
      Matrix::fac2sparse(ID[[nm]], "d", drop = FALSE))
  ZZ <- Reduce("+", Zl[-1], Zl[[1]])
  return(ZZ)
}

# Rescale function
SCALE.vector <- function(data, MIN, MAX, threshold = 1e-5) {
  if (abs(MIN - MAX) < threshold) {
    data[is.finite(data)] <- 0
    data
  } else {
    Mn = min(data)
    Mx = max(data)
    (MAX - MIN) / (Mx - Mn) * (data - Mx) + MAX
  }
}

# Define scaling function
# This will rescale from 1 to specified MAX
SCALE <- function(data, MIN, MAX, threshold = 1e-5) {
  if (abs(MIN - MAX) < threshold) {
    data[is.finite(raster::values(data))] <- 0
    data
  } else {
    Mn = cellStats(data, stat = 'min')
    Mx = cellStats(data, stat = 'max')
    (MAX - MIN) / (Mx - Mn) * (data - Mx) + MAX
  }
}

# Sample values for suggests
sv.cat <- function(levels, pop.size, min, max) {
  cat.starts <- matrix(nrow = pop.size, ncol = levels)
  for (r in 1:pop.size) {
    L <- list()
    for (i in 1:levels) {
      if (runif(1) < .5) {
        z <- runif(1)
      } else {
        z <- runif(1, min, max)
      }
      L[[i]] <- z
    }
    #   uz<-unlist(L)
    cat.starts[r, ] <- (unlist(L))
  }
  cat.starts[, 1] <- 1
  return(cat.starts)
}


# No Gaussian distribution
sv.cont.nG <- function(direction,
                       pop.size,
                       max,
                       min.scale,
                       max.scale,
                       scale = NULL, 
                       eqs) {
  inc <- c(1, 3)
  dec <- c(7, 5)
  peak <- c(2, 4, 6, 8)
  L <- list()
  
  if (!is.null(scale)) {
    cont.starts <- matrix(nrow = pop.size, ncol = 4)
    for (r in 1:pop.size) {
      scale.parm <- runif(1, min.scale, max.scale)
      if (runif(1) < .5 && direction == "Increase") {
        #       z1<-c(sample(inc,1)
        z <- Increase.starts.nG(sample(inc, 1))
        z <- c(z, scale.parm)
      } else if (runif(1) < .5 && direction == "Decrease") {
        z <- c(sample(dec, 1),
               runif(1, .01, 10),
               runif(1, 1, max),
               scale.parm)
      } else if (runif(1) < .5 && direction == "Peaked") {
        z <- c(sample(peak, 1),
               runif(1, .01, 10),
               runif(1, 1, max),
               scale.parm)
      } else {
        z <- c(sample(eqs, 1),
               # runif(1, 1, 9.99),
               runif(1, .01, 10),
               runif(1, 1, max),
               scale.parm)
      }
      cont.starts[r,] <- z
    }
  } else {
    cont.starts <- matrix(nrow = pop.size, ncol = 3)
    for (r in 1:pop.size) {
      if (runif(1) < .5 && direction == "Increase") {
        #       z1<-c(sample(inc,1)
        z <- Increase.starts.nG(sample(inc, 1))
      } else if (runif(1) < .5 && direction == "Decrease") {
        z <- c(sample(dec, 1), runif(1, .01, 10), runif(1, 1, max))
      } else if (runif(1) < .5 && direction == "Peaked") {
        z <- c(sample(peak, 1), runif(1, .01, 10), runif(1, 1, max))
      } else {
        z <- c(sample(eqs, 1),
               # runif(1, 1, 9.99), 
               runif(1, .01, 10), 
               runif(1, 1, max))
      }
      cont.starts[r,] <- z
    }
  }
  if(ncol(cont.starts) == 4) {
    rs <- sample(pop.size, floor(0.25 * pop.size), replace = F)
    cont.starts[rs, 4] <- 0.25
  }
  
  cont.starts
}

Increase.starts <- function(x) {
  if (x == 1) {
    z <- c(x, runif(1, .01, 10), runif(1, 1, 10), 1)
  } else {
    z <- c(x, runif(1, .01, 10), runif(1, 1, 100), 1)
  }
}

Increase.starts.nG <- function(x) {
  if (x == 1) {
    z <- c(x, runif(1, .01, 10), runif(1, 1, 10))
  } else {
    z <- c(x, runif(1, .01, 10), runif(1, 1, 100))
  }
}



unique <- raster::unique

eq.set <- function(include.list) {
  out <- vector(mode = "list", length = length(include.list))
  for (i in seq_along(include.list)) {
    if (include.list[[i]] == "A" | is.na(include.list[[i]])) {
      out <- 1:9
      return(out)
    } else if (include.list[[i]] == "M") {
      out <- c(1, 3, 5, 7, 9)
      return(out)
    } else if (include.list[[i]] == "R") {
      out <- c(2, 4, 6, 8, 9)
      return(out)
    } else if (!is.na(match(include.list[[i]], 1:9))) {
      out <- include.list
      return(out)
    } else {
      cat(
        "The specified transformations to assess are not valid. \n
        Please see Details of the GA.prep."
      )
    }
    }
  }

get.EQ <- function(equation) {
  # Apply specified transformation
  if (is.numeric(equation)) {
    equation = floor(equation)
    if (equation == 1) {
      EQ <- "Inverse-Reverse Monomolecular"
      
    } else if (equation == 5) {
      EQ <- "Reverse Monomolecular"
      
    } else if (equation == 3) {
      EQ <- "Monomolecular"
      
    } else if (equation == 7) {
      EQ <- "Inverse Monomolecular"
      
    } else if (equation == 8) {
      EQ <- "Inverse Ricker"
      
    } else if (equation == 4) {
      EQ <- "Ricker"
      
    } else if (equation == 6) {
      EQ <- "Reverse Ricker"
      
    } else if (equation == 2) {
      EQ <- "Inverse-Reverse Ricker"
      
    } else {
      EQ <- "Distance"
    }
    
    (EQ)
  } else {
    if (equation == "Inverse-Reverse Monomolecular") {
      EQ <- 1
      
    } else if (equation == "Reverse Monomolecular") {
      EQ <- 5
      
    } else if (equation == "Monomolecular") {
      EQ <- 3
      
    } else if (equation == "Inverse Monomolecular") {
      EQ <- 7
      
    } else if (equation == "Inverse Ricker") {
      EQ <- 8
      
    } else if (equation == "Ricker") {
      EQ <- 4
      
    } else if (equation == "Reverse Ricker") {
      EQ <- 6
      
    } else if (equation == "Inverse-Reverse Ricker") {
      EQ <- 2
      
    } else {
      EQ <- 9
    }
    
    (EQ)
  }
}

Result.txt <-
  function(GA.results,
           GA.inputs,
           method,
           k,
           Run.Time,
           fit.stats,
           MLPE.coef = NULL,
           optim,
           aic,
           AICc,
           LL) {
    summary.file <-
      paste0(GA.inputs$Results.dir, "Multisurface_Optim_Summary.txt")
    # AICc<-GA.results@fitnessValue
    # AICc<-round(AICc,digits=4)
    ELITE <- floor(GA.inputs$percent.elite * GA.inputs$pop.size)
    #   mlpe.results<-MLPE.lmm_coef(GA.inputs$Results.dir,genetic.dist=CS.inputs$response)
    
    sink(summary.file)
    cat(paste0(
      "Summary from multisurface optimization run conducted on ",
      Sys.Date()
    ),
    "\n")
    cat("\n")
    cat(paste0("Optimized using: ", method), "\n")
    cat("\n")
    cat(paste0("Objective function: ", optim), "\n")
    cat("\n")
    cat(paste0("Surfaces included in optimization:"), "\n")
    cat(GA.inputs$parm.type$name, "\n")
    cat("\n")
    cat("Genetic Algorithm optimization settings:")
    cat("\n")
    cat(paste0("Type of genetic algorithm used: ", GA.results@type),
        "\n")
    cat(paste0("popSize at each iteration: ", GA.results@popSize),
        "\n")
    cat(paste0("Maximum number of iterations: ", GA.results@maxiter),
        "\n")
    cat(paste0(
      "Number individuals retained each generation (elistism): ",
      ELITE
    ),
    "\n")
    cat(paste0("Crossover probability: ", GA.results@pcrossover),
        "\n")
    cat(paste0("Mutation probability: ", GA.results@pmutation), "\n")
    cat("\n")
    cat(
      paste0(
        "The Genetic Algorithm completed after ",
        GA.results@iter,
        " iterations"
      ),
      "\n"
    )
    cat("\n")
    cat(paste0("k =  ", k), "\n")
    cat("\n")
    cat(paste0("Minimum AIC: ", aic), "\n")
    cat("\n")
    cat(paste0("AICc: ", AICc), "\n")
    cat("\n")
    cat(paste0("Pseudo marginal R-square (R2m): ", fit.stats[[1]]), "\n")
    cat(paste0("Pseudo conditional R-square (R2c): ", fit.stats[[2]]),
        "\n")
    cat("\n")
    cat(paste0("Log Likelihood: ", LL), "\n")
    cat("\n")
    cat(paste0("Optimized values for each surface:"), "\n")
    cat(GA.results@solution, "\n")
    cat("\n")
    cat(paste0("Optimization took ", Run.Time, " seconds to complete"),
        "\n")
    sink()
  }

##########################################################################################

get.name <- function(x) {
  nm <- deparse(substitute(x))
  return(nm)
}


# Transformation Eqs ------------------------------------------------------

Monomolecular <- function(r, parm) {
  parm[3] * (1 - exp(-1 * r / parm[2])) + 1 # Monomolecular
}

Inv.Monomolecular <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    R <- parm[3] * (exp(-1 * r / parm[2]))
    (R <- (R - cellStats(R, stat = "min")) + 1)
  } else {
    R <- parm[3] * (exp(-1 * r / parm[2]))
    (R <- (R - min(R)) + 1)
  }
}

Inv.Rev.Monomolecular <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Inv.Monomolecular(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Inv.Monomolecular(rev.rast, parm)
  }
}

Rev.Monomolecular <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Monomolecular(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Monomolecular(rev.rast, parm)
  }
}


Ricker <- function(r, parm) {
  parm[3] * r * exp(-1 * r / parm[2]) + 1 # Ricker
}

Inv.Ricker <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    R <- (-1 * parm[3]) * r * exp(-1 * r / parm[2]) - 1 # Ricker
    R <-
      SCALE(R,
            MIN = abs(cellStats(R, stat = 'max')),
            MAX = abs(cellStats(R, stat = 'min'))) # Rescale
  } else {
    R <- (-1 * parm[3]) * r * exp(-1 * r / parm[2]) - 1 # Ricker
    R <- SCALE.vector(R, MIN = abs(max(R)), MAX = abs(min(R))) # Rescale
  }
}

Inv.Rev.Ricker <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Inv.Ricker(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Inv.Ricker(rev.rast, parm)
  }
}

Rev.Ricker <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Ricker(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Ricker(rev.rast, parm)
  }
}

yn.question <- function(question, add_lines_before = TRUE) {
  choices <- c("Yes", "No", "New Subdirectory")
  if(add_lines_before) cat("------------------------\n")   
  the_answer <- menu(choices, title = question)
  
  if(the_answer == 1L) {
    return(TRUE)
  } else if(the_answer == 2L) {
    return(FALSE)
  } else if(the_answer == 3L){
    return(NA)
  }
  
  # ifelse(the_answer == 1L, TRUE, FALSE)   # returns TRUE or FALSE
}

SS_optim_worker <- function(CS.inputs = NULL,
                            gdist.inputs = NULL,
                            GA.inputs,
                            nlm = FALSE,
                            dist_mod = TRUE,
                            null_mod = TRUE,
                            myproc = NULL, 
                            jobs = NULL) {
  if (!is.null(GA.inputs$scale)) {
    stop(
      "This function should NOT be used if you intend to apply kernel smoothing to your resistance surfaces"
    )
  }
  t1 <- proc.time()[3]
  RESULTS.cat <- list() # List to store categorical results within
  RESULTS.cont <- list() # List to store continuous results within
  cnt1 <- 0
  cnt2 <- 0
  k.value <- GA.inputs$k.value
  MLPE.list <- list()
  cd.list <- list()
  k.list <- list()
  
  # Optimize each surface in turn
  for (i in 1:GA.inputs$n.layers) {
    #if i not in list of job indices, skip
    if (!i %in% jobs){
      next
    }else{
      print(paste0("Proc ", myproc, " performing SS optimization on ", GA.inputs$layer.names[i]))
      r <- GA.inputs$Resistance.stack[[i]]
      names(r) <- GA.inputs$layer.names[i]
      
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
          
          single.GA@solution <-
            single.GA@solution / min(single.GA@solution)
          df <- data.frame(id = unique(r), t(single.GA@solution))
          r <- subs(r, df)
          names(r) <- GA.inputs$layer.names[i]
          
          Run_CS(CS.inputs, GA.inputs, r, EXPORT.dir = GA.inputs$Results.dir)
          
          Diagnostic.Plots(
            resistance.mat = paste0(
              GA.inputs$Results.dir,
              GA.inputs$layer.names[i],
              "_resistances.out"
            ),
            genetic.dist = CS.inputs$response,
            plot.dir = GA.inputs$Plots.dir,
            type = "categorical",
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
          
          fit.stats <-
            r.squaredGLMM(
              MLPE.lmm(
                resistance = paste0(
                  GA.inputs$Results.dir,
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
                  GA.inputs$Results.dir,
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
                  GA.inputs$Results.dir,
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
          
          RESULTS.cat[[cnt1]] <- RS
          
          MLPE.list[[i]] <- MLPE.lmm(
            resistance = paste0(
              GA.inputs$Results.dir,
              GA.inputs$layer.names[i],
              "_resistances.out"
            ),
            pairwise.genetic = CS.inputs$response,
            REML = TRUE,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
          
          cd.list[[i]] <- (read.table(paste0(
            GA.inputs$Results.dir,
            GA.inputs$layer.names[i],
            "_resistances.out"))[-1, -1])
          
          names(MLPE.list)[i] <- GA.inputs$layer.names[i]
          names(cd.list)[i] <- GA.inputs$layer.names[i]
          
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
                write.dir = GA.inputs$Write.dir
              )
            
            OPTIM <-
              Resistance.Optimization_cont.nlm(
                PARM = (Optim.nlm$estimate),
                Resistance = r,
                equation = single.GA@solution[1],
                get.best = TRUE,
                CS.inputs = CS.inputs,
                Min.Max = 'min',
                write.dir = GA.inputs$Results.dir
              )
            
            Diagnostic.Plots(
              resistance.mat = paste0(
                GA.inputs$Results.dir,
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
            RESULTS.cont[[cnt2]] <- RS
            
            MLPE.list[[i]] <- MLPE.lmm(
              resistance = paste0(
                GA.inputs$Results.dir,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              pairwise.genetic = CS.inputs$response,
              REML = TRUE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
            
            cd.list[[i]] <- (read.table(paste0(
              GA.inputs$Results.dir,
              GA.inputs$layer.names[i],
              "_resistances.out"))[-1, -1])
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
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
                   EXPORT.dir = GA.inputs$Results.dir)
            
            Diagnostic.Plots(
              resistance.mat = paste0(
                GA.inputs$Results.dir,
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
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
            fit.stats <- r.squaredGLMM(
              MLPE.lmm(
                REML = F,
                resistance = paste0(
                  GA.inputs$Results.dir,
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
                  GA.inputs$Results.dir,
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
                    GA.inputs$Results.dir,
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
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
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
            RESULTS.cont[[cnt2]] <- RS
            
            MLPE.list[[i]] <- MLPE.lmm(
              resistance = paste0(
                GA.inputs$Results.dir,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              pairwise.genetic = CS.inputs$response,
              REML = TRUE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
            
            cd.list[[i]] <- (read.table(paste0(
              GA.inputs$Results.dir,
              GA.inputs$layer.names[i],
              "_resistances.out"))[-1, -1])
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
          }
        } # Close if-else
        if (dist_mod == TRUE) {
          r <- reclassify(r, c(-Inf, Inf, 1))
          names(r) <- "dist"
          cd <- Run_CS(CS.inputs, GA.inputs, r, full.mat = T)
          
          Dist.AIC <-
            AIC(
              MLPE.lmm(
                resistance = paste0(GA.inputs$Write.dir, "dist_resistances.out"),
                pairwise.genetic = CS.inputs$response,
                REML = FALSE,
                ID = CS.inputs$ID,
                ZZ = CS.inputs$ZZ
              )
            )
          
          fit.stats <-
            r.squaredGLMM(
              MLPE.lmm(
                resistance = paste0(GA.inputs$Write.dir, "dist_resistances.out"),
                pairwise.genetic = CS.inputs$response,
                REML = FALSE,
                ID = CS.inputs$ID,
                ZZ = CS.inputs$ZZ
              )
            )
          
          LL <-
            logLik(
              MLPE.lmm(
                resistance = paste0(GA.inputs$Write.dir, "dist_resistances.out"),
                pairwise.genetic = CS.inputs$response,
                REML = FALSE,
                ID = CS.inputs$ID,
                ZZ = CS.inputs$ZZ
              )
            )
          
          MLPE.list[[i + 1]] <- MLPE.lmm(
            resistance = paste0(GA.inputs$Write.dir, "dist_resistances.out"),
            pairwise.genetic = CS.inputs$response,
            REML = TRUE,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
          
          cd.list[[i + 1]] <- cd
          # (read.table(paste0(GA.inputs$Write.dir, "dist_resistances.out"))[-1, -1])
          
          names(cd.list)[i + 1] <- 'Distance'
          
          names(MLPE.list)[i + 1] <- "Distance"
          
          if (GA.inputs$method == "AIC") {
            dist.obj <- Dist.AIC
          } else if (GA.inputs$method == "R2") {
            dist.obj <- fit.stats[[1]]
          } else {
            dist.obj <- LL[[1]]
          }
          
          k <- 2
          k.list[[i + 1]] <- k
          names(k.list)[i + 1] <- 'Distance'
          
          n <- CS.inputs$n.Pops
          AICc <-
            (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
          # (Dist.AIC)+(((2*k)*(k+1))/((CS.inputs$n.Pops)-k-1))
          
          Dist.AICc <-
            data.frame("Distance",
                       dist.obj,
                       k,
                       Dist.AIC,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Dist.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              "k",
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
        
        if (null_mod == TRUE) {
          response = CS.inputs$response
          
          dat <- data.frame(CS.inputs$ID, response = CS.inputs$response)
          colnames(dat) <- c("pop1", "pop2", "response")
          
          # Fit model
          mod <-
            lFormula(response ~ 1 + (1 | pop1),
                     data = dat,
                     REML = FALSE)
          mod$reTrms$Zt <- CS.inputs$ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          Null.AIC <-
            AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          fit.stats <-
            r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          LL <-
            logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          
          if (GA.inputs$method == "AIC") {
            null.obj <- Null.AIC
          } else if (GA.inputs$method == "R2") {
            null.obj <- fit.stats[[1]]
          } else {
            null.obj <- LL[[1]]
          }
          k <- 1
          n <- CS.inputs$n.Pops
          AICc <-
            (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
          # AICc <- (Null.AIC)+(((2*k)*(k+1))/((CS.inputs$n.Pops)-k-1))
          Null.AICc <-
            data.frame("Null",
                       null.obj,
                       k,
                       Null.AIC,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Null.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              "k",
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
        
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
            file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
            
            sep = ",",
            row.names = F,
            col.names = F
          )
          writeRaster(r,
                      paste0(GA.inputs$Results.dir, NAME, ".asc"),
                      overwrite = TRUE)
          
          # save(single.GA,
          #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
          
          saveRDS(single.GA,
                  file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
          
          Diagnostic.Plots(
            resistance.mat = cd,
            genetic.dist = gdist.inputs$response,
            plot.dir = GA.inputs$Plots.dir,
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
          
          k.list[[i]] <- k
          names(k.list)[i] <- GA.inputs$layer.names[i]
          
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
          
          RESULTS.cat[[cnt1]] <- RS
          
          MLPE.list[[i]] <-  MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            REML = TRUE,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ
          )
          
          cd.list[[i]] <- as.matrix(cd)
          names(cd.list)[i] <- GA.inputs$layer.names[i]
          
          names(MLPE.list)[i] <- GA.inputs$layer.names[i]
          
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
                write.dir = GA.inputs$Write.dir
              )
            
            r <-
              Resistance.tran(
                transformation = EQ,
                shape = Optim.nlm$estimate[1],
                max = Optim.nlm$estimate[2],
                r = r,
                out = GA.inputs$Results.dir
              )
            
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA,
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA,
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = exp(Optim.nlm$estimate),
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir,
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
            RESULTS.cont[[cnt2]] <- RS
            
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
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA,
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA,
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
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
            
            MLPE.list[[i]] <-  MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = TRUE,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
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
            RESULTS.cont[[cnt2]] <- RS
          }
        } # Close if-else
        
        if (dist_mod == TRUE) {
          r <- reclassify(r, c(-Inf, Inf, 1))
          names(r) <- "dist"
          cd <- Run_gdistance(gdist.inputs, r)
          
          write.table(
            as.matrix(cd),
            file = paste0(GA.inputs$Results.dir, 'Distance', "_", gdist.inputs$method, "_distMat.csv"),
            sep = ",",
            row.names = F,
            col.names = F
          )
          
          Dist.AIC <- suppressWarnings(AIC(
            MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ,
              REML = FALSE
            )
          ))
          
          fit.stats <- r.squaredGLMM(
            MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ,
              REML = FALSE
            )
          )
          
          LL <- logLik(
            MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ,
              REML = FALSE
            )
          )
          
          MLPE.list[[i + 1]] <-  MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            REML = TRUE,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ
          )
          
          cd.list[[i + 1]] <- as.matrix(cd)
          names(cd.list)[i + 1] <- "Distance"
          
          names(MLPE.list)[i + 1] <- "Distance"
          
          ROW <- nrow(gdist.inputs$ID)
          k <- 2
          
          k.list[[i + 1]] <- k
          names(k.list)[i + 1] <- 'Distance'
          
          if (GA.inputs$method == "AIC") {
            dist.obj <- -Dist.AIC
          } else if (GA.inputs$method == "R2") {
            dist.obj <- fit.stats[[1]]
          } else {
            dist.obj <- LL[[1]]
          }
          
          n <- gdist.inputs$n.Pops
          AICc <-
            (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
          # AICc <- (Dist.AIC)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
          Dist.AICc <- data.frame("Distance",
                                  dist.obj,
                                  k,
                                  Dist.AIC,
                                  AICc,
                                  fit.stats[[1]],
                                  fit.stats[[2]],
                                  LL[[1]])
          colnames(Dist.AICc) <- c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            'k',
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL"
          )
        }
        
        if (null_mod == TRUE) {
          dat <- data.frame(gdist.inputs$ID, response = gdist.inputs$response)
          colnames(dat) <- c("pop1", "pop2", "response")
          
          # Fit model
          mod <-
            lFormula(response ~ 1 + (1 | pop1),
                     data = dat,
                     REML = FALSE)
          mod$reTrms$Zt <- gdist.inputs$ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          Null.AIC <-
            AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          fit.stats <-
            r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          LL <-
            logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
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
          AICc <-
            (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
          # AICc <- (Null.AIC)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
          Null.AICc <-
            data.frame("Null",
                       null.obj,
                       k,
                       Null.AIC,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Null.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
      }
    }
    #write Rdata outputs
    saveRDS(RESULTS.cat, file=paste0(myproc + "_results_cat.Rdata"), compress=TRUE)
    saveRDS(RESULTS.cont, file=paste0(myproc + "_results_cont.Rdata"), compress=TRUE)
    saveRDS(MLPE.list, file=paste0(myproc + "_results_MLPE-list.Rdata"), compress=TRUE)
    saveRDS(cd.list, file=paste0(myproc + "_results_cd-list.Rdata"), compress=TRUE)
    saveRDS(k.list, file=paste0(myproc + "_results_k-list.Rdata"), compress=TRUE)
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
                            RUN_NAME = "Run1") {
  
  if(!exists('results.dir'))
    return(cat("ERROR: An empty results directory must be specified"))
  
  if(!exists('gdist.inputs'))
    return(cat("ERROR: Please specify gdist.inputs"))
  
  if(!exists('GA.inputs'))
    return(cat("ERROR: Please specify GA.inputs"))
  
  if(!is.null(GA.inputs$Results.dir) &
     !is.null(GA.inputs$Write.dir) &
     !is.null(GA.inputs$Plots.dir)) {
    return(cat("ERROR: Please correctly specify the `Results.dir` as 'all.comb' when running GA.prep"))
  }
  
  if(length(max.combination) > 2) {
    return(cat("ERROR: Please specify either a single value or a vector with the minimum and maximum value"))
  }
  
  dir.files <- list.files(results.dir)
  if(length(dir.files) != 0) {
    q <- TRUE
    if(q) { # Create subdir
      dir.NAME <- RUN_NAME
      dir.create(path = paste0(results.dir, "all_comb_", dir.NAME))
      results.dir <- paste0(results.dir, "all_comb_", dir.NAME, "/")
    }
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
  
  # Make results data frame
  Results.cat <- data.frame()
  Results.cont <- data.frame()
  # cnt1<-0
  # cnt2<-0
  for (i in 1:GA.inputs$n.layers) {
    if (GA.inputs$surface.type[i] == 'cat') {
      #     cnt1 <- cnt1+1
      #     RS <- data.frame(GA.inputs$layer.names[i], -(RESULTS.cat[[i]]@fitnessValue),RESULTS[[i]]@solution)
      Results.cat <- do.call(rbind.fill, RESULTS.cat)
    } else {
      #   cnt2 <-cnt2+1
      #   RS <- data.frame(GA.inputs$layer.names[i], -(RESULTS.cont[[i]]@fitnessValue), Cont.Param(RESULTS[[i]]@solution))
      Results.cont <- do.call(rbind, RESULTS.cont)
    }
  }
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
      paste0(GA.inputs$Results.dir, "CategoricalResults.csv"),
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
      paste0(GA.inputs$Results.dir, "ContinuousResults.csv"),
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
  
  if (dist_mod == TRUE)
    Results.All <- rbind(Results.All, Dist.AICc)
  if (null_mod == TRUE)
    Results.All <- rbind(Results.All, Null.AICc)
  
  Results.All <- Results.All[order(Results.All$AICc), ]
  
  cat("\n")
  cat("\n")
  write.table(
    Results.All,
    paste0(GA.inputs$Results.dir, "All_Results_Table_", gdist.inputs$method,".csv"),
    
    sep = ",",
    col.names = T,
    row.names = F
  )
  
  # Get parameter estimates
  if (!is.null(CS.inputs)) {
    MLPE.results <- MLPE.lmm_coef(
      resistance = GA.inputs$Results.dir,
      genetic.dist = CS.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "cs",
      ID = CS.inputs$ID,
      ZZ = CS.inputs$ZZ
    )
    
  } else {
    MLPE.results <- MLPE.lmm_coef(
      resistance = GA.inputs$Results.dir,
      genetic.dist = gdist.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "gd",
      ID = gdist.inputs$ID,
      ZZ = gdist.inputs$ZZ
    )
  }
  
  k.list <- plyr::ldply(k.list)
  colnames(k.list) <- c("surface", "k")
  
  rt <- proc.time()[3] - t1
  # Full Results
  if (nrow(Results.cat) > 0 & nrow(Results.cont) > 0) {
    RESULTS <-
      list(
        ContinuousResults = Results.cont,
        CategoricalResults = Results.cat,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
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
                                 RUN_NAME = "Run1") {
  
  if(!exists('results.dir')) 
    return(cat("ERROR: An empty results directory must be specified"))
  
  if(!exists('gdist.inputs')) 
    return(cat("ERROR: Please specify gdist.inputs"))
  
  if(!exists('GA.inputs')) 
    return(cat("ERROR: Please specify GA.inputs"))
  
  if(!is.null(GA.inputs$Results.dir) & 
     !is.null(GA.inputs$Write.dir) &
     !is.null(GA.inputs$Plots.dir)) {
    return(cat("ERROR: Please correctly specify the `Results.dir` as 'all.comb' when running GA.prep"))
  }
  
  if(length(max.combination) > 2) {
    return(cat("ERROR: Please specify either a single value or a vector with the minimum and maximum value"))
  }
  
  dir.files <- list.files(results.dir)
  if(length(dir.files) != 0) {
    q <- TRUE
    
    if(q) { # Create subdir
      dir.NAME <- RUN_NAME
      dir.create(path = paste0(results.dir, "all_comb_", dir.NAME))
      results.dir <- paste0(results.dir, "all_comb_", dir.NAME, "/")
    } 
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
  
  # Begin Replicate Loop --------------------------------------------------
  num_jobs <- replicate * GA.inputs$n.layers
  if (num_jobs < procs){
    procs = num_jobs
    floorjobs = 1
    remainsjobs = 0
  }else{
    floorjobs <- floor(num_jobs / procs)
    remainjobs <- num_jobs %% procs #last proc has to take these too
  }
  if (procnum > procs){
    exit(paste0("Proc number ", procnum, " not needed (only ", numjobs, " and ", procs, " available processes)"))
  }
  
  assignments <- list()
  indices <- list()
  jcount<-1
  pcount<-1
  for(i in 1:replicate){
    temp<-vector()
    itemp<-vector()
    for( j in 1:GA.inputs$n.layers){
      itemp[[j]] <- j
      if (pcount == procs){
        #print(1)
        temp[[j]] <- pcount
      }else if (jcount == floorjobs){
        #print(2)
        temp[[j]] = pcount
        jcount = 1
        pcount= pcount + 1
      }else{
        #print(3)
        temp[[j]] <- pcount
        jcount = jcount + 1
      }
    }
    assignments[[i]]<-temp
    indices[[i]]<-itemp
  }
  print(assignments)
  print(indices)
  
  print(paste0("My process number is " , procnum , " of " , procs , " processes!\n"))
  if (procnum != procs){
    print(paste0("I'm responsible for " , toString(floorjobs) , " single-surface optimizations!\n"))
  }else{
    print(paste0("I'm responsible for " , toString(floorjobs+remainjobs) , " single-surface optimizations!\n"))
  }
  
  for(i in 1:replicate){
    assigns <- assignments[[i]]
    ind       <- indices[[i]]
    #skip if my proc has no jobs in this replicate
    if (!procnum %in% assigns){ next }
    #print(assigns)
    myname <- paste0("rep",i,"_",procnum)
    myjobs <- ind[assigns == procnum]
    #print(myname)
    print(paste(c("My job indices in Replicate", i, ": ", myjobs), collapse=" "))
    
    # Skip if min combination > 1
    if(ss == FALSE) {
      ss.results <- NULL
      AICc.tab <- NULL
      dir.create(paste0(results.dir,'rep_',i))
    } else {  # Do single surface optimization
      dir.create(paste0(results.dir,'rep_',i))
      dir.create(paste0(results.dir,'rep_',i, "/", "single.surface"))
      dir.create(paste0(results.dir,'rep_',i, "/", "single.surface/", "Plots"))
      
      # Single Surface optimization -----------------------------------------------------
      if(!is.null(GA.inputs$scale)) {
        stop(
          "all_comb_pipeline_SS  NOT WORKING YET FOR SCALED SURFACE OPTIMIZATION. EMAIL tkchafin@uark.edu FOR MORE INFORMATION"
        )
        # * Single Surface: scaled --------------------------------------------------------
        
        # Update GA.input directories
        GA.inputs$Plots.dir <- paste0(results.dir,
                                      'rep_',i,
                                      "/",
                                      "single.surface/",
                                      "Plots/")
        
        GA.inputs$Results.dir <- paste0(results.dir,
                                        'rep_',i,
                                        "/",
                                        "single.surface/")
        
        ss.results <- SS_optim.scale(gdist.inputs = gdist.inputs,
                                     GA.inputs = GA.inputs,
                                     nlm = nlm,
                                     dist_mod = dist_mod,
                                     null_mod = null_mod)
        
        AICc.tab <- ss.results$AICc
      } else {
        # * Single Surface --------------------------------------------------------
        
        # Update GA.input directories
        GA.inputs$Plots.dir <- paste0(results.dir,
                                      'rep_',i,
                                      "/",
                                      "single.surface/",
                                      "Plots/")
        
        GA.inputs$Results.dir <- paste0(results.dir,
                                        'rep_',i,
                                        "/",
                                        "single.surface/")
        
        ss.results <- SS_optim_worker(gdist.inputs = gdist.inputs,
                                      GA.inputs = GA.inputs,
                                      nlm = nlm,
                                      dist_mod = dist_mod,
                                      null_mod = null_mod,
                                      jobs = myjobs,
                                      myproc = myname)
        #NEXT: Call SS_optim_gather
        #THEN: Continue pipeline for MS optimizations
        AICc.tab <- ss.results$AICc
      }
    }
  } # End function
}
