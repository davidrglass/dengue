###################################################################################
#
# Script: batch_normalize.R
# Project: Dengue
# Author: David Glass
# Date: 8-14-20
#
# This script batch normalizes fcs files from different barcoded plates. Input and
# output are both debarcoded fcs files from 8 different runs.
#
# Pseudocode:
#   Read in control samples
#   Asinh transform
#   Visualize
#   Identify median of each marker of each control sample
#   Subtract downward so every sample has the same mode as the lowest sample
#   Convert negative numbers to zero + random normal distribution
#   Identify percentiles of each control sample
#   Percentile normalize
#   Visualize
#  
#  Read in samples
#  Asinh transform
#  Visualize
#  Subtract downward
#  Convert negative numbers to zero
#  Percentile normalize
#  Reverse transform
#  Write out samples
#
# Instructions:
# Update USER INPUTS with directory locations
# Run panel editor GUI to remove beadDist channel and place
#   resulting fcs files into input.fcs.path
# Run script
#
# Versions:
# R 3.6.3
# Rstudio 1.3.1093
# premessa 0.2.6
# flowCore 1.52.1
# tidyverse 1.3.1
# data.table 1.14.0
# 
##################################################################################



###### LIBRARIES ######

require(premessa)
require(flowCore)
require(tidyverse)
require(data.table)



##### INPUTS #####

# directory with fcs files to read in
input.fcs.path <- "~/Research/dengue/fcs/pre_batch_correction/bead_dist_removed/"
# directory to write and read tables
tables.path <- "~/Research/dengue/tables/"
# directory to write out fcs files
output.fcs.path <- "~/Research/dengue/fcs/post_batch_correction/"
# directory for creation of batch normalization images
images.path <- "~/Research/dengue/images/batch_normalization/"



##### FUNCTIONS #####

combineFiles <- function(frames=fs) {
  # reads in a flowset and ouputs a list of expression matrices
  # Inputs:
  #   frames - a flowset
  # Outputs:
  #   matrices - a list of expression matrices
  print("Combining files")
  matrices <- list()
  for (i in seq(frames)) {
    matrices[[i]] <- frames[[i]]@exprs
    colnames(matrices[[i]]) <- pData(frames[[i]]@parameters)$desc
  }
  names(matrices) <- frames@phenoData@data$name
  return(matrices)
}



asinTransform <- function(matrices, channels=ch) {
  # asinh transforms
  # Inputs:
  #   matrices - list of expression matrices
  #   channels - vector of column number
  # Outputs:
  #   matrices
  print("Asinh transforming")
  for (i in seq(matrices)) {
    matrices[[i]][, channels] <- asinh(matrices[[i]][, channels]/5)
  }
  return(matrices)
}


plotViolinsControl <- function(matrices, channels=ch, path=images.path, suffix, percentile=0.99) {
  # Plots violins of markers to be transformed
  # Inputs:
  #   matrices - list of expression matrices
  #   channels - vector of column number
  #   path - string of directory for images
  #   suffix - subdirectory name and appended to each file
  #   percentile - percentile to visualize
  # Outputs:
  #   png of violins
  print("Plotting violins")
  if (!dir.exists(path)) dir.create(path)
  path <- paste0(path, "/", suffix, "/")
  if (!dir.exists(path)) dir.create(path)
  for (i in channels) {
    temp <- data.table(matrix(nrow=0, ncol=2)) %>%
      setnames(c("expression", "file"))
    for (j in seq(matrices)) {
      temp <- data.table(expression=matrices[[j]][, i], file=j) %>%
        rbind(temp)
    }
    temp[, file:=as.factor(file)]
    ggplot(temp, aes(x=file, y=expression, fill=file)) + geom_violin(scale="width") +
      theme_bw() + theme(text=element_text(size=20), legend.position="none") +
      stat_summary(fun=median, geom="point", shape=18, size=2, position=position_dodge()) +
      stat_summary(fun=getPeak, geom="point", color="blue", shape=18, size=2, position=position_dodge()) +
      stat_summary(fun=quantile, fun.args=list(probs=c(percentile)), geom="point",
                   color="red", shape=18, size=2, position=position_dodge()) +
      labs(title=colnames(matrices[[1]])[i])
    ggsave(paste0(path, colnames(matrices[[1]])[i], "_", suffix, ".png"), height=4, width=5)
  }
  return(matrices)
}


plotViolinsTotal <- function(matrices, channels=ch, path=images.path, suffix, percentile=0.99) {
  # Plots violins of markers to be transformed
  # Inputs:
  #   matrices - list of expression matrices
  #   channels - vector of column number
  #   path - string of directory for images
  #   suffix - subdirectory name and appended to each file
  #   percentile - percentile to visualize
  # Outputs:
  #   png of violins
  print("Plotting violins")
  if (!dir.exists(path)) dir.create(path)
  path <- paste0(path, "/", suffix, "/")
  if (!dir.exists(path)) dir.create(path)
  batch <- names(matrices) %>%
    substring(1,1) %>%
    as.integer()
  for (i in channels) {
    temp <- data.table(matrix(nrow=0, ncol=3)) %>%
      setnames(c("expression", "file", "batch"))
    n <- 0
    for (j in seq(matrices)) {
      n <- n + 1
      temp <- data.table(expression=matrices[[j]][, i], file=n, batch=batch[j]) %>%
        rbind(temp)
      if(j==length(matrices)) next
      if (batch[j] != batch[j+1]) n <- 0
    }
    temp[, `:=`(file=as.factor(file), batch=as.factor(batch))]
    ggplot(temp, aes(x=batch, y=expression, fill=batch)) + geom_violin(scale="width") +
      theme_bw() + theme(text=element_text(size=20), legend.position="none") +
      stat_summary(fun=median, geom="point", shape=18, size=2, position=position_dodge()) +
      stat_summary(fun=getPeak, geom="point", color="blue", shape=18, size=2, position=position_dodge()) +
      stat_summary(fun=quantile, fun.args=list(probs=c(percentile)), geom="point",
                   color="red", shape=18, size=2, position=position_dodge()) +
      facet_wrap(~file, nrow=5) +
      labs(title=colnames(matrices[[1]])[i])
    ggsave(paste0(path, colnames(matrices[[1]])[i], "_", suffix, ".png"), height=16, width=10)
  }
  return(matrices)
}


getPeak <- function(x) {
  # Returns the mode of a continuous distribution
  # Inputs:
  #   x - a numeric vector
  # Outputs:
  #   the numeric value of the mode
  d <- density(x)
  return(d$x[which.max(d$y)])
}


writeMedian <- function(matrices, channels=ch, output.path=tables.path, filename) {
  # Writes a csv of the median of each marker for each matrix
  # Inputs:
  #   matrices - list of expression matrices
  #   channels - vector of column number
  #   output.path - directory to write csv
  #   filename - csv filename
  # Outputs:
  #   matrices
  #   csv file
  print("Calculating median")
  output <- sapply(matrices, function(x) apply(x[, channels], 2, median)) %>%
    t() %>%
    as.data.table()
  fwrite(output, paste0(output.path, filename))
  return(matrices)
}


subtractMedian <- function(matrices, channels=ch, input.path=tables.path, filename, sdev=asinh(0.1)) {
  # Aligns medians of each marker between batches through subtraction
  # Inputs:
  #   matrices - list of expression matrices
  #   channels - vector of column number
  #   input.path - directory to read csv
  #   filename - csv filename
  #   sdev - standard deviation of random normal distribution to replace marker negative cells
  # Outputs:
  #   matrices
  print("Subtracting median")
  medians <- fread(paste0(input.path, filename)) %>%
    as.matrix() %>%
    apply(2, function(x) x-min(x))
  batch <- names(matrices) %>%
    substring(1,1) %>%
    as.integer()
  for (i in seq(matrices)) {
    matrices[[i]][, channels] <- sweep(matrices[[i]][, channels], 2, medians[batch[i], ])
  
    # replace negative numbers with random distribution around zero
    count <- length(matrices[[i]][matrices[[i]]<=0])
    set.seed(666)
    matrices[[i]][matrices[[i]]<=0] <- abs(rnorm(count, mean=0, sd=sdev))
  }
  return(matrices)
}


writePercentileNormalize <- function(matrices, channels=ch, percentile=0.99, output.path=tables.path, filename) {
  # percentile normalizes channels
  # Inputs:
  #   matrices - list of expression matrices
  #   channels - vector of column number
  #   percentile - percentile to normalize to
  #   output.path - directory to write csv
  #   filename - csv filename
  # Outputs:
  #   matrices
  output <- mapply(function(x) apply(x, 2, quantile, probs=c(percentile)), matrices) %>%
    t() %>%
    as.data.table() %>%
    .[, channels, with=F]
  fwrite(output, paste0(output.path, filename))
  return(matrices)
}


percentileNormalize <- function(matrices, channels=ch, input.path=tables.path, filename) {
  # percentile normalizes channels
  # Inputs:
  #   matrices - list of expression matrices
  #   channels - vector of column number
  #   input.path - directory to read csv
  #   filename - csv filename
  # Outputs:
  #   matrices
  print("Percentile normalizing")
  multipliers <- fread(paste0(input.path, filename)) %>%
    as.matrix() %>%
    apply(2, function(x) x/min(x))
  batch <- names(matrices) %>%
    substring(1,1) %>%
    as.integer()
  for (i in seq(matrices)) {
    matrices[[i]][, channels] <- sweep(matrices[[i]][, channels], 2, multipliers[batch[i], ], FUN = "/")
  }
  return(matrices)
}


reverseAsinTransform <- function(matrices, channels=ch) {
  # reverse asinh transforms
  # Inputs:
  #   matrices - list of expression matrices
  #   channels - vector of column number
  # Outputs:
  #   matrices
  print("Reverse asinh transforming")
  for (i in seq(matrices)) {
    matrices[[i]][, channels] <- sinh(matrices[[i]][, channels])*5
  }
  return(matrices)
}


updateFS <- function(matrices, frames=fs) {
  # updates the exprs matrices of the flowset with the batch-normalized values
  # Inputs:
  #   matrices - list of expression matrices
  #   frames - a flowset
  # Outputs:
  #   frames
  print("Writing back to flow set")
  for (i in seq(frames)) {
    colnames(matrices[[i]]) <- pData(frames[[i]]@parameters)$name
    frames[[i]]@exprs <- matrices[[i]]
  }
  return(frames)
}



##### MAIN #####

### Adjust panels to remove beadDist since run 6 lacks this column
paneleditor_GUI()

### Run control samples
fs <- read.flowSet(path=input.fcs.path, pattern="C-PA1", transformation=F)
# View channel information:
pData(fs[[1]]@parameters)
ch <- c(3, 15:52, 66) # protein channels

control.matrices <- combineFiles() %>%
  asinTransform() %>%
  plotViolinsControl(suffix="pre_control") %>%
  writeMedian(filename="control_medians.csv") %>%
  subtractMedian(filename="control_medians.csv") %>%
  plotViolinsControl(suffix="subtracted_control") %>%
  writePercentileNormalize(filename="control_percentiles.csv") %>%
  percentileNormalize(filename="control_percentiles.csv") %>%
  plotViolinsControl(suffix="post_control")


### Run samples
fs <- read.flowSet(path=input.fcs.path, pattern=".fcs", transformation=F)
pData(fs[[1]]@parameters)
ch <- c(3, 15:52, 66) # protein channels

updated.fs <- combineFiles() %>%
  asinTransform() %>%
  plotViolinsTotal(suffix="pre_total") %>%
  subtractMedian(filename="control_medians.csv") %>%
  plotViolinsTotal(suffix="subtracted_total") %>%
  percentileNormalize(filename="control_percentiles.csv") %>%
  plotViolinsTotal(suffix="post_total") %>%
  reverseAsinTransform() %>%
  updateFS()

if (!dir.exists(output.fcs.path)) dir.create(output.fcs.path)
write.flowSet(updated.fs, output.fcs.path)