################################################################################
# Set Workspace and Load in Libraries                                          #
################################################################################
library(Cardinal)
library(devtools)
load_all('C:/Users/gordon/Data/code/cardinalscripts')
library(stringr)
library(fuzzyjoin)
library(dplyr)
library(ggplot2)
library(reticulate)
kneed <- import('kneed')
numpy <- import('numpy')

################################################################################
# Initialize MultiProcessing for Windows and Set Working Directory             #
################################################################################
register(SnowParam(workers=1, log=TRUE, logdir=getwd()))
setwd('F:/ims_workflow_testing')

################################################################################
# Load in Datasets, Spot Lists, and ROI Lists                                  #
################################################################################
#fleximaging <- readMSIData('imzml/mini_ims_fleximaging.imzML', resolution=2000, units='ppm')
fleximaging <- as(readMSIData('imzml/mini_ims_fleximaging.imzML', resolution=2000, units='ppm'), 'MSContinuousImagingExperiment')
#scils <- readMSIData('imzml/mini_ims_scils_raw_import.imzML', resolution=2000, units='ppm')
scils <- as(readMSIData('imzml/mini_ims_scils_raw_import.imzML', resolution=2000, units='ppm'), 'MSContinuousImagingExperiment')

################################################################################
# Everything Else                                                              #
################################################################################
centroided(fleximaging) <- FALSE
centroided(scils) <- FALSE

# Normalization, Signal Smoothing, and Baseline Subtraction
preprocessed <- normalize(fleximaging, method='tic')
preprocessed <- smoothSignal(preprocessed, method='sgolay')
preprocessed <- reduceBaseline(preprocessed, method='median')
preprocessed <- process(preprocessed)
# Peak Picking, Alignment, and Filtering
peaks <- peakPick(preprocessed, method='simple', SNR=4)
peaks <- peakAlign(peaks, tolerance=NA, units='mz')
peaks <- peakFilter(peaks, freq.min=0.05)
processed_fleximaging <- process(peaks)
rm(preprocessed, peaks)

# Normalization, Signal Smoothing, and Baseline Subtraction
preprocessed <- normalize(scils, method='tic')
preprocessed <- smoothSignal(preprocessed, method='sgolay')
preprocessed <- reduceBaseline(preprocessed, method='median')
preprocessed <- process(preprocessed)
# Peak Picking, Alignment, and Filtering
peaks <- peakPick(preprocessed, method='simple', SNR=4)
peaks <- peakAlign(peaks, tolerance=NA, units='mz')
peaks <- peakFilter(peaks, freq.min=0.05)
processed_scils <- process(peaks)
rm(preprocessed, peaks)

set.seed(1120)
ssc_fleximaging <- spatialShrunkenCentroids(processed_fleximaging, r=1, k=4, s=3)
ssc_scils <- spatialShrunkenCentroids(processed_scils, r=1, k=4, s=3)
