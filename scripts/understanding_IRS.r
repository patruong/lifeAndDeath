
# load libraries (this gets ggplot2 and related libraries)
library(tidyverse)

# these are from Bioconductor
library(limma) 
library(edgeR) 
library(sva)

# we use this for the pairwise correlation plots
library(psych)

# read the Supplemental 01 file (saved as a CSV export from XLSX file)
data_start <- read_csv("iovs-58-13-55_s01.csv")

# filter out proteins not seen in all three runs
data_no_na <- na.omit(data_start)

# fix the column headers
col_headers <- colnames(data_no_na) 
col_headers <- str_replace(col_headers, " {2,3}", " ") 
col_headers <- str_replace(col_headers, "Reporter ion intensities ", "") 
colnames(data_no_na) <- col_headers

# save the annotations (gene symbol and protein accession) and remove from data frame
annotate_df <- data_no_na[1:2] 
data_raw <- as.data.frame(data_no_na[3:20]) 
row.names(data_raw) <- annotate_df$`Protein Accession No.`

head(data_raw)

# let's see what the starting data look like
boxplot(log2(data_raw), col = rep(c('red', 'green', 'blue'), each = 6), 
        notch = TRUE, main = 'RAW data: Exp1 (red), Exp2 (green), Exp3 (blue)',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity')

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_raw), col = rep(c('red', 'green', 'blue'), 6), 
              main = 'Raw data')


# check the column totals (per channel sums)
format(round(colSums(data_raw), digits = 0), big.mark = ",")

# separate the TMT data by experiment
# we do not need to do this for the normalization factor calculation here,
# but we will need these data frames for the IRS step below.
exp1_raw <- data_raw[c(1:6)]
exp2_raw <- data_raw[c(7:12)]
exp3_raw <- data_raw[c(13:18)]

# first basic normalization is to adjust each TMT experiment to equal signal per channel
# figure out the global scaling value
target <- mean(c(colSums(exp1_raw), colSums(exp2_raw), colSums(exp3_raw)))

# do the sample loading normalization before the IRS normalization
# there is a different correction factor for each column
norm_facs <- target / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3_raw)
exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")

# make a pre-IRS data frame after sample loading normalizations
data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl)

# see what the SL normalized data look like
boxplot(log2(data_sl), col = rep(c("red", "green", "blue"), each = 6), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nExp1 (red), Exp2 (green), Exp3 (blue)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity')

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_sl), col = rep(c("red", "green", "blue"), 6), main = "SL normalization")


raw_tmm <- calcNormFactors(data_raw)
data_raw_tmm <- sweep(data_raw, 2, raw_tmm, FUN = "/") # raw data after TMM on original scale

# perform TMM on the SL-normed data and visualize resulting distributions
sl_tmm <- calcNormFactors(data_sl)
data_sl_tmm <- sweep(data_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

boxplot(log2(data_sl_tmm), notch = TRUE, col = rep(c("red", "green", "blue"), each = 6), 
        main = "TMM normalization of SL data\nExp1 (red), Exp2 (green), Exp3 (blue)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity')

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_sl_tmm), col = rep(c("red", "green", "blue"), 6), main = "SL/TMM normalization")

# make new data frame with row sums from each frame
irs <- tibble(rowSums(exp1_sl), rowSums(exp2_sl), rowSums(exp3_sl))
colnames(irs) <- c("sum1", "sum2", "sum3")

# get the geometric average intensity for each protein
irs$average <- apply(irs, 1, function(x) exp(mean(log(x))))

# compute the scaling factor vectors
irs$fac1 <- irs$average / irs$sum1
irs$fac2 <- irs$average / irs$sum2
irs$fac3 <- irs$average / irs$sum3

# make new data frame with IRS normalized data
data_irs <- exp1_sl * irs$fac1
data_irs <- cbind(data_irs, exp2_sl * irs$fac2)
data_irs <- cbind(data_irs, exp3_sl * irs$fac3)

# see what the IRS data look like
boxplot(log2(data_irs), col = rep(c("red", "green", "blue"), each = 6), 
        main = "Internal Reference Scaling (IRS) normalized data: \nExp1 (red), Exp2 (green), Exp3 (blue)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity', notch = TRUE)

# can also look at density plots (like a distribution histogram)    
plotDensities(log2(data_irs), col = rep(c("red", "green", "blue"), 6), main = "IRS data")

