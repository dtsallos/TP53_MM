# AllDSRT_data_loop_TP53

#load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggprism)
library(broom)
library(ggsignif)
library(DT)

# folders
main = "/MyelomaData/"
code = "/MyelomaData/TP53project/"
fdata = "/MyelomaData/data/"
here = paste0(code, "1_DSRT/significance/")


#### import tables ####
df2 <- read.csv( paste0(here, "DSRT.csv"  )) 

##### Analyse results #####
results <- data.frame(matrix(ncol = 18, nrow = 0))

colnames(results) <- c("drug", "functional_class", "mechanism_target",
                       "avg_DSS", "avg_DSS_TP53mut", "avg_DSS_monoTP53", "avg_DSS_other", 
                       "median_DSS", "median_DSS_TP53mut", "median_DSS_monoTP53", "median_DSS_other", 
                       "sd_DSS", "sd_DSS_TP53mut", "sd_DSS_monoTP53", "sd_DSS_other", 
                       "num_samples", "num_TP53mut", "num_monoTP53")

drugs <- unique(df2$drug)
for(i in 1:length(drugs)) {  
  temp <- df2[which(df2$drug == drugs[i]),]
  
  # Extract functional class and mechanism target
  functional_class <- unique(temp$functional.class)[1]
  mechanism_target <- unique(temp$mechanism.target)[1]
  
  # Count total samples, TP53_mutation samples, and monoTP53 samples for each drug
  total_samples <- nrow(temp)
  TP53mut_samples <- sum(temp$TP53_inclusive_status == "TP53_mutation")
  monoTP53_samples <- sum(temp$TP53_inclusive_status == "TP53_deletion")
  
  # Compute statistics
  avg_DSS <- mean(temp$DSS, na.rm = TRUE)
  avg_DSS_TP53mut <- mean(temp$DSS[temp$TP53_inclusive_status == "TP53_mutation"], na.rm = TRUE)
  avg_DSS_monoTP53 <- mean(temp$DSS[temp$TP53_inclusive_status == "TP53_deletion"], na.rm = TRUE)
  avg_DSS_other <- mean(temp$DSS[!(temp$TP53_inclusive_status %in% c("TP53_mutation", "TP53_deletion"))], na.rm = TRUE)
  
  median_DSS <- median(temp$DSS, na.rm = TRUE)
  median_DSS_TP53mut <- median(temp$DSS[temp$TP53_inclusive_status == "TP53_mutation"], na.rm = TRUE)
  median_DSS_monoTP53 <- median(temp$DSS[temp$TP53_inclusive_status == "TP53_deletion"], na.rm = TRUE)
  median_DSS_other <- median(temp$DSS[!(temp$TP53_inclusive_status %in% c("TP53_mutation", "TP53_deletion"))], na.rm = TRUE)
  
  sd_DSS <- sd(temp$DSS, na.rm = TRUE)
  sd_DSS_TP53mut <- sd(temp$DSS[temp$TP53_inclusive_status == "TP53_mutation"], na.rm = TRUE)
  sd_DSS_monoTP53 <- sd(temp$DSS[temp$TP53_inclusive_status == "TP53_deletion"], na.rm = TRUE)
  sd_DSS_other <- sd(temp$DSS[!(temp$TP53_inclusive_status %in% c("TP53_mutation", "TP53_deletion"))], na.rm = TRUE)
  
  # Store results
  results[i,] <- c(drugs[i], functional_class, mechanism_target, avg_DSS, avg_DSS_TP53mut, avg_DSS_monoTP53, avg_DSS_other, 
                   median_DSS, median_DSS_TP53mut, median_DSS_monoTP53, median_DSS_other, 
                   sd_DSS, sd_DSS_TP53mut, sd_DSS_monoTP53, sd_DSS_other, 
                   total_samples, TP53mut_samples, monoTP53_samples)
}

# Convert numeric columns back from character type
results[, 4:18] <- lapply(results[, 4:18], as.numeric)

write.csv(results, paste0(here, "/drug_DSS_statistics.csv"), row.names = FALSE)
