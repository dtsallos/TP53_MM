#### AllDSRT_data_functional_coding_style ######
## I edited the previous code using a functional coding style learnt during a course taken after the last update of the manuscript ##


#load libraries
library(tidyr)
library(dplyr)

# folders
main = "/MyelomaData/"
code = "/MyelomaData/TP53project/"
fdata = "/MyelomaData/data/"
here = paste0(code, "1_DSRT/significance/")


#### import tables ####
df <- read.csv( paste0(here, "DSRT.csv"  )) 

##### Analyse results #####
results <- data.frame(matrix(ncol = 18, nrow = 0))

colnames(results) <- c("drug", "functional_class", "mechanism_target",
                       "avg_DSS", "avg_DSS_TP53mut", "avg_DSS_monoTP53", "avg_DSS_other", 
                       "median_DSS", "median_DSS_TP53mut", "median_DSS_monoTP53", "median_DSS_other", 
                       "sd_DSS", "sd_DSS_TP53mut", "sd_DSS_monoTP53", "sd_DSS_other", 
                       "num_samples", "num_TP53mut", "num_monoTP53")

## Split data by drug
df_by_drug <- split(df, df$drug)

## Function that computes the summary for one drug subset
summarise_drug <- function(temp) {
  # Extract functional class and mechanism target
  functional_class   <- unique(temp$functional.class)[1]
  mechanism_target   <- unique(temp$mechanism.target)[1]
  
  # Count samples
  total_samples      <- nrow(temp)
  TP53mut_samples    <- sum(temp$TP53_inclusive_status == "TP53_mutation")
  monoTP53_samples   <- sum(temp$TP53_inclusive_status == "TP53_deletion")
  
  # Predefine logicals once 
  is_mut   <- temp$TP53_inclusive_status == "TP53_mutation"
  is_del   <- temp$TP53_inclusive_status == "TP53_deletion"
  is_other <- !(temp$TP53_inclusive_status %in% c("TP53_mutation", "TP53_deletion"))
  
  # Means
  avg_DSS           <- mean(temp$DSS, na.rm = TRUE)
  avg_DSS_TP53mut   <- mean(temp$DSS[is_mut],   na.rm = TRUE)
  avg_DSS_monoTP53  <- mean(temp$DSS[is_del],   na.rm = TRUE)
  avg_DSS_other     <- mean(temp$DSS[is_other], na.rm = TRUE)
  
  # Medians
  median_DSS           <- median(temp$DSS, na.rm = TRUE)
  median_DSS_TP53mut   <- median(temp$DSS[is_mut],   na.rm = TRUE)
  median_DSS_monoTP53  <- median(temp$DSS[is_del],   na.rm = TRUE)
  median_DSS_other     <- median(temp$DSS[is_other], na.rm = TRUE)
  
  # SDs
  sd_DSS           <- sd(temp$DSS, na.rm = TRUE)
  sd_DSS_TP53mut   <- sd(temp$DSS[is_mut],   na.rm = TRUE)
  sd_DSS_monoTP53  <- sd(temp$DSS[is_del],   na.rm = TRUE)
  sd_DSS_other     <- sd(temp$DSS[is_other], na.rm = TRUE)
  
  # Return a data.frame row with correct types
  data.frame(
    drug               = temp$drug[1],
    functional_class   = functional_class,
    mechanism_target   = mechanism_target,
    
    avg_DSS            = avg_DSS,
    avg_DSS_TP53mut    = avg_DSS_TP53mut,
    avg_DSS_monoTP53   = avg_DSS_monoTP53,
    avg_DSS_other      = avg_DSS_other,
    
    median_DSS         = median_DSS,
    median_DSS_TP53mut = median_DSS_TP53mut,
    median_DSS_monoTP53= median_DSS_monoTP53,
    median_DSS_other   = median_DSS_other,
    
    sd_DSS             = sd_DSS,
    sd_DSS_TP53mut     = sd_DSS_TP53mut,
    sd_DSS_monoTP53    = sd_DSS_monoTP53,
    sd_DSS_other       = sd_DSS_other,
    
    total_samples      = total_samples,
    TP53mut_samples    = TP53mut_samples,
    monoTP53_samples   = monoTP53_samples,
    
    stringsAsFactors = FALSE
  )
}

## Apply over groups and bind rows
results <- do.call(rbind, lapply(split(df, df$drug), summarise_drug))
row.names(results) <- NULL

write.csv(results, paste0(here, "/drug_DSS_statistics.csv"), row.names = FALSE)
