# brains acclimation temperature differentially expressed genes ----------------

directory_results <- here::here("cahan/results")


# after running the SARtools script, save the deg results to .csvs
write.csv(as.data.frame(out.DESeq2$results$hot_vs_ctrl), 
          "brain_deg_hot_results.csv")
write.csv(as.data.frame(out.DESeq2$results$cold_vs_ctrl), 
          "brain_deg_cold_results.csv")

# load the differentially expressed genes results
setwd(directory_results)
res_cold <- read.csv(list.files(pattern = "cold_results"))
res_hot <- read.csv(list.files(pattern = "hot_results"))


