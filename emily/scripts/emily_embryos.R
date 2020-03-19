# Emily filtering

# load packages
require(dplyr)

# set wd
setwd(here::here("emily"))

# import and round data like you normally would
countsTable <- read.table("Dm_countsMatrix.txt", header = TRUE)
countsTableRound <- round(countsTable)

# get a list of the genes that are in all samples (only about 6500 transcripts)
# or whatever genes meet the cutoff you decide
keep_genes <- countsTableRound %>%
  rownames_to_column("gene") %>%
  pivot_longer(contains("_"), names_to = "sample", values_to = "counts") %>%
  group_by(gene) %>%
  filter(counts != 0) %>%
  count() %>%
  # this number below: n == 100 you can adjust if you want 
  # say you want transripts that appear in 90% of the samples 
  # you can say n >= 99 or something like that (brings it up to 9000 transcripts)
  filter(n > 27) 

# I think this prefiltering is an interesting idea, 
# I don't think it is too typical to require that a transcript was found in all samples
# DESeq2 normalization and independent filtering try to address these sorts of problems
# https://www.youtube.com/watch?v=UFB993xufUU
# https://www.youtube.com/watch?v=Gi0JdrxRq5s
# Maybe do it so it so the transcript is found in at least half or 25% of the samples 
# 25% or > 27 out of the 110 takes it to 18,002 transrctipts which seems reasonable
# I just worry that if your hypothesis has to do with the MZT happening earlier 
# in the tropical flies then those zygotic transcripts wouldn't be present at 
# all in a lot of the temperate embryos so you would lose a lot of potenitally cool differential expression

# then you need to filter the genes in the counts matrix 
# so that they are only genes that match the list created above
countsTableRound <- countsTableRound %>%
  rownames_to_column("gene") %>%
  filter(gene %in% keep_genes$gene) %>%
  column_to_rownames("gene") 



