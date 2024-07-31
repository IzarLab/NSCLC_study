library(ggplot2)

####dotplot
mylist <- read.csv(file="interactions.csv")
mylist <- as.data.frame(mylist)
all_rows = c()

clusters <- c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6", "MP7", "MP8")

for (i in  1:ncol(mylist)) {
  
  
  val = paste0("MP",i)
  genes <- mylist[[val]]
  
  # Load Enrichr and check website status
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Set to human genes
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE
  if (websiteLive) head(dbs)
  
  # read signature
  dbs <- c("MSigDB_Hallmark_2020", "GO_Biological_Process_2021")
  if (websiteLive) {
    enriched <- enrichr(genes, dbs)
    mut_enr_ch <- mutate(enriched[[1]], qscore = -log10(Adjusted.P.value))
    mut_enr_go <- mutate(enriched[[2]], qscore = -log10(Adjusted.P.value))
    
    mut_enr <- rbind.data.frame(mut_enr_ch, mut_enr_go)
    filtered_mut_enr <- subset(mut_enr[1:20, ], P.value <0.05)
    
    all_rows <- c(all_rows,filtered_mut_enr$Term )
    
  }
}
all_rows <- unique(all_rows)

all_dfs = list()
for (i in  1:ncol(mylist)) {
  
 # i=1
  val = paste0("MP",i)
  genes <- mylist[[val]]
  
  # Load Enrichr and check website status
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Set to human genes
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE
  if (websiteLive) head(dbs)
  
  # read signature
  dbs <- c("MSigDB_Hallmark_2020", "GO_Biological_Process_2021")
  if (websiteLive) {
    enriched <- enrichr(genes, dbs)
    mut_enr_ch <- mutate(enriched[[1]], qscore = -log10(Adjusted.P.value))
    mut_enr_go <- mutate(enriched[[2]], qscore = -log10(Adjusted.P.value))
    
    mut_enr <- rbind.data.frame(mut_enr_ch, mut_enr_go)
    filtered_mut_enr <- subset(mut_enr, Term %in% all_rows)
    filtered_mut_enr$interaction <- val
    
    all_dfs <- append(all_dfs, list(filtered_mut_enr))
    print(filtered_mut_enr$Adjusted.P.value)
    
    
  }
}

combined_df_all <- bind_rows(all_dfs)

pdf("bubbleplot_NMF.pdf", width = 7, height = 7)
combined_df_all$interaction <- as.factor(combined_df_all$interaction)

desired_order <- c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6", "MP7", "MP8")

# Convert the Category column to a factor with the specified order
#combined_df_all$interaction <- factor(combined_df_all$interaction, levels = desired_order)

# Adjusted ggplot code for the combined_df DataFrame
ggp <- ggplot(combined_df_all, aes(x = interaction, y = Term, size = qscore, 
                                   colour = ifelse(P.value < 0.05, "red", "grey"))) +
  geom_point() + 
  scale_colour_identity() +  # Use exact colors specified in the data
  theme_minimal() +
  labs(title = "Bubble Plot of Interaction Terms",
       x = "Interaction",
       y = "Term",
       size = "Q Score",
       color = "P-Value Significance") +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5))

# Remove legend for color if not needed
print(ggp)
dev.off()
