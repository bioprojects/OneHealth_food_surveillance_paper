d_human <- read.table("human_human_mge_combined.sorted.fna.len95_iden95.presence_absence.rowInd.txt", header=T, sep="\t")
d_animal <- read.table("animal_human_mge_combined.sorted.fna.len95_iden95.presence_absence.rowInd.txt", header=T, sep="\t")
d_food   <- read.table("food_human_mge_combined.sorted.fna.len95_iden95.presence_absence.rowInd.txt", header=T, sep="\t")

c_human_freq <- apply(d_human[,-1], 2, sum, na.rm=T)/261
c_animal_freq <- apply(d_animal[,-1], 2, sum, na.rm=T)/364
c_food_freq   <- apply(d_food[,-1], 2, sum, na.rm=T)/57

names(c_human_freq) == names(c_food_freq)
names(c_animal_freq) == names(c_food_freq)

df <- data.frame(human=c_human_freq, food=c_food_freq, animal=c_animal_freq)

library(pheatmap)
pdf("mge_freq_heatmap.pdf")
pheatmap(
  t(as.matrix(df)),
  fontsize_col = 3, 
  fontsize_row = 8, 
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  treeheight_col = 50,
  clustering_distance_rows = "euclidean",
  breaks = seq(0, 1, length.out = 101)
)
dev.off()

write.csv(df, file="mge_freq.csv")
