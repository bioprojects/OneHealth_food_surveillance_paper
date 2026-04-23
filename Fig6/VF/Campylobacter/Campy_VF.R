
## Campy 

d_meta <- read.csv("Campy682_phandango.csv", header=T)       # ­‚È‚­‚Æ‚àstrainID‚Æhost‚Æmajor_serovar‚Ì3—ñ
d_vir_pre_abs0  <- read.table("Campy_VF_join.txt", header=T, sep="\t")

c_col_order <- c(
"flaB__Cj1338c"
,"flgE__Cj1729c"
,"Cj0887c"
,"fliD__Cj0548"
,"pseD__Cj1333"
)

c_target_lineages <- c(
"ST-21 complex",
"ST-22 complex",
"ST-48 complex",
"ST922",
"ST-443 complex",
"ST-42 complex",
"ST-45 complex",
"ST-464 complex",
"ST-61 complex",
"ST-52 complex",
"ST-828 complex"
)



dd <- merge(d_meta, d_vir_pre_abs0, by.x="strainID", by.y="indID")
dim(dd)

d_vir_pre_abs <- merge(d_meta[,c("strainID","major_serovar")], d_vir_pre_abs0, by.x="strainID", by.y="indID")

d_out <- data.frame()
for (each_vir in names(d_vir_pre_abs)[-c(1,2)]) {
  res_table <- table(dd$host, dd[, each_vir])

  print(each_vir)
  print(res_table)

  if (ncol(res_table) == 1) {

	  d_line <- data.frame(
	            vir_gene=each_vir
	            , type="all strains carry it"
	            , p_animal_vs_food=NA
	            , OR_animal_vs_food=NA
	            , p_animal_vs_human=NA
	            , OR_animal_vs_human=NA
	            , p_food_vs_human=NA
	            , OR_food_vs_human=NA
	            , freq_in_animal=NA
	            , freq_in_food=NA
	            , freq_in_human=NA
	            )

	} else {
		
		stopifnot (nrow(res_table) == 3) 
		
    res_fisher_animal_vs_food  <- fisher.test(res_table[c(1,2), ])
    res_fisher_animal_vs_human <- fisher.test(res_table[c(1,3), ])
    res_fisher_food_vs_human   <- fisher.test(res_table[c(2,3), ])
    
    freq_in_animal <- res_table[1, 2]/sum(res_table[1,])
    freq_in_food   <- res_table[2, 2]/sum(res_table[2,])
    freq_in_human  <- res_table[3, 2]/sum(res_table[3,])

    d_line <- data.frame(
              vir_gene=each_vir
              , type="polymorphic across 3 groups"
              , p_animal_vs_food=res_fisher_animal_vs_food$p.value
              , OR_animal_vs_food=res_fisher_animal_vs_food$estimate
              , p_animal_vs_human=res_fisher_animal_vs_human$p.value
              , OR_animal_vs_human=res_fisher_animal_vs_human$estimate
              , p_food_vs_human=res_fisher_food_vs_human$p.value
              , OR_food_vs_human=res_fisher_food_vs_human$estimate
              , freq_in_animal=freq_in_animal
              , freq_in_food=freq_in_food
              , freq_in_human=freq_in_human
              )

	}

  d_out <- rbind(d_out, d_line)

}
dim(d_out)

d_out[d_out$vir_gene=="safC__STM0301", ]

order_by_p_food_vs_human <- order(d_out$p_food_vs_human)
d_out <- d_out[order_by_p_food_vs_human, ]

d_out$Padjust_food_vs_human   <- p.adjust(d_out$p_food_vs_human,   method="bonferroni")
d_out$Padjust_animal_vs_human <- p.adjust(d_out$p_animal_vs_human, method="bonferroni")
#d_out$Padjust_animal_vs_food  <- p.adjust(d_out$p_animal_vs_food,  method="bonferroni")

sum(d_out$Padjust_food_vs_human < 0.05, na.rm=T)
sum(d_out$Padjust_animal_vs_human < 0.05, na.rm=T)


d_out_signif <- rbind(
   d_out[!is.na(d_out$Padjust_food_vs_human)   & d_out$Padjust_food_vs_human   < 0.05, ]
  ,d_out[!is.na(d_out$Padjust_animal_vs_human) & d_out$Padjust_animal_vs_human < 0.05, ]
)

d_out_signif_unique <- unique(d_out_signif)

d_out_signif_numeric         <- data.frame( 
                                  t(d_out_signif_unique[, c("freq_in_human","freq_in_food", "freq_in_animal")])
                                )
names(d_out_signif_numeric)  <- d_out_signif_unique[, "vir_gene"]

d_out_signif_numeric_ordered <- d_out_signif_numeric[, c_col_order]

library(pheatmap)

pdf("vir_freq_matrix_3categories_heatmap.pdf")
pheatmap(
  as.matrix(d_out_signif_numeric_ordered),
  fontsize_col = 4, 
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  treeheight_col = 50,
  clustering_distance_rows = "euclidean"
)
dev.off()
write.csv(d_out_signif_numeric, file="vir_freq_matrix_3categories_heatmap.csv")



library(eulerr)

A  <- d_out[!is.na(d_out$Padjust_food_vs_human)   & d_out$Padjust_food_vs_human   < 0.05, 1]
B  <- d_out[!is.na(d_out$Padjust_animal_vs_human) & d_out$Padjust_animal_vs_human < 0.05, 1]
#C  <- d_out[!is.na(d_out$Padjust_animal_vs_food)  & d_out$Padjust_animal_vs_food < 0.05, 1]


c_ABuniq <- A 
#c_ABCuniq <- A
for (i in 1:length(B)) {
  if (B[i] %in% c_ABuniq) {
    next
  } else {
    c_ABuniq <- c(c_ABuniq, B[i])
  }
}

dim(
  d_vir_pre_abs[, c_ABuniq]
)

d_vir_pre_abs_signif <- d_vir_pre_abs[, c("strainID", "major_serovar", c_ABuniq)]

d_vir_freq_out <- data.frame()
for (i in 3:ncol( d_vir_pre_abs_signif ) ) {
  res_table <- table( d_vir_pre_abs_signif[,c(2,i)] )

  c_target_order <- c()
  for (each_target in c_target_lineages) {
    #print(each_target)
    c_target_order <- c(c_target_order, which(rownames(res_table) == each_target) )
  }
  
  res_table_sorted <- res_table[c_target_order, ]
  
  res_freq_vector <- res_table_sorted[,2] / (res_table_sorted[,1] + res_table_sorted[,2])
  
  df_freq_vector <- data.frame(
                      c(names(d_vir_pre_abs_signif)[i]
                      , res_freq_vector
                      )
                    )

  if (ncol(d_vir_freq_out) == 0) {
    d_vir_freq_out <- df_freq_vector
  } else {
    d_vir_freq_out <- cbind(d_vir_freq_out, df_freq_vector)
  }
  
}

names(d_vir_freq_out) <- as.character(d_vir_freq_out[1, ])
d_vir_freq_out        <- d_vir_freq_out[-1, ]



d_vir_freq_out_numeric           <- data.frame(lapply(d_vir_freq_out, as.numeric))
rownames(d_vir_freq_out_numeric) <- rownames(d_vir_freq_out)

d_vir_freq_out_numeric_ordered <- d_vir_freq_out_numeric[, c_col_order]

library(pheatmap)
pdf("vir_freq_matrix_heatmap.pdf")
pheatmap(
  as.matrix(d_vir_freq_out_numeric_ordered),
  fontsize_col = 4, 
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  treeheight_col = 50,
  clustering_distance_rows = "euclidean"
)
dev.off()
write.csv(d_vir_freq_out_numeric, file="vir_freq_matrix_heatmap.csv")