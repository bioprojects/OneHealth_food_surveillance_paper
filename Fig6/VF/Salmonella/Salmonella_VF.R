d_meta <- read.csv("Se1613_phandango.csv", header=T) 
d_vir_pre_abs0  <- read.table("Se_VF_join.txt", header=T, sep="\t")

c_col_order <- c(
"lpfA__STM3640"
,"lpfB__STM3639"
,"lpfC__STM3638"
,"lpfD__STM3637"
,"lpfE__STM3636"
,"pefA__PSLT018"
,"pefB__PSLT019"
,"pefC__PSLT017"
,"pefD__PSLT016"
,"yehD__SPA_RS03495"
,"SPA_RS03500"
,"SPA_RS03505"
,"SPA_RS03510"
,"ratB__STM2514"
,"safA__STM0299"
,"safB__STM0300"
,"safC__STM0301"
,"safD__STM0302"
,"sefA__STY_RS23155"
,"sefB__STY_RS23160"
,"sefC__STY_RS23165"
,"sefD__SEN_RS22085"
,"shdA__STM2513"
,"STY_RS00955"
,"STY_RS00950"
,"staC__STY_RS00945"
,"STY_RS00940"
,"staE__STY_RS00935"
,"STY_RS00930"
,"STY_RS00925"
,"stcA__STM2152"
,"stcB__STM2151"
,"stcC__STM2150"
,"stcD__STM2149"
,"SPA_RS14285"
,"STY_RS14630"
,"STY_RS14635"
,"STY_RS14640"
,"STY_RS14645"
,"STY_RS14650"
,"stfA__STM0195"
,"stfC__STM0196"
,"stfD__STM0197"
,"stfE__STM0198"
,"stfF__STM0199"
,"stfG__STM0200"
,"sthD__STM4592"
,"stiA__STM0177"
,"stiB__STM0176"
,"stiC__STM0175"
,"stiH__STM0174"
,"SEHA_RS24060"
,"STM4571"
,"SNSL254_RS24290"
,"stjB__STM4572"
,"stjC__STM4573"
,"SPA_RS00935"
,"SPA_RS00930"
,"SPA_RS00925"
,"SPA_RS00920"
,"SPA_RS00915"
,"SPA_RS00910"
,"SPA_RS00905"
,"tcfD__STY_RS01595"
,"avrA__STM2865"
,"sopE__STY_RS22015"
,"gogB__STM2584"
,"spvC__PSLT038"
,"spvD__PSLT037"
,"sseI__STM1051"
,"STM4157"
,"STM2137"
,"sspH2__STM2241"
,"spvB__PSLT039"
,"cdtB__STY_RS08895"
,"pltA__STY_RS08910"
,"pltB__STY_RS08915"
,"PSLT046"
,"rcK__PSLT009"
,"sodC__STM1044"
)

c_target_lineages <- c(
"Schwarzengrund", 
"Infantis", 
"Manhattan", 
"Heidelberg", 
"Enteritidis", 
"Blockley", 
"Minnesota", 
"O4:i:-", 
"Thompson", 
"Typhimurium", 
"Saintpaul", 
"Newport", 
"Stanley", 
"Hadar"
)



## Campy 
#d_meta <- read.csv("Campy682_phandango.csv", header=T)       # ­‚È‚­‚Æ‚àstrainID‚Æhost‚Æmajor_serovar‚Ì3—ñ
#d_vir_pre_abs0  <- read.table("Campy_VF_join.txt", header=T, sep="\t")
#
#c_col_order <- c(
#"flaB__Cj1338c"
#,"flgE__Cj1729c"
#,"Cj0887c"
#,"fliD__Cj0548"
#,"pseD__Cj1333"
#)
#
#c_target_lineages <- c(
#"ST-21 complex",
#"ST-22 complex",
#"ST-48 complex",
#"ST922",
#"ST-443 complex",
#"ST-42 complex",
#"ST-45 complex",
#"ST-464 complex",
#"ST-61 complex",
#"ST-52 complex",
#"ST-828 complex"
#)



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