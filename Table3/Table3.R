# This script was executed in the parent directory of the directory containing the TSV files of SNP results ("Table3_data").


curr_path = getwd()
file_names = list.files("./Table3_data")

CC_or_ST = c("ST-21_complex","ST-22_complex","ST-48_complex","ST-42_complex","ST-45_complex","ST-61_complex","ST922","ST-464_complex","ST-443_complex")
source_pair = c("human_vs_human", "human_vs_food", "human_vs_animal")

# An object that stores results, where each CC or ST is added as a row in the loop below.
temp = paste0("CC or ST\t", "human- vs human-derived isolates\t", "human- vs food-derived isolates\t", "human- vs animal-derived isolates\n")

# Writing the median and 25th/75th percentiles for each sectoral combination within each serovar
for (i in CC_or_ST){
	temp = paste0(temp, i, "\t")
file_name_CCorST = file_names[grep(i, file_names)]
	for (j in source_pair){
		target_file = file_name_CCorST[grep(j,file_name_CCorST)]
		if (length(target_file) == 0){ 
			temp = paste0(temp,"-")
			}
		else{
			df = read.table(paste0(curr_path, "/Table3_data/", target_file), header=T, sep="\t",comment.char="",quote ='"')
			if (nrow(df) == 0){
				temp = paste0(temp,"-")
				}
			else{
				SNPs = as.integer(df["dist"][,1])
				mdn = quantile(SNPs, type=1)[3]
				snp_25perc = quantile(SNPs, type=1)[2]
				snp_75perc = quantile(SNPs, type=1)[4]
				input_val = paste0(mdn, " (", snp_25perc, "-", snp_75perc, ")")
				temp = paste0(temp, input_val)
				}
			}
		if (j == "human_vs_animal"){
			temp = paste0(temp, "\n")
			}
		else{
			temp = paste0(temp, "\t")
			}
		}
	}

write(temp, "Campyl_SNPs_Table3.txt")
