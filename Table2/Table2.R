# This script was executed in the parent directory of the directory containing the TSV files of SNP results ("Table2_data").


curr_path = getwd()
file_names = list.files("./Table2_data")

serovars = c("Schwarzengrund","Infantis","Manhattan","Heidelberg","Enteritidis","Blockley","Minnesota","O4i-","Thompson","Typhimurium","Saintpaul","Newport","Stanley","Agona")
source_pair = c("human_vs_human", "human_vs_food", "human_vs_animal")

# An object that stores results, where each serovar is added as a row in the loop below.
temp = paste0("Serovar\t", "human- vs human-derived isolates\t", "human- vs food-derived isolates\t", "human- vs animal-derived isolates\n")

# Writing the median and 25th/75th percentiles for each sectoral combination within each serovar
for (i in serovars){
	temp = paste0(temp, i, "\t")
	file_name_sero = file_names[grep(i, file_names)]
	for (j in source_pair){
		target_file = file_name_sero[grep(j,file_name_sero)]
		if (length(target_file) == 0){
			temp = paste0(temp,"-")
			}
		else{
			df = read.table(paste0(curr_path, "/Table2_data/", target_file), header=T, sep="\t",comment.char="",quote ='"')
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

write(temp, "Sal_SNPs_Table2.txt")
