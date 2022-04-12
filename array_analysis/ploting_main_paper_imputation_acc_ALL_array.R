###

library(ggplot2)

setwd("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/Array_analysis/Impute_acc")

#setwd("/media/datn/data/DatProjects/vn_array/4_paper/collected_all_res")

# Ploting imputation accuracy
load_res <- function(data){
	load(data)

	cutt_off = list(
	  c(0.01, 0.025),
	  c(0.025 ,0.05),
	  c(0.05, 0.075),
	  c(0.075, 0.125),
	  c(0.125, 0.2),
	  c(0.2, 0.3),
	  c(0.3, 0.4),
	  c(0.45, 0.5)
	)	

	#(0.01, 0.02],  (0.02, 0.03],  (0.03 ,0.04],  (0.04, 0.05],  (0.05, 0.075],  (0.075, 0.1],  (0.1, 0.125],  (0.125, 0.15],  (0.15, 0.2],  (0.2, 0.25],  (0.25, 0.3],  (0.3, 0.35],  (0.35, 0.4],  (0.4, 0.45],  (0.45, 0.5]	
	

	res = c()
	for (i in 1:length(cutt_off) ){
	  if(cutt_off[[i]][1] == 0.01){
	      pick = data$MAF >= cutt_off[[i]][1] &  data$MAF <= cutt_off[[i]][2]
	    }else{
	      pick = data$MAF > cutt_off[[i]][1] &  data$MAF <= cutt_off[[i]][2]
	    }
	  
	  tem = data[pick,]
	  tem_res = mean(tem$r_2, na.rm = T) * 100
	  names(tem_res) = paste( cutt_off[[i]][1],  cutt_off[[i]][2], sep = ":")
	  res = c(res, tem_res)
	}
	df = data.frame(cut_off = names(res), r_2 = res)
	return(df)
}




read_result <- function(array, pop){
	fn = paste0(array, "_", pop, "_chr10.Rdata")
	tem = load_res(fn)
	tem$array = array
	tem$pop = pop
	return(tem)
}

get_all_res <- function(array_list, pop_list){
	res = list()
	for(array in array_list){
		for(pop in pop_list){
			flag = paste0(array, ":", pop)
			res[[flag]] = read_result(array, pop)
		}
	}
	df = do.call("rbind", res)
	return(df)
}

pop_list = c("EAS", "EUR", "SAS", "KHV")
array_list = c("Affymetrix_6.0", "Axiom_GW_ASI", "Axiom_GW_EUR", "Axiom_JAPONICA", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB", "Infinium_GDA", "Infinium_GSA", "Multi-Ethnic_EUR_EAS_SAS", "LmTag_K200_24000", "LmTag_K200_28000", "LmTag_K200_32000" )


array_list = c("Affymetrix_6.0", "Axiom_GW_ASI", "Axiom_GW_EUR", "Axiom_JAPONICA", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB", "Infinium_GSA", "LmTag_K200_28000", "LmTag_K200_32000" )
#array_list = c("Affymetrix_6.0", "Axiom_GW_ASI", "Axiom_GW_EUR", "Axiom_JAPONICA", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB", "Infinium_GSA", "LmTag_K200_32000" )

pop_list = c("KHV")
array_list = c("Axiom_GW_ASI", "Axiom_GW_EUR", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB", "Infinium_GSA", "LmTag_K200_32000" )


df = get_all_res(array_list, pop_list)
df$r_2 = as.numeric(df$r_2)
df$pop[df$pop == "KHV"] = "VNP"
#df =x
df$TagSNPs = df$array
pick = grepl( "LmTag", df$array)
df$Types = "Arrays"
df$Types[pick] = "LmTag_K200"

df$SNP_arrays = df$TagSNPs
df$SNP_arrays [ df$TagSNPs == "LmTag_K200_32000"] = "VinGenChip"

library(ggplot2)
# p = ggplot(data= df, aes(x=cut_off, y=r_2, group = TagSNPs, col = TagSNPs)) + geom_line(aes(linetype= Types))+ geom_point() + guides(x = guide_axis(angle = 45)) + theme_light() + scale_y_continuous(breaks=seq(0,100,2.5)) + theme_light() + facet_wrap(~ pop, nrow = 3) + ylab("Imputation accuracy (%)") + xlab("MAF bins") 
#+ theme(legend.position="top")

library(ggplot2)
p = ggplot(data= df, aes(x=cut_off, y=r_2, group = SNP_arrays, col = SNP_arrays)) + geom_line() + geom_point() + guides(x = guide_axis(angle = 45)) + theme_light() + scale_y_continuous(breaks=seq(0,100,2.5)) + theme_light() + facet_wrap(~ pop, nrow = 3) + ylab("Imputation accuracy (%)") + xlab("MAF bins") 
#+ theme(legend.position="top")

pdf(file="Imputation_accuracy_chr10_ARRAYs_VingenChip.pdf",  width=6, height=4)
p
dev.off()



load_res_table <- function(data){
	load(data)

	cutt_off = list(
	  c(0.01, 0.025),
	  c(0.025 ,0.05),
	  c(0.05, 0.075),
	  c(0.075, 0.125),
	  c(0.125, 0.2),
	  c(0.2, 0.3),
	  c(0.3, 0.4),
	  c(0.45, 0.5)
	)	

	#(0.01, 0.02],  (0.02, 0.03],  (0.03 ,0.04],  (0.04, 0.05],  (0.05, 0.075],  (0.075, 0.1],  (0.1, 0.125],  (0.125, 0.15],  (0.15, 0.2],  (0.2, 0.25],  (0.25, 0.3],  (0.3, 0.35],  (0.35, 0.4],  (0.4, 0.45],  (0.45, 0.5]	
	

	res = c()
	for (i in 1:length(cutt_off) ){
	  if(cutt_off[[i]][1] == 0.01){
	      pick = data$MAF >= cutt_off[[i]][1] &  data$MAF <= cutt_off[[i]][2]
	    }else{
	      pick = data$MAF > cutt_off[[i]][1] &  data$MAF <= cutt_off[[i]][2]
	    }
	  
	  tem = data[pick,]
	  tem_res = mean(tem$r_2, na.rm = T)
	  names(tem_res) = paste( cutt_off[[i]][1],  cutt_off[[i]][2], sep = ":")
	  res = c(res, tem_res)
	}
	tem_res = mean(data$r_2, na.rm = T)
	names(tem_res) = "all_snp"
	res = c(res, tem_res)
	df = data.frame(cut_off = names(res), r_2 = res)
	return(df)
}




read_result_table <- function(array, pop){
	fn = paste0(array, "_", pop, "_chr10.Rdata")
	tem = load_res_table(fn)
	tem$array = array
	tem$pop = pop
	return(tem)
}

get_statistics <- function(array_list, pop_list){
	res = list()
	for(pop in pop_list){
		for(array in array_list){
			flag = paste0(pop, ":", array)
			tem = read_result_table(array, pop)
			tem = as.data.frame(tem)
			tem$r_2 = as.numeric(tem$r_2) *100
			tem$r_2 = round(tem$r_2, 2)
			tem2 = t(tem[,c(1,2)])
			tem2 = as.data.frame(tem2)
			tem2 = tem2[-1,]
			tem = tem2
			tem2$Population = pop
			tem2$Array = array
			tem = cbind(tem2[, c("Population", "Array")], tem)
			res[[flag]] = tem
		}
	}
	df = do.call("rbind", res)
	return(df)
}

df = get_statistics(array_list, pop_list)
df$Population[df$Population == "KHV"] = "VNP"

fwrite(df, "imp_acc_arrays_statistics.csv", sep = "\t", row.names = F)


#################

setwd("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/Array_analysis/array_tag_snp")
require(data.table)

array_list = c("Affymetrix_6.0", "Axiom_GW_ASI", "Axiom_GW_EUR", "Axiom_JAPONICA", "Axiom_PMDA", "Axiom_PMRA", "Axiom_UKB", "Infinium_GSA" )
# "LmTag_K200_32000"
fun = fread("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/GWAS_ClinVar_DB/GWAS_CLINVAR_ALL.txt")


get_stat <- function(array, fun){
	# path_res="/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/LmTag_K50_KHV.Rdata"
  fun_gwas = fun[fun$gwas_catalog == 1,]
  fun_cl = fun[fun$clinvar == 1,]
  df = fread(paste0(array, "_chr10.txt"), header = F)
  df$id = paste(df$V1, df$V2, sep = ":")
  n_tag = nrow(df)
  n_gwas = sum(df$id %in% fun_gwas$id)
  n_cl = sum(df$id %in% fun_cl$id)

  res = data.frame(array = array, n_tag = n_tag, n_gwas = n_gwas, n_cl = n_cl)
}


pop_list = c("EAS", "EUR", "SAS", "KHV")

get_stat_LmTag <- function(pop, fun){
  df = fread(paste0("LmTag_K200_", pop, "_32000_effect_score.txt"), header = T)	
  fun_gwas = fun[fun$gwas_catalog == 1,]
  fun_cl = fun[fun$clinvar == 1,]
  df$id = paste(df$chr, df$pos, sep = ":")
  n_tag = nrow(df)
  n_gwas = sum(df$id %in% fun_gwas$id)
  n_cl = sum(df$id %in% fun_cl$id)
  res = data.frame(array = paste0("LmTag_K200_", pop), n_tag = n_tag, n_gwas = n_gwas, n_cl = n_cl)
}


count_snp = list()
for(array in array_list){
	count_snp[[array]] = get_stat(array, fun)
}

for(pop in pop_list){
	count_snp[[pop]] = get_stat_LmTag(pop, fun)
}


df = do.call("rbind", count_snp)

fwrite(df, "count_chr10_array.csv", sep = "\t")










