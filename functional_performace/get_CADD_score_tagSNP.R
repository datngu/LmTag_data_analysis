#!/usr/bin/env Rscript

setwd("/functional_performace")

get_cadd <- function(cadd, tag, outname){
	t = read.table(tag, header = F)
	names(t) = c("chr", "pos")
	c = read.table(cadd, header = T)	
	t$id = paste(t[,1], t[,2], sep = ":")
	c$id = paste(c[,1], c[,2], sep = ":")	
	t$effect_score = c$eff_score[match(t$id, c$id)]	
	write.table(t, file = outname, sep = "\t", quote = F, row.names = F)
}

get_cadd_LmTag <- function(tag, outname){
	t = read.table(tag, header = T)
	t = t[, c("chr","pos", "id", "effect_score")]
	write.table(t, file = outname, sep = "\t", quote = F, row.names = F)
}



# LmTag
pop = c("VNP", "EAS", "SAS", "EUR")
K_list = c(1,5,10,20,30,50,100,200)

for( p in pop ){
	for(k in K_list){
		#cadd = paste0("/media/datn/data/1st_DC_sever/res_tagSNP/CADD_scores/mapped_", p, "_chr10_CCAD.txt")
		#tag =  paste0("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/LmTag_32000_raw/", m , "_", p, "_32000.tag.txt")
		tag =  paste0("top_32000_raw/LmTag_K", k ,"_", p, "_32000.tag.txt")
		out = paste0("top_32000_CADD/LmTag_K", k , "_", p , "_32000_effect_score.txt")
		get_cadd_LmTag(tag, out )
	}
}


pop = c("VNP", "EAS", "SAS", "EUR")
K_list = c(500, 1000, 1500, 2000)

for( p in pop ){
	for(k in K_list){
		#cadd = paste0("/media/datn/data/1st_DC_sever/res_tagSNP/CADD_scores/mapped_", p, "_chr10_CCAD.txt")
		#tag =  paste0("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/LmTag_32000_raw/", m , "_", p, "_32000.tag.txt")
		tag =  paste0("top_32000_raw/LmTag_K", k ,"_", p, "_32000.tag.txt")
		out = paste0("top_32000_CADD/LmTag_K", k , "_", p , "_32000_effect_score.txt")
		get_cadd_LmTag(tag, out )
	}
}

# All method
pop = c("VNP", "EAS", "SAS", "EUR")
method = c("TagIt", "EQ_MAF", "EQ_uniform", "FastTagger")

for( p in pop ){
	for(m in method){
		cadd = paste0("functional_scoring/mapped_", p, "_chr10_CCAD.txt")
		tag =  paste0("top_32000_raw/", m , "_", p, "_32000.tag.txt")
		out = paste0("top_32000_CADD/", m , "_", p , "_32000_effect_score.txt")
		get_cadd(cadd, tag, out )
	}
}



## Get CADD score - Baseline - all input SNPs

get_cadd_control <- function(cadd, impute_res, outname, out_Rdata){
	load(impute_res)
	t = data[,c(1,2)]
	data$r_2 = 1
	names(t) = c("chr", "pos")
	c = read.table(cadd, header = T)	
	t$id = paste(t[,1], t[,2], sep = ":")
	c$id = paste(c[,1], c[,2], sep = ":")	
	t$effect_score = c$eff_score[match(t$id, c$id)]	
	write.table(t, file = outname, sep = "\t", quote = F, row.names = F)
	save(data, file = out_Rdata)
}

pop = c("VNP", "EAS", "SAS", "EUR")

for( p in pop ){
	cadd = paste0("functional_scoring/mapped_", p, "_chr10_CCAD.txt")
	impute_res =  paste0("top_32000_CADD/EQ_MAF_", p, "_32000.Rdata")
	out = paste0("top_32000_CADD/", "Baseline_", p , "_32000_effect_score.txt")
	out_Rdata = paste0("top_32000_CADD/", "Baseline_", p , "_32000.Rdata")
	get_cadd_control(cadd, impute_res, out, out_Rdata )
}








