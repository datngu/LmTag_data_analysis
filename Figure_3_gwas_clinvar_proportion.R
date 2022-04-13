
setwd("functional_performace/top_32000_CADD")

library(ggplot2)
library(data.table)

fun = fread("../functional_scoring/GWAS_CLINVAR_ALL.txt")

cout_fun_snp <- function(fun, pop, method){
  path = paste0(method, "_", pop, "_32000_effect_score.txt")
  path_res = paste0(method, "_", pop, "_32000.Rdata")
  # path ="/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/LmTag_K50_KHV_32000_effect_score.txt"
  # path_res="/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/LmTag_K50_KHV.Rdata"
  fun_gwas = fun[fun$gwas_catalog == 1,]
  fun_cl = fun[fun$clinvar == 1,]
 
  tag = read.table(path, header = T)
  
  load(path_res)
  data$id = paste(data$X.CHROM, data$POS, sep = ":")
  tag$id = paste(tag$chr, tag$pos, sep = ":")
  all_count = nrow(data)
  all_tag = nrow(tag)
  fun_all_count = sum(data$id %in% fun$id)
  fun_gwas_count = sum(data$id %in% fun_gwas$id)
  fun_cl_count = sum(data$id %in% fun_cl$id)
  tag_fun_count = sum(tag$id %in% fun$id)
  tag_gwas_count = sum(tag$id %in% fun_gwas$id)
  tag_cl_count = sum(tag$id %in% fun_cl$id)
  res = data.frame(pop = pop, method = method , all_count = all_count, all_tag = all_tag, fun_all_count = fun_all_count, fun_gwas_count = fun_gwas_count, fun_cl_count = fun_cl_count, tag_fun_count = tag_fun_count, tag_gwas_count = tag_gwas_count, tag_cl_count = tag_cl_count)
  return(res)
}

#cout_fun_snp(fun, "KHV", "TagIt")
pop_list = c("SAS", "EAS", "EUR", "VNP")

#methods = c("EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "Baseline")
methods =  c("Baseline", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "LmTag_K500", "LmTag_K1000", "LmTag_K1500", "LmTag_K2000")

res = NULL
for(pop in pop_list){
  for(k in methods){
    tem = cout_fun_snp(fun, pop, k)
    if(is.null(res)){
      res = tem
    }else{
      res = rbind(res, tem)
    }
  }
}

##################################################################

df = res

size = 32000

df$pc_fun_all = (df$tag_fun_count / df$fun_all_count * 100) / df$all_tag * size
df$pc_gwas = (df$tag_gwas_count / df$fun_gwas_count * 100) / df$all_tag * size
df$pc_cl = (df$tag_cl_count / df$fun_cl_count * 100) / df$all_tag * size

# methods = c("Control", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200")

df2 = df[, c("method", "pop", "tag_gwas_count", "tag_cl_count", "pc_gwas", "pc_cl")]

l = list()
for( p in c("EAS", "EUR", "SAS", "KHV")){
  tem = df2[df2$pop == p,]
  order = match(methods, tem$method)
  tem = tem[order,]
  l[[p]] = tem
}

df3 = do.call("cbind", l)

write.csv(df2, file = "percentage_GWAS_Clinvar.csv")

df1 = df[, c("pop", "method", "pc_fun_all")]
colnames(df1) = c("pop", "method", "Percentages")
df1$Databases = "Clinvar + GWAS"
df2 = df[, c("pop", "method", "pc_gwas")]
colnames(df2) = c("pop", "method", "Percentages")
df2$Databases = "GWAS"
df3 = df[, c("pop", "method", "pc_cl")]
colnames(df3) = c("pop", "method", "Percentages")
df3$Databases = "Clinvar"

df = rbind(df1,df2,df3)
df = rbind(df2,df3)

#df$Methods = factor(df$method, levels = c("Baseline", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200"))
df$Methods = factor(df$method,  levels = c("Baseline", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "LmTag_K500", "LmTag_K1000", "LmTag_K1500", "LmTag_K2000"))


#df$pop[df$pop == "KHV"] = "VNP"

p = ggplot(data=df, aes(x= Methods, y= Percentages, fill= Databases)) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_brewer(palette="Paired") + theme_bw() + facet_wrap(~ pop, nrow = 2) + guides(x = guide_axis(angle = 45))

# + theme_minimal()


pdf(file="../../Figure_3_gwas_clinvar_proportion.pdf", width=8, height=7)
p
dev.off()






# counting clinvar and gwas catalog SNPs

library(data.table)
library(dplyr)

#setwd("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/")

fun = fread("/media/datn/data/db_VingenChip/GWAS_CLINVAR_ALL.txt")

cout_fun_snp <- function(fun, pop, method){
  path = paste0(method, "_", pop, "_32000_effect_score.txt")
  path_res = paste0(method, "_", pop, "_32000.Rdata")
  # path ="/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/LmTag_K50_KHV_32000_effect_score.txt"
  # path_res="/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/LmTag_K50_KHV.Rdata"
  fun_gwas = fun[fun$gwas_catalog == 1,]
  fun_cl = fun[fun$clinvar == 1,]
 
  tag = read.table(path, header = T)
  
  load(path_res)
  data$id = paste(data$X.CHROM, data$POS, sep = ":")
  tag$id = paste(tag$chr, tag$pos, sep = ":")
  all_count = nrow(data)
  all_tag = nrow(tag)
  fun_all_count = sum(data$id %in% fun$id)
  fun_gwas_count = sum(data$id %in% fun_gwas$id)
  fun_cl_count = sum(data$id %in% fun_cl$id)
  tag_fun_count = sum(tag$id %in% fun$id)
  tag_gwas_count = sum(tag$id %in% fun_gwas$id)
  tag_cl_count = sum(tag$id %in% fun_cl$id)
  res = data.frame(pop = pop, method = method , all_count = all_count, all_tag = all_tag, fun_all_count = fun_all_count, fun_gwas_count = fun_gwas_count, fun_cl_count = fun_cl_count, tag_fun_count = tag_fun_count, tag_gwas_count = tag_gwas_count, tag_cl_count = tag_cl_count)
  return(res)
}

#cout_fun_snp(fun, "KHV", "TagIt")
pop_list = c("SAS", "EAS", "EUR", "KHV")

#methods = c("EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "Control")

res = NULL
for(pop in pop_list){
  for(k in methods){
    tem = cout_fun_snp(fun, pop, k)
    if(is.null(res)){
      res = tem
    }else{
      res = rbind(res, tem)
    }
  }
}

##################################################################

df = res
df_bk = res

size = 32000

df$pc_fun_all = (df$tag_fun_count / df$fun_all_count * 100) / df$all_tag * size
df$pc_gwas = (df$tag_gwas_count / df$fun_gwas_count * 100) / df$all_tag * size
df$pc_cl = (df$tag_cl_count / df$fun_cl_count * 100) / df$all_tag * size

methods = c("Control", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200")

df2 = df[, c("method", "pop", "tag_gwas_count", "tag_cl_count", "pc_gwas", "pc_cl")]

l = list()
for( p in c("EAS", "EUR", "SAS", "KHV")){
  tem = df2[df2$pop == p,]
  order = match(methods, tem$method)
  tem = tem[order,]
  l[[p]] = tem
}

df3 = do.call("cbind", l)

write.csv(df2, file = "percentage_GWAS_Clinvar.csv")