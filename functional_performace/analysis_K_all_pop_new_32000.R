


setwd("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/")
pop_list = c("SAS", "EAS", "EUR", "KHV")

methods = c("LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "EQ_uniform", "TagIt", "EQ_MAF", "FastTagger", "Control")

#methods ="EQ"
for(pop in pop_list){
  for(k in methods){
    path = paste0(k, "_", pop, "_32000_effect_score.txt")
    tem = read.table(path, header = T)
    if(nrow(tem) < 2) next
    tem$K = as.character(k)
    tem$pop = pop
    if(pop == pop_list[1] && k == methods[1]){
      all = tem
    }else{
      all = rbind(all, tem)
    }
  }
}






library(ggplot2)
#all$K <- factor(all$K, levels =c("EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200"))
all$Methods = factor(all$K, levels = c("Control", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200"))

all$pop[all$pop == "KHV"] = "VNP"

compute_p <- function(q){
  return(10^(-q/10))
}

#all$effect_score 
all$rank = (1 - compute_p(all$effect_score))*100

p = ggplot(all, aes(x = Methods, y = rank, color = Methods)) + geom_boxplot() + guides(x = guide_axis(angle = 0)) + theme_bw() + theme(legend.position = "none")   + stat_summary(fun = mean, geom = "point", col = "black") + facet_wrap(~ pop, nrow = 2) + theme_light() + ylab("Percentile ranking by CCAD scores") + xlab("Methods") + theme(legend.position = "NONE") + guides(x = guide_axis(angle = 45)) + scale_y_continuous( breaks= seq(0,100,10))

#p = ggplot(all, aes(x = Methods, y = rank)) + geom_boxplot() + guides(x = guide_axis(angle = 0)) + theme_bw() + theme(legend.position = "none")   + stat_summary(fun = mean, geom = "point", col = "black") + facet_wrap(~ pop, nrow = 2) + theme_light() + ylab("Percentile ranking by CCAD scores") + xlab("Methods") + theme(legend.position = "NONE") + guides(x = guide_axis(angle = 45)) + scale_y_continuous( breaks= seq(0,100,10))

#aggregate(all$rank, list(all$Methods), FUN=mean)

pdf(file="Percentile_rank_by_K_values.pdf", width=8, height=7)
p
dev.off()



####

# SHELL CODES
# in=/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/impute_acc

# out=/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score

# for p in EAS SAS EUR KHV
# do
#     cp ${in}/LmTag_${p}_K1.Rdata $out/LmTag_K1_${p}.Rdata
#     cp ${in}/LmTag_${p}_K10.Rdata $out/LmTag_K10_${p}.Rdata
#     cp ${in}/LmTag_${p}_K20.Rdata $out/LmTag_K20_${p}.Rdata
#     cp ${in}/LmTag_${p}_K30.Rdata $out/LmTag_K30_${p}.Rdata
#     cp ${in}/LmTag_${p}_K50.Rdata $out/LmTag_K50_${p}.Rdata
#     cp ${in}/LmTag_${p}_K100.Rdata $out/LmTag_K100_${p}.Rdata
#     cp ${in}/LmTag_${p}_K200.Rdata $out/LmTag_K200_${p}.Rdata
# done


# for m in EQ_MAF EQ_uniform TagIt FastTagger
# do
#     cp ${in}/${m}_KHV_32067.Rdata $out/${m}_KHV.Rdata
#     cp ${in}/${m}_SAS_44866.Rdata $out/${m}_SAS.Rdata
#     cp ${in}/${m}_EUR_41114.Rdata $out/${m}_EUR.Rdata
#     cp ${in}/${m}_EAS_32970.Rdata $out/${m}_EAS.Rdata
# done

setwd("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/")

pop_list = c("EAS", "EUR", "SAS", "KHV")

methods = c("Control", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200")

n = length(pop_list)
m = length(methods)
mat = matrix(0, nrow = m, ncol = n)

number_tag = as.data.frame(mat)
colnames(number_tag) = pop_list
row.names(number_tag) = as.character(methods)


imputed_acc = total_snp = mean_functional_score = prob_score = number_tag

load_res = function(path){
  load(path)
  res = mean(data$r_2)
  return(res)
}

compute_p <- function(q){
  return(10^(-q/10))
}



#methods ="EQ"
for(pop in pop_list){
  for(k in methods){
    path = paste0(k, "_", pop, "_32000_effect_score.txt")
    tem = read.table(path, header = T)
    x = nrow(tem)
    number_tag[as.character(k),pop] = x
    tem$pctl = (1 - compute_p(tem$effect_score))*100
    mean_functional_score[as.character(k),pop] = mean(tem$effect_score)
    prob_score[as.character(k),pop] = mean(tem$pctl)
    path_res = paste0(k, "_", pop, "_32000.Rdata")
    load(path_res)
    imputed_acc[as.character(k),pop] = mean(data$r_2)
    total_snp[as.character(k),pop] = nrow(data)
  }
}





write.csv(imputed_acc, file = "imputed_acc.csv")
write.csv(total_snp, file = "total_snp.csv")
write.csv(prob_score, file = "pctl_functional_score.csv")
write.csv(mean_functional_score, file = "mean_CADD_score.csv")
write.csv(number_tag, file = "number_tag.csv")


# counting clinvar and gwas catalog SNPs

library(data.table)
library(dplyr)

setwd("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/")

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

methods = c("EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "Control")

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

df$Methods = factor(df$method, levels = c("Control", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200"))


df$pop[df$pop == "KHV"] = "VNP"

p = ggplot(data=df, aes(x= Methods, y= Percentages, fill= Databases)) + geom_bar(stat="identity", position=position_dodge()) + scale_fill_brewer(palette="Paired") + theme_bw() + facet_wrap(~ pop, nrow = 2) + guides(x = guide_axis(angle = 45))

# + theme_minimal()







pdf(file="Percentages_funcional_markers.pdf", width=8, height=7)
p
dev.off()




