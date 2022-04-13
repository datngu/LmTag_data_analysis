
setwd("functional_performace/top_32000_CADD")

library(ggplot2)

pop_list = c("SAS", "EAS", "EUR", "VNP")

#methods = c("LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "LmTag_K500", "LmTag_K1000", "LmTag_K1500", "LmTag_K2000", "EQ_uniform", "TagIt", "EQ_MAF", "FastTagger", "Baseline")

methods =  c("Baseline", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "LmTag_K500", "LmTag_K1000", "LmTag_K1500", "LmTag_K2000")

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







#all$K <- factor(all$K, levels =c("EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200"))
all$Methods = factor(all$K, levels = c("Baseline", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "LmTag_K500", "LmTag_K1000", "LmTag_K1500", "LmTag_K2000"))


compute_p <- function(q){
  return(10^(-q/10))
}

#all$effect_score 
all$rank = (1 - compute_p(all$effect_score))*100

p = ggplot(all, aes(x = Methods, y = rank, color = Methods)) + geom_boxplot() + guides(x = guide_axis(angle = 0)) + theme_bw() + theme(legend.position = "none")   + stat_summary(fun = mean, geom = "point", col = "black") + facet_wrap(~ pop, nrow = 2) + theme_light() + ylab("CADD score percentile") + xlab("Methods") + theme(legend.position = "NONE") + guides(x = guide_axis(angle = 45)) + scale_y_continuous( breaks= seq(0,100,10))

#p = ggplot(all, aes(x = Methods, y = rank)) + geom_boxplot() + guides(x = guide_axis(angle = 0)) + theme_bw() + theme(legend.position = "none")   + stat_summary(fun = mean, geom = "point", col = "black") + facet_wrap(~ pop, nrow = 2) + theme_light() + ylab("Percentile ranking by CCAD scores") + xlab("Methods") + theme(legend.position = "NONE") + guides(x = guide_axis(angle = 45)) + scale_y_continuous( breaks= seq(0,100,10))

#aggregate(all$rank, list(all$Methods), FUN=mean)

pdf(file="../../Figure_2_CADD_score_percentile.pdf", width=8, height=7)
p
dev.off()



#setwd("/media/datn/data/1st_DC_sever/LmTag_GWAS_CLINVAR/tags_32000/tag_snp_with_cadd_score/")

pop_list = c("EAS", "EUR", "SAS", "KHV")

#methods = c("Control", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200")
methods =  c("Baseline", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "LmTag_K500", "LmTag_K1000", "LmTag_K1500", "LmTag_K2000")


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
