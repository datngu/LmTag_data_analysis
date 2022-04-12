
setwd("functional_performace/top_32000_CADD")

library(ggplot2)

pop_list = c("SAS", "EAS", "EUR", "VNP")

methods = c("LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200", "EQ_uniform", "TagIt", "EQ_MAF", "FastTagger", "Baseline")

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
all$Methods = factor(all$K, levels = c("Baseline", "EQ_uniform", "EQ_MAF", "TagIt", "FastTagger", "LmTag_K1", "LmTag_K10", "LmTag_K20", "LmTag_K30", "LmTag_K50", "LmTag_K100", "LmTag_K200"))


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
