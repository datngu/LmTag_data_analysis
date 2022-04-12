
setwd("tagging_performance")

library(ggplot2)

compute_mean_r2 <- function(path){
	load(path)
	return(mean(data$r_2))
}

mean_r2 <- function( method, pop){
	if(pop == "EUR") cutoff = c(8000, 12000, 16000, 20000, 24000, 28000, 32000, 36000, 40000, 41114)
	if(pop == "SAS") cutoff = c(8000, 12000, 16000, 20000, 24000, 28000, 32000, 36000, 40000, 44000, 44866)
	if(pop == "EAS") cutoff = c(8000, 12000, 16000, 20000, 24000, 28000, 32000, 32970)
	if(pop == "VNP") cutoff = c(8000, 12000, 16000, 20000, 24000, 28000, 32000, 32067)
	res = c()
	for( i in cutoff){
		file = paste0( method, "_",  pop, "_", i, ".Rdata")
		r2 = compute_mean_r2(file)
		res = c(res, r2)
	}
	df = data.frame(cutoff = cutoff, mean_r2  = res, Methods = method, Pop = pop)
	return(df)
}



read_pop <- function( pop_list = c("VNP", "EAS", "EUR", "SAS"), method_list = c("EQ_MAF", "EQ_uniform", "LmTag", "TagIt", "FastTagger"))
{
	df = NULL
	for( p in pop_list){
		for (m in method_list){
			tem = mean_r2(m, p)
			if(is.null(df)){
				df = tem}
			else{
				df = rbind(df, tem)
			}
		}
	}
	return(df)
}


get_cut_stat <- function(cutoff, df2){
	methods = c("LmTag_K200","EQ_uniform", "EQ_MAF", "TagIt", "FastTagger")
	df2 = df2[df2$cutoff == cutoff,]

	l = list()
	for( p in c("EAS", "EUR", "SAS", "VNP")){
	  tem = df2[df2$Pop == p,]
	  order = match(methods, tem$Methods)
	  tem = tem[order,]
	  l[[p]] = tem
	}

	res = do.call("cbind", l)
}


tem_df = read_pop()
tem_df$mean_r2 = round(tem_df$mean_r2*100,2)
tem_df$Methods[tem_df$Methods == "LmTag"] = "LmTag_K200"


## PLOTING

#cutoffs = rev(c(8000, 12000, 16000, 20000, 24000, 28000, 32000))
breaks = c(8000, 12000, 16000, 20000, 24000, 28000, 32000)
p = ggplot(data = tem_df[tem_df$cutoff <= max(breaks),], aes( x = cutoff, y = mean_r2  , col = Methods)) + geom_line() + labs(x = "Number of top Tag SNP", y = "Mean imputation accuracy of all SNPs (%)") + scale_y_continuous(breaks = seq(0,100,5))  + theme_light() + facet_wrap(~ Pop, nrow = 2) + scale_x_continuous(breaks = breaks) + guides(x = guide_axis(angle = 45))
# + theme(legend.position="top")

out = "../Figure_4_imputation_accuracy.pdf"
pdf(file=out,  width=8, height=6.5)
p
dev.off()


#### GETTING STATISTIC TABLES

cutoffs = rev(c(8000, 12000, 16000, 20000, 24000, 28000, 32000))


l = list()
for( c in cutoffs){
	l[[ as.character(c)]] = get_cut_stat(c, tem_df)
}

df2 = do.call("rbind", l)

write.csv(df2, file = "../Figure_4_imputation_accuracy_statistics.csv")
