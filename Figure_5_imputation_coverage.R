
setwd("tagging_performance")

library(ggplot2)



compute_coverage <- function(path, imp_cutoff = 0.8){
	load(path)
	cover = (sum(data$r_2>= imp_cutoff)/ nrow(data))
	return(cover)
}

get_coverage <- function( method, pop){
	if(pop == "EUR") cutoff = c(8000, 12000, 16000, 20000, 24000, 28000, 32000, 36000, 40000, 41114)
	if(pop == "SAS") cutoff = c(8000, 12000, 16000, 20000, 24000, 28000, 32000, 36000, 40000, 44000, 44866)
	if(pop == "EAS") cutoff = c(8000, 12000, 16000, 20000, 24000, 28000, 32000, 32970)
	if(pop == "VNP") cutoff = c(8000, 12000, 16000, 20000, 24000, 28000, 32000, 32067)
	res = c()
	for( i in cutoff){
		file = paste0( method, "_",  pop, "_", i, ".Rdata")
		r2 = compute_coverage(file)
		res = c(res, r2)
	}
	df = data.frame(cutoff = cutoff, coverage  = res, Methods = method, Pop = pop)
	return(df)
}



read_pop <- function( pop_list = c("VNP", "EAS", "EUR", "SAS"), method_list = c("EQ_MAF", "EQ_uniform", "LmTag", "TagIt", "FastTagger"))
{
	df = NULL
	for( p in pop_list){
		for (m in method_list){
			tem = get_coverage(m, p)
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
tem_df$coverage = tem_df$coverage*100
breaks = c(8000, 12000, 16000, 20000, 24000, 28000, 32000)
tem_df$Methods[tem_df$Methods == "LmTag"] = "LmTag_K200"

## PLOTING

p = ggplot(data = tem_df[tem_df$cutoff <= max(breaks),], aes( x = cutoff, y = coverage  , col = Methods)) + geom_line() + labs(x = "Number of top Tag SNP", y = "Imputation coverage (%)") + scale_y_continuous(breaks = seq(0,100,5))  + theme_light() + facet_wrap(~ Pop, nrow = 2) + scale_x_continuous(breaks = breaks) + guides(x = guide_axis(angle = 45))

out = "../Figure_5_imputation_coverage.pdf"
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

write.csv(df2, file = "../Figure_5_imputation_coverage_statistics.csv")
