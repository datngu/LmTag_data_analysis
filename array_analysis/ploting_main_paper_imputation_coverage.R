###

library(ggplot2)

setwd("/media/datn/data/DatProjects/vn_array/4_paper/collected_all_res")


process_coverage <- function(data_path){
	load(data_path)
	res = compute_coverage(data)
	return(res)
}

compute_coverage <- function(df){
  #df <- df[df$MAF>=s & df$MAF <e,]
  c <- c()
  n <- c()
  for (i in c(1:100)){
    a <- i/100
    cover <- sum(df$r_2>=a, na.rm = T)/ nrow(df)
    c <- c(c,cover)
    n = c(n, a)
  }
  names(c) = n
  return(c)
}

read_pop <- function(pop, method_list = c("EQ", "Naive", "LmTag", "TagIt", "FastTagger")){

	EQ_path = paste0(pop, "_EQ_imputation_final_chr10.Rdata")
	Naive_path = paste0(pop, "_Naive_imputation_final_chr10.Rdata")
	LmTag_path = paste0(pop, "_LmTag_imputation_final_chr10.Rdata")
	TagIt_path = paste0(pop, "_TagIt_imputation_final_chr10.Rdata")
	FastTagger_path = paste0(pop, "_FastTagger_imputation_final_chr10.Rdata")

	for(m in method_list){
		if(m == "EQ") EQ = process_coverage(EQ_path)
		if(m == "Naive") Naive = process_coverage(Naive_path)
		if(m == "LmTag") LmTag = process_coverage(LmTag_path)
		if(m == "TagIt") TagIt = process_coverage(TagIt_path)
		if(m == "FastTagger") FastTagger = process_coverage(FastTagger_path)
	}

	plot_list = list()
	for(m in method_list){
		plot_list[[m]] = get(m)
	}

	all = data.frame()
		for(i in 1: length(plot_list)){
  	x = data.frame(x = names(plot_list[[i]]), y = plot_list[[i]] )
  	x$method  = names(plot_list[i])
  	all = rbind(all, x)
	}


	all$flag = grepl("IGSR", all$method)
	all$Reference_panel = "LOOCV_VN_504"
	all$Reference_panel[all$flag] = "IGSR"
	all$Methods = gsub("IGSR_", "", all$method)
	all$Pop = pop

	#p = ggplot(data=all, aes(x=MAF, y=r_2, group = Methods, col = Methods)) + geom_line()+ geom_point() + guides(x = guide_axis(angle = 45)) + theme_light() + scale_y_continuous(breaks=seq(0,1,0.025)) + theme_light() + expand_limits(y = c(0.5, 0.95)) + xlab("") + ylab("") + ggtitle(pop) + theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

  return(all)
}

pop_list = c("VNP", "EAS", "EUR", "AFR", "SAS", "AMR")
for(i in 1:length(pop_list)){
	#tem_df = data.frame()
	
	if(pop_list[i] == "AFR"){
		tem1 = read_pop("AFR", method_list = c("EQ", "Naive", "LmTag", "TagIt"))
	}else{
		tem1 = read_pop(pop_list[i])
	}
	
	if(i == 1){
		tem_df = tem1
	}else{
		tem_df = rbind(tem_df, tem1)
	}
}

tem_df$x = as.numeric(tem_df$x)

#p = ggplot(data=tem_df, aes(x = x, y = y, group = Methods, col = Methods)) + geom_line() + guides(x = guide_axis(angle = 45)) + theme_light() + scale_y_continuous(breaks=seq(0,1,0.025)) + theme_light() + facet_wrap(~ Pop, nrow = 3) + ylab("Imputation coverage") + xlab("Imputation cutoff")


p = ggplot(data = tem_df, aes( x = x, y = y , col = Methods)) + geom_line() + labs(x = "Imputation cutoff", y = "Imputation coverage") + scale_y_continuous(breaks=seq(0,1,0.05)) + scale_x_continuous(breaks=seq(0,1, 0.1)) + theme_light() + facet_wrap(~ Pop, nrow = 3)
# + theme(legend.position="top")



pdf(file="Imputation_coverage_chr10_all_pop.pdf", width=7, height=9)
p
dev.off()





# # Ploting imputation accuracy

# setwd("/media/datn/data/DatProjects/vn_array/4_paper/collected_all_res")

# pop = "VNP"


# EQ_path = paste0(pop, "_EQ_imputation_final_chr10.Rdata")
# LmTag_path = paste0(pop, "_LmTag_imputation_final_chr10.Rdata")
# TagIt_path = paste0(pop, "_TagIt_imputation_final_chr10.Rdata")
# FastTagger_path = paste0(pop, "_FastTagger_imputation_final_chr10.Rdata")



# #Naive_IGSR_path = "/media/datn/data/DatProjects/vn_array/4_paper/p_VN504/IGSR_2504_imputation_res/naiveimputed_2504_IGSR.Rdata"
# EQ_IGSR_path = paste0("IGSR_", pop, "_EQ_imputation_final_chr10.Rdata")
# LmTag_IGSR_path = paste0("IGSR_", pop, "_LmTag_imputation_final_chr10.Rdata")
# TagIt_IGSR_path = paste0("IGSR_", pop, "_TagIt_imputation_final_chr10.Rdata")
# FastTagger_IGSR_path = paste0("IGSR_", pop, "_FastTagger_imputation_final_chr10.Rdata")

# #Naive = load_res(Naive_path)
# EQ = load_res(EQ_path)
# LmTag = load_res(LmTag_path)
# TagIt = load_res(TagIt_path)
# FastTagger = load_res(FastTagger_path)


# EQ_IGSR = load_res(EQ_IGSR_path)
# LmTag_IGSR = load_res(LmTag_IGSR_path)
# TagIt_IGSR = load_res(TagIt_IGSR_path)
# FastTagger_IGSR = load_res(FastTagger_IGSR_path)


# #plot_list = list( Naive = Naive, LmTag = LmTag, TagIt = TagIt, Hybrid = Hybrid, FastTagger = FastTagger, a_IGSR_LmTag = LmTag_IGSR, a_IGSR_TagIt = TagIt_IGSR,  a_IGSR_Hybrid = Hybrid_IGSR, a_IGSR_FastTagger = FastTagger_IGSR, a_IGSR_Naive = Naive_IGSR)

# plot_list = list( LmTag = LmTag, TagIt = TagIt, FastTagger = FastTagger, EQ = EQ, IGSR_LmTag = LmTag_IGSR, IGSR_TagIt = TagIt_IGSR, IGSR_FastTagger = FastTagger_IGSR, IGSR_EQ = EQ_IGSR)

# all = data.frame()
# for(i in 1: length(plot_list)){
#   x = data.frame(MAF = names(plot_list[[i]]), r_2 = plot_list[[i]] )
#   x$method  = names(plot_list[i])
#   all = rbind(all, x)
# }

# all$flag = grepl("IGSR", all$method)
# all$Reference_panel = "LOOCV_VN_504"
# all$Reference_panel[all$flag] = "IGSR"
# all$Methods = gsub("IGSR_", "", all$method)
# # all$flag = as.factor(all$flag)
# # all$flag[ all$flag == "FALSE"] = "LOOCV_VN_504"
# # all$flag[ all$flag == "TRUE"] = "IGSR_panel"#


# library(ggplot2)
# P1 = ggplot(data = all, aes(x = MAF, y=r_2, group = method, col = Methods)) + geom_line(aes(linetype = Reference_panel)) +
#   geom_point() + guides(x = guide_axis(angle = 45)) + theme_light() + scale_y_continuous(breaks=seq(0,1,0.025)) + theme_light() + expand_limits(y = c(0.6, 1))

# ## nohup bash ./bash.sh > log.txt 2>&1

