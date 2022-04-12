setwd("model_optimization")

model_list = c('model_linear.Rdata', 'model0.1.Rdata', 'model0.2.Rdata', 'model0.3.Rdata', 'model0.4.Rdata', 'model1.1.Rdata', 'model1.2.Rdata', 'model1.3.Rdata', 'model2.2.Rdata', 'model2.3.Rdata', 'model2.4.Rdata', 'model3.2.Rdata', 'model3.3.Rdata', 'model3.4.Rdata', 'model4.2.Rdata', 'model4.3.Rdata', 'model4.4.Rdata', 'model5.2.Rdata', 'model5.3.Rdata', 'model5.4.Rdata')

res = data.frame(model = model_list)
res$model_R_square = 0.0
res$imputation_r2 = 0.0

load_model <- function(model_name){
	load(model_name)
	return(round(summary(model)$r.squared*100,2))
}

compute_mean_r2 <- function(model_name){
	fn = paste0(model_name, "_imputed_accuracy.Rdata")
	load(fn)
	return(mean(data$r_2))
}

for(i in 1:length(model_list) ){
	model_name = model_list[i]
	res$model_R_square[i] = load_model(model_name)
	res$imputation_r2[i] = compute_mean_r2(model_name)
}

write.table(res, file = "imputation_r2_all_models.txt", sep = "\t", quote = F, row.names = F)
