#setwd("/Users/datn/n.dat/LmTag_revise_model")
#load("chr10_EAS_model.Rdata")
load("EAS_model.Rdata")

save(model, file = "model_linear.Rdata")

## variable selection

model <- lm(imputed_R2 ~ R2 , data = train)
round(summary(model)$r.squared*100,2)
save(model, file = "model0.1.Rdata")


model <- lm(imputed_R2 ~ tag_AF, data = train)
round(summary(model)$r.squared*100,2)
save(model, file = "model0.2.Rdata")


model <- lm(imputed_R2 ~ tagged_AF, data = train)
round(summary(model)$r.squared*100,2)
save(model, file = "model0.3.Rdata")

model <- lm(imputed_R2 ~ distance, data = train)
round(summary(model)$r.squared*100,2)
save(model, file = "model0.4.Rdata")



      
## interation models


model <- lm(imputed_R2 ~ R2*tag_AF + R2 + tag_AF + tagged_AF + distance, data = train)
round(summary(model)$r.squared*100,2)
save(model, file = "model1.1.Rdata")


model <- lm(imputed_R2 ~ R2*tagged_AF + R2 + tag_AF + tagged_AF + distance, data = train)
round(summary(model)$r.squared*100,2)
save(model, file = "model1.2.Rdata")


model <- lm(imputed_R2 ~ R2*distance + R2 + tag_AF + tagged_AF + distance, data = train)
round(summary(model)$r.squared*100,2)
save(model, file = "model1.3.Rdata")


#model1.4 <- lm(imputed_R2 ~ tag_AF*distance + R2 + tag_AF + tagged_AF + distance, data = train)
#round(summary(model1.4)


##################################################
# Polynomial regression
## R2
train2 = train

train2$R2_2 = train$R2^2
train2$R2_3 = train$R2^3
train2$R2_4 = train$R2^4


round(summary(model)

model <- lm(imputed_R2 ~ R2_2 + R2 + tag_AF + tagged_AF + distance, data = train2)
round(summary(model)$r.squared*100,2)
save(model, file = "model2.2.Rdata")


model <- lm(imputed_R2 ~ R2_3 + R2_2 + R2 + tag_AF + tagged_AF + distance, data = train2)
round(summary(model)$r.squared*100,2)
save(model, file = "model2.3.Rdata")

model <- lm(imputed_R2 ~ R2_4 + R2_3 + R2_2 + R2 + tag_AF + tagged_AF + distance, data = train2)
round(summary(model)$r.squared*100,2)
save(model, file = "model2.4.Rdata")



# Polynomial regression
## tag_AF
train3 = train

train3$tag_AF_2 = train$tag_AF^2
train3$tag_AF_3 = train$tag_AF^3
train3$tag_AF_4 = train$tag_AF^4




model <- lm(imputed_R2 ~ tag_AF_2 + R2 + tag_AF + tagged_AF + distance, data = train3)
round(summary(model)$r.squared*100,2)
save(model, file = "model3.2.Rdata")

model <- lm(imputed_R2 ~ tag_AF_3 + tag_AF_2 + R2 + tag_AF + tagged_AF + distance, data = train3)
round(summary(model)$r.squared*100,2)
save(model, file = "model3.3.Rdata")

model <- lm(imputed_R2 ~ tag_AF_4 + tag_AF_3 + tag_AF_2 + R2 + tag_AF + tagged_AF + distance, data = train3)
round(summary(model)$r.squared*100,2)
save(model, file = "model3.4.Rdata")


# Polynomial regression
## tagged_AF
train4 = train

train4$tagged_AF_2 = train$tagged_AF^2
train4$tagged_AF_3 = train$tagged_AF^3
train4$tagged_AF_4 = train$tagged_AF^4



model <- lm(imputed_R2 ~ tagged_AF_2 + R2 + tag_AF + tagged_AF + distance, data = train4)
round(summary(model)$r.squared*100,2)
save(model, file = "model4.2.Rdata")

model <- lm(imputed_R2 ~ tagged_AF_3 + tagged_AF_2 + R2 + tag_AF + tagged_AF + distance, data = train4)
round(summary(model)$r.squared*100,2)
save(model, file = "model4.3.Rdata")	

model <- lm(imputed_R2 ~ tagged_AF_4 + tagged_AF_3 + tagged_AF_2 + R2 + tag_AF + tagged_AF + distance, data = train4)
round(summary(model)$r.squared*100,2)
save(model, file = "model4.4.Rdata")


# Polynomial regression
## distance
train5 = train

train5$distance_2 = train$distance^2
train5$distance_3 = train$distance^3
train5$distance_4 = train$distance^4



model <- lm(imputed_R2 ~ distance_2 + R2 + tag_AF + tagged_AF + distance, data = train5)
round(summary(model)$r.squared*100,2)
save(model, file = "model5.2.Rdata")


model <- lm(imputed_R2 ~ distance_3 + distance_2 + R2 + tag_AF + tagged_AF + distance, data = train5)
round(summary(model)$r.squared*100,2)
save(model, file = "model5.3.Rdata")	

model <- lm(imputed_R2 ~ distance_4 + distance_3 + distance_2 + R2 + tag_AF + tagged_AF + distance, data = train5)
round(summary(model)$r.squared*100,2)
save(model, file = "model5.4.Rdata")





