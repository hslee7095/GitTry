library(DiceKriging)
library(RobustGaSP)
library(lhs)

#citation(package = "RobustGaSP", lib.loc = NULL)


set.seed(1)
input <- 10*maximinLHS(n=15, k=1) 
#input=10*seq(0,1,1/9)
output<-higdon.1.data(input)
model<- rgasp(design = input, response = output)
model

plot(model)

testing_input = as.matrix(seq(0,10,1/100))
model.predict<-predict(model,testing_input)
names(model.predict)

#########plot predictive distribution
testing_output=higdon.1.data(testing_input)
plot(testing_input,model.predict$mean,type='l',col='blue',
     xlab='input',ylab='output')
polygon( c(testing_input,rev(testing_input)),c(model.predict$lower95,
     rev(model.predict$upper95)),col = "grey80", border = F)
lines(testing_input, testing_output)
lines(testing_input,model.predict$mean,type='l',col='blue')
lines(input, output,type='p')


##########plot sampling distribution
model.sample=simulate(model,testing_input,num_sample=10)
matplot(testing_input,model.sample, type='l',xlab='input',ylab='output')
lines(input,output,type='p')


##########Borehole examples to test findInertInputs()
set.seed(1)
input <- maximinLHS(n=40, k=8)  # maximin lhd sample
# rescale the design to the domain of the Borehole function
LB=c(0.05,100,63070,990,63.1,700,1120,9855)
UB=c(0.15,50000,115600,1110,116,820,1680,12045)
range=UB-LB
for(i in 1:8){
  input[,i]=LB[i]+range[i]*input[,i]
}
num_obs=dim(input)[1]
output=matrix(0,num_obs,1)
for(i in 1:num_obs){
  output[i]=borehole(input[i,])
}
m<- rgasp(design = input, response = output, lower_bound=FALSE)
P=findInertInputs(m)

##plot the Borehole one input
mean_value=(UB+LB)/2
record_value=matrix(0,8,100)
##the boundary
for(i in 1:8){
  input_i=seq(LB[i],UB[i],(UB[i]-LB[i])/99)
  for(j in 1:100){
    plot_input=mean_value
    plot_input[i]=input_i[j]
    record_value[i,j]=borehole(plot_input)
  }
}

max_val=max(record_value)
min_val=min(record_value)
pdf(file="borehole.pdf",height=7,width=7)
par(mfrow=c(2,4))
for(i in 1:8){
  plot(seq(LB[i],UB[i],(UB[i]-LB[i])/99),record_value[i,],
       type='l',ylim=c(min_val,max_val),cex=1.5,cex.axis=1.5,cex.lab=1.5,
       xlab=paste('input',i),ylab='output')
}
dev.off()
#######

m.subset<- rgasp(design = input[,c(1,4,6,7,8)], response = output, nugget.est=T)

m.subset
###
###

m.full<- rgasp(design = input, response = output)
m.subset<- rgasp(design = input[,c(1,4,6,7,8)], response = output, nugget.est=T)
dk.full<- km(design = input, response = output)
dk.subset<- km(design = input[,c(1,4,6,7,8)], response = output,nugget.estim=T)


#####create testing
set.seed(1)
dim_inputs=dim(input)[2]
num_testing_input <- 100    
testing_input <- matrix(runif(num_testing_input*dim_inputs),
                        num_testing_input,dim_inputs)
# resale the points to the region to be predict
for(i in 1:8){
  testing_input[,i]=LB[i]+range[i]*testing_input[,i]
}

# Perform prediction
m.full.predict<-predict(m.full, testing_input)
m.subset.predict<-predict(m.subset, testing_input[,c(1,4,6,7,8)])
dk.full.predict<-predict(dk.full, newdata = testing_input,type = 'UK')
dk.subset.predict<-predict(dk.subset, 
                           newdata = testing_input[,c(1,4,6,7,8)],type = 'UK')


# generate real outputs
testing_output <- matrix(0,num_testing_input,1)
for(i in 1:num_testing_input){
  testing_output[i]<-borehole(testing_input[i,])
}
#see the prediction mean squared error
m.full.error=abs(m.full.predict$mean-testing_output)
m.subset.error=abs(m.subset.predict$mean-testing_output)
dk.full.error=abs(dk.full.predict$mean-testing_output)
dk.subset.error=abs(dk.subset.predict$mean-testing_output)

###plot the error
# plot(m.full.error,pch=0,col='blue',ylim=c(0,max(m.full.error,m.subset.error,
#      dk.full.error,dk.subset.error)),xlab='number',ylab='absolute error' )
# lines(m.subset.error,type='p',pch=1,col='green')
# lines(dk.full.error,type='p',pch=2,col='red')
# lines(dk.subset.error,type='p',pch=3,col='brown')
#pdf(file='abs_errors_Borehole.pdf',height=5,width=7)
record_error=cbind(m.full.error,m.subset.error,dk.full.error,dk.subset.error)
colnames(record_error)=c('RobustGasP \n full set','5','7','8')
par(mgp = c(3, 1.8,0))#mgp sets position of axis label, tick labels and axis
boxplot(record_error,cex.axis=1.1,cex.lab=1.1,
        names=c('RobustGaSP \n full set','RobustGaSP \n  subset',
                'DiceKriging \n full set','DiceKriging \n  subset'),
        ylab='absolute error')
#dev.off()


RMSE.m.full <- sqrt(sum((m.full.predict$mean-testing_output)^2)/(num_testing_input))  
RMSE.m.subset <- sqrt(sum((m.subset.predict$mean-testing_output)^2)/(num_testing_input))  

RMSE.dk.full <- sqrt(sum((dk.full.predict$mean-testing_output)^2)/(num_testing_input))  
RMSE.dk.subset <- sqrt(sum((dk.subset.predict$mean-testing_output)^2)/(num_testing_input))  

RMSE.m.full
RMSE.m.subset
RMSE.dk.full
RMSE.dk.subset



#####################examples

sinewave<-function(x){
  3*sin(5*pi*x)*x+cos(7*pi*x)
}


#########12 design points
input=as.matrix(seq(0,1,1/11))
output=sinewave(input)

m<-rgasp(design=input, response=output)
m
dk<- km(design = input, response = output)
dk

testing_input <- as.matrix(seq(0,1,1/99))
m.predict <- predict(m,testing_input)
dk.predict <- predict(dk,testing_input,type='UK')

testing_output=sinewave(testing_input)

plot(testing_input,dk.predict$mean,type='l',col='red',
     ylim=c(min(m.predict$lower95), max(m.predict$upper95)),
     xlab='input',ylab='output')
polygon( c(testing_input,rev(testing_input)),c(m.predict$lower95,
                                               rev(m.predict$upper95)),col = "grey80", border = F)
lines(testing_input,dk.predict$mean,type='l',col='red')
lines(testing_input,m.predict$mean,col='blue')
lines(testing_input,testing_output)
lines(input,output,type='p')




#########13 design points
set.seed(1)
input=as.matrix(seq(0,1,1/12))
output=sinewave(input)

m<-rgasp(design=input, response=output)
m

dk1<- km(design = input, response = output)
dk1

dk2<- km(design = input, response = output)
dk2


testing_input <- as.matrix(seq(0,1,1/99))
m.predict <- predict(m,testing_input)
dk.predict1 <- predict(dk1,testing_input,type='UK')
dk.predict2 <- predict(dk2,testing_input,type='UK')

testing_output=sinewave(testing_input)



####plot
par(mfrow=c(1,2))
plot(testing_input,dk.predict1$mean,type='l',col='red',
     ylim=c(min(m.predict$lower95), max(m.predict$upper95)),
     xlab='input',ylab='output')
polygon( c(testing_input,rev(testing_input)),c(m.predict$lower95,
                                               rev(m.predict$upper95)),col = "grey80", border = F)
lines(testing_input,dk.predict1$mean,type='l',col='red')
lines(testing_input,m.predict$mean,col='blue')
lines(testing_input,testing_output)
lines(input,output,type='p')


plot(testing_input,dk.predict2$mean,type='l',col='red',
     ylim=c(min(m.predict$lower95), max(m.predict$upper95)),
     xlab='input',ylab='output')
polygon( c(testing_input,rev(testing_input)),c(m.predict$lower95,
                                               rev(m.predict$upper95)),col = "grey80", border = F)
lines(testing_input,dk.predict2$mean,type='l',col='red')
lines(testing_input,m.predict$mean,col='blue')
lines(testing_input,testing_output)
lines(input,output,type='p')
dev.off()

########friedman function
set.seed(1)
input <- maximinLHS(n=40, k=5)  # maximin lhd sample
num_obs=dim(input)[1]
output=rep(0,num_obs)
for(i in 1:num_obs){
  output[i]=friedman.5.data(input[i,])
}

m<-rgasp(design=input, response=output)

dk<- km(design = input, response = output)

set.seed(1)
dim_inputs=dim(input)[2]
num_testing_input <- 200    
testing_input <- matrix(runif(num_testing_input*dim_inputs),
                        num_testing_input,dim_inputs)


m.predict <- predict(m,testing_input)
dk.predict <- predict(dk,testing_input,type='UK')


testing_output <- matrix(0,num_testing_input,1)
for(i in 1:num_testing_input){
  testing_output[i]<-friedman.5.data(testing_input[i,])
}

m.rmse=sqrt(mean( (m.predict$mean-testing_output)^2))
m.rmse
dk.rmse=sqrt(mean( (dk.predict$mean-testing_output)^2))
dk.rmse

#######plot
plot(m.predict$mean,testing_output,xlim=c(min(dk.predict$mean,m.predict$mean),
          max(dk.predict$mean,m.predict$mean) ), xlab='prediction',ylab='real output')
lines(dk.predict$mean,testing_output,type='p',col='red')
abline(a=0,b=1)


####emultion with trend function 
# set.seed(1)
# input <- maximinLHS(n=40, k=5)  # maximin lhd sample
# num_obs=dim(input)[1]
# output=rep(0,num_obs)
# for(i in 1:num_obs){
#   output[i]=fried(input[i,])
# }
# set.seed(1)
# dim_inputs=dim(input)[2]
# num_testing_input <- 200    
# testing_input <- matrix(runif(num_testing_input*dim_inputs),
#                         num_testing_input,dim_inputs)
# testing_output <- matrix(0,num_testing_input,1)
# for(i in 1:num_testing_input){
#   testing_output[i]<-fried(testing_input[i,])
# }
# 

colnames(input)=c("x1","x2","x3","x4","x5")
trend.rgasp=cbind(rep(1,num_obs),input)
m.trend<-rgasp(design=input, response=output, trend=trend.rgasp)
dk.trend<- km(formula~x1+x2+x3+x4+x5,design = input, response = output,)
colnames(testing_input)=c("x1","x2","x3","x4","x5")
trend.test.rgasp=cbind(rep(1,num_testing_input),testing_input)
m.trend.predict <- predict(m.trend,testing_input,testing_trend=trend.test.rgasp)
dk.trend.predict <- predict(dk.trend,testing_input,type='UK')



m.trend.rmse=sqrt(mean( (m.trend.predict$mean-testing_output)^2))
m.trend.rmse
dk.trend.rmse=sqrt(mean( (dk.trend.predict$mean-testing_output)^2))
dk.trend.rmse

#######plot trend 
plot(m.trend.predict$mean,testing_output,xlim=c(min(dk.trend.predict$mean,m.trend.predict$mean),
                                          max(dk.trend.predict$mean,m.trend.predict$mean) ), xlab='prediction',ylab='real output')
lines(dk.trend.predict$mean,testing_output,type='p',col='red')
abline(a=0,b=1)


#####see the proportion 
prop.m <- length(which((m.predict$lower95<=testing_output)
                       &(m.predict$upper95>=testing_output)))/num_testing_input
length.m <- sum(m.predict$upper95-m.predict$lower95)/num_testing_input
prop.m
length.m

prop.dk <- length(which((dk.predict$lower95<=testing_output)
                        &(dk.predict$upper95>=testing_output)))/num_testing_input
length.dk <- sum(dk.predict$upper95-dk.predict$lower95)/num_testing_input
prop.dk
length.dk



prop.m.trend <- length(which((m.trend.predict$lower95<=testing_output)
                &(m.trend.predict$upper95>=testing_output)))/num_testing_input
length.m.trend <- sum(m.trend.predict$upper95-
                      m.trend.predict$lower95)/num_testing_input
prop.m.trend
length.m.trend

prop.dk.trend <- length(which((dk.trend.predict$lower95<=testing_output)
                &(dk.trend.predict$upper95>=testing_output)))/num_testing_input
length.dk.trend <- sum(dk.trend.predict$upper95-
                         dk.trend.predict$lower95)/num_testing_input
prop.dk.trend
length.dk.trend

# 


#####real data



###humanity model
data(humanity_model)


###PP GaSP Emulation
m.ppgasp=ppgasp(design=humanity.X,response=humanity.Y,
                nugget.est= TRUE,num_initial_values = 3)
m_pred=predict(m.ppgasp,humanity.Xt)
sqrt(mean((m_pred$mean-humanity.Yt)^2))
sd(humanity.Yt)
mean(m_pred$upper95>humanity.Yt & humanity.Yt>m_pred$lower95)
mean(m_pred$upper95-m_pred$lower95)

##pp gasp with trend
n<-dim(humanity.Y)[1]
n_testing=dim(humanity.Yt)[1]
H=cbind(matrix(1,n,1),humanity.X$foodC)
H_testing=cbind(matrix(1,n_testing,1),humanity.Xt$foodC)
m.ppgasp_trend=ppgasp(design=humanity.X,response=humanity.Y,trend=H,
                      nugget.est= TRUE,num_initial_values = 3)
m_pred_trend=predict(m.ppgasp_trend,humanity.Xt,testing_trend=H_testing)
sqrt(mean((m_pred_trend$mean-humanity.Yt)^2))
mean(m_pred_trend$upper95>humanity.Yt & humanity.Yt>m_pred_trend$lower95)
mean(m_pred_trend$upper95-m_pred_trend$lower95)

sqrt( mean( (mean(humanity.Y)-humanity.Yt)^2 ))


##separate emulator by DiceKriging

k=dim(humanity.Y)[2]
predict_mean_ind_GP_km=matrix(0,n_testing,k)
predict_ind_GP_km_lower_95=matrix(0,n_testing,k)
predict_ind_GP_km_upper_95=matrix(0,n_testing,k)


predict_mean_ind_GP_km_trend=matrix(0,n_testing,k)
predict_ind_GP_km_lower_95_trend=matrix(0,n_testing,k)
predict_ind_GP_km_upper_95_trend=matrix(0,n_testing,k)


for(i_k in 1:k){

  m.km=km(formula=~1,design=humanity.X,response=humanity.Y[,i_k],nugget.estim= TRUE)
  
  m_pred_km=predict(m.km,humanity.Xt,type='UK')
  predict_mean_ind_GP_km[,i_k]=m_pred_km$mean
  predict_ind_GP_km_lower_95[,i_k]=m_pred_km$lower95
  predict_ind_GP_km_upper_95[,i_k]=m_pred_km$upper95
  
  
  m.km_trend=km(formula=~1+foodC,design=humanity.X,response=humanity.Y[,i_k],nugget.estim= TRUE)
  
  m_pred_km_trend=predict(m.km_trend,humanity.Xt,type='UK')
  predict_mean_ind_GP_km_trend[,i_k]=m_pred_km_trend$mean
  predict_ind_GP_km_lower_95_trend[,i_k]=m_pred_km_trend$lower95
  predict_ind_GP_km_upper_95_trend[,i_k]=m_pred_km_trend$upper95
  
}

sqrt(mean( (predict_mean_ind_GP_km-humanity.Yt)^2))
mean(predict_ind_GP_km_upper_95>humanity.Yt & humanity.Yt>predict_ind_GP_km_lower_95)
mean(predict_ind_GP_km_upper_95-predict_ind_GP_km_lower_95)

sqrt(mean( (predict_mean_ind_GP_km_trend-humanity.Yt)^2))
mean(predict_ind_GP_km_upper_95_trend>humanity.Yt & humanity.Yt>predict_ind_GP_km_lower_95_trend)
mean(predict_ind_GP_km_upper_95_trend-predict_ind_GP_km_lower_95_trend)


##TITAN2D data
library(repmis)
source_data("https://github.com/MengyangGu/TITAN2D/blob/master/TITAN2D.rda?raw=True")

input=input_variables[1:50,]
testing_input=input_variables[51:683,]
output=pyroclastic_flow_heights[1:50,which(loc_index[3,]==1)]

testing_output=pyroclastic_flow_heights[51:683,which(loc_index[3,]==1)]
n=dim(output)[1]
n_testing=dim(testing_output)[1]


##delete those location where all output are zero

index_all_zero=NULL
for(i_loc in 1: dim(output)[2]){
  if(sum(output[,i_loc]==0)==50){
    index_all_zero=c(index_all_zero,i_loc)
  }
}

##transforming the output
output_log_1=log(output+1)
k=dim(output_log_1[,-index_all_zero])[2]

system.time(
  for(i in 1:1){
m.ppgasp=ppgasp(design=input[,1:3],response=as.matrix(output_log_1[,-index_all_zero]),trend=cbind(rep(1,n),input[,1]),
                nugget.est=T,max_eval=100,num_initial_values=3)
pred_ppgasp=predict.ppgasp(m.ppgasp,testing_input[,1:3],testing_trend=cbind(rep(1,n_testing),testing_input[,1]))
##transforming back for prediction

m_pred_ppgasp_mean=exp(pred_ppgasp$mean)-1
m_pred_ppgasp_LB=exp(pred_ppgasp$lower95)-1
m_pred_ppgasp_UB=exp(pred_ppgasp$upper95)-1
  }
)
##BV
#user  system elapsed 
#4.346   0.067   4.416 

testing_output_nonallzero=as.matrix(testing_output[,-index_all_zero]) 
sqrt(mean(( (m_pred_ppgasp_mean-testing_output_nonallzero)^2)))
sum(  testing_output_nonallzero>m_pred_ppgasp_LB &  testing_output_nonallzero< m_pred_ppgasp_UB)/length(m_pred_ppgasp_LB)
mean(m_pred_ppgasp_UB-m_pred_ppgasp_LB)


##　

# > sqrt(mean(( (m_pred_ppgasp_mean-testing_output_nonallzero)^2)))
# [1] 0.2999377
# > sum(  testing_output_nonallzero>m_pred_ppgasp_LB &  testing_output_nonallzero< m_pred_ppgasp_UB)/length(m_pred_ppgasp_LB)
# [1] 0.9375439
# > mean(m_pred_ppgasp_UB-m_pred_ppgasp_LB)
# [1] 0.5947367

##try the separate emulator by DK, another package, but it is slow
pred_DK_record=matrix(0,n_testing,k)
lower_95_DK_record=matrix(0,n_testing,k)
upper_95_DK_record=matrix(0,n_testing,k)

set.seed(1)
system.time(
for(i_k in 1:k){
   m.km=km(formula~1+volume, design=input[,1:3],response=output_log_1[,-index_all_zero][,i_k],
           nugget.estim=T)
   m.km_pred=predict(m.km,testing_input[,1:3],type='UK')
   pred_DK_record[,i_k]=exp(m.km_pred$mean)-1
   lower_95_DK_record[,i_k]=exp(m.km_pred$lower95)-1
   upper_95_DK_record[,i_k]=exp(m.km_pred$upper95)-1
}
)
##BV time
#user  system elapsed 
#268.702  25.433 294.430 


sqrt(mean( (pred_DK_record-testing_output_nonallzero)^2))
sum(testing_output_nonallzero>lower_95_DK_record &testing_output_nonallzero< upper_95_DK_record)/length(as.matrix(testing_output)[,-index_all_zero])
mean(upper_95_DK_record-lower_95_DK_record)

#> sqrt(mean( (pred_DK_record-testing_output_nonallzero)^2))
#[1] 0.3016556
#> sum(testing_output_nonallzero>lower_95_DK_record &testing_output_nonallzero< upper_95_DK_record)/length(as.matrix(testing_output)[,-index_all_zero])
#[1] 0.9109993
#> mean(upper_95_DK_record-lower_95_DK_record)
#[1] 0.5295695


###noncrater area
output=pyroclastic_flow_heights[1:50,which(loc_index[1,]==0)]
testing_output=pyroclastic_flow_heights[51:683,which(loc_index[1,]==0)]
n=dim(output)[1]
n_testing=dim(testing_output)[1]


index_all_zero=NULL

##delete those location where all output are zero

for(i_loc in 1: dim(output)[2]){
  if(sum(output[,i_loc]==0)==50){
    index_all_zero=c(index_all_zero,i_loc)
  }
}


output_log_1=log(output+1)
k=dim(output_log_1[,-index_all_zero])[2]

record_MSE=matrix(0,dim(output)[2],3)
record_prop=matrix(0,dim(output)[2],3)

record_length=matrix(0,dim(output)[2],3)

system.time(
  for(i in 1:1){
    m.ppgasp=ppgasp(design=input[,1:3],response=as.matrix(output_log_1[,-index_all_zero]),trend=cbind(rep(1,n),input[,1]),
                    nugget.est=T,num_initial_values=3)
    pred_ppgasp=predict.ppgasp(m.ppgasp,testing_input[,1:3],testing_trend=cbind(rep(1,n_testing),testing_input[,1]))
    
    m_pred_ppgasp_mean=exp(pred_ppgasp$mean)-1
    m_pred_ppgasp_LB=exp(pred_ppgasp$lower95)-1
    m_pred_ppgasp_UB=exp(pred_ppgasp$upper95)-1
  }
)
##non crater
#user  system elapsed 
#19.704   0.534  20.281 


testing_output_nonallzero=as.matrix(testing_output[,-index_all_zero]) 
sqrt(mean(( (m_pred_ppgasp_mean-testing_output_nonallzero)^2)))
sum(  testing_output_nonallzero>m_pred_ppgasp_LB &  testing_output_nonallzero< m_pred_ppgasp_UB)/length(m_pred_ppgasp_LB)
mean(m_pred_ppgasp_UB-m_pred_ppgasp_LB)

# > sqrt(mean(( (m_pred_ppgasp_mean-testing_output_nonallzero)^2)))
# [1] 0.3251579
# > sum(  testing_output_nonallzero>m_pred_ppgasp_LB &  testing_output_nonallzero< m_pred_ppgasp_UB)/length(m_pred_ppgasp_LB)
# [1] 0.9485486
# > mean(m_pred_ppgasp_UB-m_pred_ppgasp_LB)
# [1] 0.6043169

##try the separate emulator by DK
pred_DK_record=matrix(0,n_testing,k)
lower_95_DK_record=matrix(0,n_testing,k)
upper_95_DK_record=matrix(0,n_testing,k)


set.seed(1)
system.time(
  for(i_k in 1:k){
    m.km=km(formula~1+volume, design=input[,1:3],response=output_log_1[,-index_all_zero][,i_k],
            nugget.estim=T)
    m.km_pred=predict(m.km,testing_input[,1:3],type='UK')
    pred_DK_record[,i_k]=exp(m.km_pred$mean)-1
    lower_95_DK_record[,i_k]=exp(m.km_pred$lower95)-1
    upper_95_DK_record[,i_k]=exp(m.km_pred$upper95)-1
  }
)
##noncrater time
#user   system  elapsed 
#1287.875  107.501 1402.044 


sqrt(mean( (pred_DK_record-testing_output_nonallzero)^2))
sum(testing_output_nonallzero>lower_95_DK_record &testing_output_nonallzero< upper_95_DK_record)/length(as.matrix(testing_output)[,-index_all_zero])
mean(upper_95_DK_record-lower_95_DK_record)

#> sqrt(mean( (pred_DK_record-testing_output_nonallzero)^2))
#[1] 0.3337396
#> sum(testing_output_nonallzero>lower_95_DK_record &testing_output_nonallzero< upper_95_DK_record)/length(as.matrix(testing_output)[,-index_all_zero])
#[1] 0.9140655
#> mean(upper_95_DK_record-lower_95_DK_record)
#[1] 0.5345407

