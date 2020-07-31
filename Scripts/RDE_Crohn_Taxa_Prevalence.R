ind_filter = as.numeric(which(apply(X>0,2,sum)>=2))

prev_H = as.numeric(apply(X[Y==0,ind_filter]>0,2,mean))
prev_S = as.numeric(apply(X[Y==1,ind_filter]>0,2,mean))

summary(prev_H)
summary(prev_S)

boxplot(prev_H,prev_S)
