
#This script is used for producing figure 2 in section 3 (explains why dropping samples based on 
#observed number of counts in subvector is wrong)
#See manuscript for description of graph


pdf(file = paste0('../../Results/NO_FILTER_GRAPH.pdf'),width = 7,height = 2)
par(mar = c(2,2,1,1))
layout(matrix(c(1,1, 2,3), nrow = 1, ncol = 4, byrow = TRUE))
cex_param = 0.8
n_X = n_Y = 16
rate_X = 20
rate_Y = 30
set.seed(1)
pj = sample(x = c(0.9,0.5), size = n_X, replace = T) #sample the overdispersion pattern
qj = sample(x = c(0.9,0.5), size = n_Y, replace = T)
N_X = rpois(n_X,lambda = rate_X + 1*rate_X*(pj>0.5)) # sample subvector depths
N_Y = rpois(n_X,lambda = rate_Y)

X_j = rbinom(n_X,size = N_X,prob = pj)  #sample counts in subvectors
Y_j = rbinom(n_Y,size = N_Y,prob = qj)


X_r = N_X - X_j
Y_r = N_Y - Y_j

AXIS_LIMIT = max(max(c(X_j,Y_j)),max(c(X_r,Y_r)))

#PLOT1 - counts in vectors - prior to subsampling
plot(Y_j,Y_r,pch=1,col = 'black',xlim = c(0,40),ylim = c(0,20),main = 'A',xlab = "Counts in taxon j",ylab = "Counts in reference",asp = 1 ,cex = cex_param)
points(X_j,X_r,pch=4,col = 'black' ,cex = cex_param)
clip(0,40, 0, 20)
abline(0,1,col = 'blue',lwd = 1)
abline(0,0.1,col = 'blue',lwd = 1)
FILTER_FROM = 25
abline(FILTER_FROM,-1,col = 'black',lwd = 1,lty = 2)

#PLOT2 - counts in vector after subsampling, without filtration
lambda_rar = min(N_X,N_Y)
X_tilde_j = rhyper(n_X,m = X_j,n = X_r,k = lambda_rar)
Y_tilde_j = rhyper(n_Y,m = Y_j,n = Y_r,k = lambda_rar)
X_tilde_r = lambda_rar - X_tilde_j
Y_tilde_r = lambda_rar - Y_tilde_j


plot(jitter(Y_tilde_j,0),jitter(Y_tilde_r,6),pch=1,col = 'black',xlim = c(0,lambda_rar),ylim = c(0,lambda_rar),main = 'B',xlab = "Counts in taxon j",ylab = "" ,cex = cex_param)
points(jitter(X_tilde_j,0),jitter(X_tilde_r,6),pch=4,col = 'black' ,cex = cex_param)
clip(0,lambda_rar, 0, lambda_rar)
abline(0,1,col = 'blue',lwd = 1)
abline(0,0.1,col = 'blue',lwd = 1)
arrows(0+0.3,lambda_rar-0.3,x1 = lambda_rar-0.3,y1 = 0.3,col = 'black',lty = 2,lwd = 1,code = 0)


#PLOT3 - counts in vector after subsampling, with filtration

X_j_filtered = X_j[(X_j + X_r)>=FILTER_FROM]
X_r_filtered = X_r[(X_j + X_r)>=FILTER_FROM]
Y_j_filtered = Y_j[(Y_j + Y_r)>=FILTER_FROM]
Y_r_filtered = Y_r[(Y_j + Y_r)>=FILTER_FROM]

lambda_rar_filtered = min(c(X_j_filtered+X_r_filtered,Y_j_filtered+Y_r_filtered))
X_tilde_j_filtered = rhyper(length(X_j_filtered),m = X_j_filtered, n = X_r_filtered, k = lambda_rar_filtered)
Y_tilde_j_filtered = rhyper(length(Y_j_filtered),m = Y_j_filtered, n = Y_r_filtered, k = lambda_rar_filtered)
X_tilde_r_filtered = lambda_rar_filtered - X_tilde_j_filtered
Y_tilde_r_filtered = lambda_rar_filtered - Y_tilde_j_filtered


plot(jitter(Y_tilde_j_filtered,0),jitter(Y_tilde_r_filtered,6),pch=1,col = 'black',xlim = c(0,lambda_rar_filtered),ylim = c(0,lambda_rar_filtered),main = 'C',xlab = "Counts in taxon j",ylab = "" ,cex = cex_param)
points(jitter(X_tilde_j_filtered,0),jitter(X_tilde_r_filtered,6),pch=4,col = 'black' ,cex = cex_param)
clip(0,25, 0, 25)


abline(0,1,col = 'blue',lwd = 1)
abline(0,0.1,col = 'blue',lwd = 1)
arrows(0+0.3,lambda_rar_filtered-0.3,x1 = lambda_rar_filtered-0.3,y1 = 0.3,col = 'black',lty = 2,lwd = 1,code = 0)


dev.off()
par(mfrow=c(1,1))