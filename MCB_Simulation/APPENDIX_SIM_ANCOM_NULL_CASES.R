library(ancom.R)
library(doRNG)
library(doParallel)


cases_mat = expand.grid(lambda = c(30,60),m = c(50,100),sample_size = c(50,100))
B = 200


compute_FWER_for_case = function(case_id,B){
  library(ancom.R)
  rejections = c(0,0,0)  
  n_X = n_Y = cases_mat$sample_size[case_id]
  m = cases_mat$m[case_id]
  lambda = cases_mat$lambda[case_id]
  for(b in 1:B){
    sim_counts = matrix(data = rpois((n_X+n_Y)*m,30),ncol = m)
    
    ANCOM_MAT = cbind(sim_counts,c(rep(0,n_X),rep(1,n_Y)))
    ANCOM_MAT = data.frame(ANCOM_MAT)
    for(multcorr_ind in 1:3){
      ancom_res  = ANCOM(OTUdat = ANCOM_MAT , multcorr = multcorr_ind)
      if(ancom_res$detected!= "No significant OTUs detected"){
        rejections[multcorr_ind] = rejections[multcorr_ind]+1
      }  
    }  
  }
  FWER = rejections/B
  return(FWER)

}
print('Start Time')
print(Sys.time())
cl <- makeCluster(8)
registerDoParallel(cl)
res <- foreach(i=1:nrow(cases_mat), .options.RNG=1,.combine = rbind) %dorng% { compute_FWER_for_case(i,B) }
stopCluster(cl)
print('Stop Time')
print(Sys.time())


res_mat = cbind(cases_mat,res)
colnames(res_mat)[ncol(cases_mat)+(1:3)] = paste0('multcorr = ',1:3)
for(i in 1:3){
  res_mat[,ncol(cases_mat)+i] = round(res_mat[,ncol(cases_mat)+i],2)  
}


save(res_mat,file = '../../Results/ANCOM_T1E_SIM.Rdata')

library(xtable)
print(xtable(res_mat), include.rownames=FALSE)
