
#include <Rcpp.h>
#include <math.h>   
using namespace Rcpp;

//function for computing the Wilcoxon rank sum test statistic
double Compute_Wilcoxon_Stat(NumericVector X, IntegerVector Y){
  double stat=0.0;
  for(int i=0;i<X.length();i++){
    if(Y(i) == 1)
      stat += X(i);
  }
  return stat;
}
//function for performing the wilcoxon rank sum test. See entry point from subzero::rcpp_Wilcoxon_PermTest()
// [[Rcpp::export]]
List rcpp_Wilcoxon_PermTest(NumericVector X, IntegerVector Y,IntegerVector B,IntegerVector DoWald,IntegerVector ReportPerm) {
  
  int _perm_vec_length = 1;
  int reportperm = (ReportPerm(0) == 1);
  if(reportperm)
    _perm_vec_length = B(0);
  NumericVector perms(_perm_vec_length);
  
  double stat = Compute_Wilcoxon_Stat(X,Y);
  
  double pval = 0;
  double dist_bigger_or_equal_to = 0;
  int b = 0;
  int dowald = (DoWald(0) == 1);
  bool flag = true;
  IntegerVector sampled_Y;
  double current_stat;
  double Perm_Sum = 0;
  double Perm_Sum_Squared = 0;
  int Perm_Nr_performed = 0;
  //we compute the sum of the vector, and the mean of the distribution
  double _sum = 0.0;
  int _nr_1 = 0;
  double _dist_H0_mean = 0.0;
  for(int i=0;i<X.length();i++){
    _sum += X(i);
    if(Y(i)==1)
      _nr_1++;
  }
  _dist_H0_mean = _sum*((double)_nr_1)/((double)X.length());
  double stat_dist_from_mean = std::abs(stat - _dist_H0_mean);
  double _current_distance = stat_dist_from_mean;
  //NumericVector Perm_Stats(B(0)+2);
  
  if(B(0) == 0){
    flag = false;
    pval = 1.0;
  }
    
  
  while(flag){

    sampled_Y = sample(Y,Y.length(),false);
    current_stat = Compute_Wilcoxon_Stat(X,sampled_Y);
    
    _current_distance = std::abs(current_stat - _dist_H0_mean);

    if(_current_distance >= stat_dist_from_mean)
      dist_bigger_or_equal_to += 1.0;
    pval = (dist_bigger_or_equal_to + 1.0)/ (((double)b) + 2.0); // note that b starts at zero
    
    b++;
    if(b>=B(0))
      flag = false;
    if(dowald && b >= 1000 && pval>0.1)
      flag = false;
    if(dowald && b >= 10000 && pval>0.001)
      flag = false;
    if(dowald && b >= 20000 && pval>0.0005)
      flag = false;
    if(b>=1){
      Perm_Sum += current_stat;
      Perm_Sum_Squared += current_stat*current_stat;
      Perm_Nr_performed ++;
      if(reportperm)
        perms(b-1) = current_stat;
    }
  }
  
  
  List z  = List::create( stat, pval,b,Perm_Sum,Perm_Sum_Squared,Perm_Nr_performed,perms);
  return z ;
}

// function for performing the wilcoxon rank sum test, for a given matrix of permutations
// [[Rcpp::export]]
List rcpp_Wilcoxon_PermTest_Given_Permutations(NumericVector X, IntegerMatrix Y) {
  
  NumericVector stats(Y.ncol());
  for(int i=0;i<Y.ncol();i++){
    stats(i) = Compute_Wilcoxon_Stat(X,Y(_,i));
  }
  List z  = List::create( stats);
  return z ;
}