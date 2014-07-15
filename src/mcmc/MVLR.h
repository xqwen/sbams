#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <vector>



class MVLR {
  
 private:
  
  int n; // sample size
  int s; // subgroup size
  int q; // control size
  int p; // genotype size

  int m; // wishart prior default m = s-1

  int ep; //effective genotype size

  gsl_matrix *Y; // phenotype matrix nxs
  gsl_matrix *Xg; // genotype matrix nxp
  gsl_matrix *Xc; // control matrix nxq
  gsl_matrix *H; // wishart prior sxs 


  // used throughout
  gsl_matrix *T;
  gsl_matrix *Sigma0;
  gsl_matrix *Sigma0_inv;
  
  // used per configuration
  gsl_matrix *eVb; 
  gsl_matrix *eVg_inv;
  gsl_matrix *Gamma;

  // used per (phi, omg) value
  gsl_matrix *Wg; // effect prior  

  // maybe used in either case (depends on option)
  gsl_matrix *Sigma; // residual covariance estimate sxs
  gsl_matrix *Sigma_inv; // Sigma^{-1}
 
  
  vector<double> omg2_vec; //effect size grid
  vector<double> phi2_vec; //effect size grid


 private:
  // options
  double sigma_option;  // 0 to 1, mixture fraction of mle of Sigma under the alternative model, default 0
  int prior_option;  // 1-- meta prior      2-- diagonal prior
  int no_corr;   // ignore the correlation between subgroups is set to 1, default 0
  int no_incpt;  // regression without intercept term
  
 private:
  
  void compute_common();
  
  // utilites for computing residual error cov
  void compute_Sigma(vector<vector<int> >& config);
  void compute_Sigma_null();
  void compute_Sigma_mle(vector<vector<int> >& config);
  void invert_Sigma();
  gsl_matrix *compute_residual(gsl_matrix *y, gsl_matrix *X, int size, double &factor);

  
  // utilites for configuration specific computation
  
  void construct_Gamma(vector<vector<int> >& config, vector<int> &noz_vec);
  void construct_meta_Gamma(vector<vector<int> >& config, vector<int> &noz_vec);
  void construct_diag_Gamma(vector<vector<int> >& config, vector<int> &noz_vec);
  
  void set_Wg(double phi2, double omg2);
  
  // compute stats common for a configuration
  void compute_stats(vector<int> &noz_vec);
  
  // evaluating ABF
  double compute_log10_ABF(gsl_matrix *Wg);
  
  gsl_matrix *vec(gsl_matrix *M, int a, int b);
  gsl_matrix *kron (gsl_matrix *M, gsl_matrix *L, int a, int b);
  gsl_matrix *kron2 (gsl_matrix *M, int mr, int mc, gsl_matrix *L, int lr, int lc);

  void print_matrix(gsl_matrix *M, int a, int b);
  
  
 public:
  
  // interface
  // empty constructor, assign default options
  MVLR(){
    sigma_option = 0.0; 
    prior_option=1;
    no_corr = 0;
    no_incpt = 0;
  }


  // init
  void init(vector<vector<double> > & Y_in, vector<vector<double> > & Xg_in, vector<vector<double> > & Xc_in);
  
  // options
  void set_IW_prior(gsl_matrix *H_in, int m_in);
  void set_effect_vec(vector<double> &phi2,vector<double>& omg2_vec);
  void set_sigma_option(double option){
    sigma_option = option;
  }
  void set_prior_option(int option){
    prior_option = option;
  }

  void set_no_incpt(){
    no_incpt = 1;
  }

  void set_no_corr(){
    no_corr=1;
  }
  
  double log10_weighted_sum(vector<double> &vec, vector<double> &wts);

  double compute_log10_ABF(vector<vector<int> > &indicator);
  
  double compute_log10_ABF(vector<vector<int> >& indicator, double phi2, double omg2);
 
  ~MVLR();
  


};


