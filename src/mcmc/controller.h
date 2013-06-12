#include "parser.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include "MVLR.h"
#include <set>
#include <map>


class SNP {

 public:

  string name;
  int index;
  double max_abf;
  double incl_prob;
  int max_config;
  
  vector<double> abfv;

  SNP(string name_in, int index_in){
    name = name_in;
    index = index_in;
  }


};

bool sort_snp_dec_by_ip(const SNP &lhs, const SNP &rhs);
bool sort_snp_dec_by_abf(const SNP &lhs, const SNP &rhs);

class model {
 public:

  string id;
  double prob;

  map<int,int> msum;
  double post;
  double lik;

  model(string id_, double prob_, double post_,double lik_, map<int,int> msum_){
    id= id_;
    prob = prob_;
    post = post_;
    lik = lik_;
    msum = msum_;
  }
};


bool sort_model_dec(const model &lhs, const model &rhs);


class controller {
  
 private:

  int option; // marginal scan or MCMC


  parser pars;
  string gene;

  vector<double> phi2_vec;
  vector<double> omg2_vec;

  void load_grid(char *grid_file);
  
  MVLR mvlr;

  int p;
  int s;

  vector<vector<int> > config;
  map<int,int> cfg_map;
  vector<int> select_snp_vec;

  
  vector<double> pi_vec; //prior for configs


  vector<double> abf_vec;
  int total;
  int burnin;

  // proposal frequency
  double remove_prob;
  double change_prob;

  int prop_size;
  double *prop_wts;
  double *prop_prob;
  vector<double*> prop_prob_vec;
  
  double **config_prop_prob;
  
      
  double mvlr_sigma_option;
  int mvlr_no_corr;
  int mvlr_prior_option;
  

  const gsl_rng_type * T;
  gsl_rng * r;

  double curr_lik;
  double curr_pos;
  double log10_q;


  map <string,int> mtable;
  vector<model> model_vec;
  vector<SNP> snp_vec;
  
  

 public:
  
  controller(char *data_file, char *grid_file);
    
  void set_prior(char *prior_file);
  void set_prior(double pes=1, double lambda=0.5);


  
  void run_mcmc(int burnin_in, int rep_in);
  void prep_single_bf();
  void scan();
  void gene_level_analysis();

  
  void set_gene(string gname){
    gene=gname;
  }
  
  void init_mcmc();
  
  
  
  void set_mvlr_option(double sigma_option, int no_corr){
    mvlr_sigma_option= sigma_option;
    mvlr_no_corr = no_corr; 
  }




  void summarize_posterior();
  
 private:

  //MCMC sampler
  int propose();
  int propose_snp();
  int propose_config(int index);
  //int propose_removal(double *rp, int size);
  
  double* build_marginal_proposal(map<int,int> &control_map);
  
  double compute_log10_prior(map<int,int>& mcfg_map);
  double compute_log10_q(int index, int curr, int target);
  void init();
  vector<int> get_config(int c);
  double *get_weights(vector<double>& vec);


  double get_conditional_prob(int index, vector<vector<int> >& mconfig, map<int,int>& mcfg_map);  
  vector<double> get_rb_incl_prob(map<int,int>& mcfg_map);


};  


