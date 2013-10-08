using namespace std;

#include "controller.h"
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <string.h>


controller::controller(char *data_file, char *grid_file){

  load_grid(grid_file);
  pars.process_data(data_file);
  
  
  p = pars.geno_vec.size();
  s = pars.pheno_vec.size();



  // gsl random number generator
  gsl_rng_env_setup();
     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  long seed = time (NULL) * getpid();    
  gsl_rng_set (r, seed);       

  // default option: running MCMC
  option = 1; 
  mvlr_no_corr = 0;
}


void controller::load_grid(char *grid_file){

  ifstream gfile(grid_file);
  string line;
  istringstream ins;

  while(getline(gfile,line)){

    ins.clear();
    ins.str(line);
    double phi;
    double omg;
    if(ins>>phi>>omg){
      phi2_vec.push_back(pow(phi,2));
      omg2_vec.push_back(pow(omg,2));
    }
  }
  gfile.close();

}



void controller::set_outfile(char *outfile){
  
  if(strlen(outfile)==0)
    outfd = stdout;
  else
    outfd = fopen(outfile,"w");
}



void controller::init(){  
  
  vector<vector<double> > Xc;
  vector<vector<double> > Xg;

   
  for(int i=0;i<p;i++){
    Xg.push_back(pars.geno_vec[i]);
  }

  if(Xg.size() == 0)
    exit(0);
  
  
  // set MVLR parameters and options
  mvlr.set_sigma_option(mvlr_sigma_option);
  if(mvlr_no_corr)
    mvlr.set_no_corr();

  
  mvlr.init(pars.pheno_vec,Xg,Xc);
  


  
  //mvlr.set_prior_option(mvlr_prior_option);
  mvlr.set_effect_vec(phi2_vec,omg2_vec);
  

  vector<int> vec0;
  for(int j=0;j<s;j++)
    vec0.push_back(0);
  
  //vector<int> cfg0;
  int total_cfg = (1<<s);
  
  //start with empty table
  for(int i=0;i<p;i++){
    config.push_back(vec0);
    cfg_map[i] = 0;
  }

}


void controller::prep_single_bf(){

  init();

  int totalc = (1<<s)-1;
 
  for(int i=0;i<p;i++){

    SNP snp(pars.geno_map[i],i);
    vector<vector<int> > mcfg = config;

    for(int j=1;j<=totalc;j++){
      mcfg[i]= get_config(j); 
      fprintf(outfd,"%s_%s  %d  ",snp.name.c_str(), gene.c_str(),j);
      
      for(int k=0;k<phi2_vec.size();k++){
	fprintf(outfd, " %7.4e ", mvlr.compute_log10_ABF(mcfg,phi2_vec[k], omg2_vec[k]));
      }
      fprintf(outfd,"\n");
    }
  }
  
}

// option 2: scan single SNP


void controller::scan(){

  init();
  int totalc = (1<<s)-1;

  for(int i=0;i<p;i++){
  
    SNP snp(pars.geno_map[i],i);
    vector<vector<int> > mcfg = config;
    map<int,int> mcfg_map = cfg_map;
    
    vector<double> wv;
    vector<double> abfv;
    vector<double> config_wv;
    wv.push_back(0);
    
    for(int j=1;j<=totalc;j++){

      mcfg[i]=get_config(j);
      mcfg_map[i] = j;
      double rst = mvlr.compute_log10_ABF(mcfg);
      abfv.push_back(rst);
      config_wv.push_back(pi_vec[j]/(1-pi_vec[0]));
      //wv.push_back(rst+compute_log10_prior(mcfg_map));
      
    }
    
    //double *pp = get_weights(wv);
    
    
    // for(int j=0;j<=totalc;j++){
    //   printf("%10s   %d     %7.4f  %7.5f     ",pars.geno_map[i].c_str(),j, log10(pp[j]),pp[j]);
      
    //  double abf = 0;
    //  if(j!=0)
    //	abf = abfv[j-1];
      
    //  printf("%7.4f\n",abf);
      
    // }
   
    
    fprintf(outfd,"%15s   %9.3f          ",pars.geno_map[i].c_str(), mvlr.log10_weighted_sum(abfv,config_wv) );
    for(int j=1;j<=totalc;j++){
      fprintf(outfd,"(%d) %7.3f   ",j, abfv[j-1]);
    }
    fprintf(outfd,"\n");

    
    
    //delete[] pp;
  }

}
      


// option 3: gene level analysis -- assuming at most single SNP per gene
void controller::gene_level_analysis(){
  
  init();
  int totalc = (1<<s) -1;
  
  vector<double> wts_vec;
  vector<double> abf_vec;
  
  vector<double> configwv;
  for(int i=1;i<=totalc;i++){
    configwv.push_back(pi_vec[i]/(1-pi_vec[0]));
  }
  
  double max_snp_abf = -9999;
  string max_snp;

  for(int i=0;i<p;i++){

    SNP snp(pars.geno_map[i],i);
    vector<vector<int> > mcfg = config;
    map<int,int> mcfg_map = cfg_map;

    vector<double> abfv;
    
    for(int j=1;j<=totalc;j++){
      mcfg[i]=get_config(j);
      mcfg_map[i] = j;
      double rst = mvlr.compute_log10_ABF(mcfg);
      abfv.push_back(rst);
    }
    
    double snp_abf = mvlr.log10_weighted_sum(abfv,configwv);  
    wts_vec.push_back(1.0/p);
    abf_vec.push_back(snp_abf);
    
    if(snp_abf > max_snp_abf){
      max_snp_abf = snp_abf;
      max_snp = snp.name;
    }
      
  }

  double gene_log10_abf = mvlr.log10_weighted_sum(abf_vec,wts_vec);
  double gene_abf = pow(10,gene_log10_abf);
  double gene_pi0 = pow(pi_vec[0],p);
  double gene_post_prob = gene_abf*(1-gene_pi0)/(gene_pi0+  gene_abf*(1-gene_pi0));
  
  fprintf(outfd,"%s\t%d\t%7.3f\t%7.4f\t%s\t%7.3f\n",gene.c_str(),p,gene_log10_abf,gene_post_prob,max_snp.c_str(),max_snp_abf);
  

}




// option 1: run MCMC 



void controller::run_mcmc(int burnin_in, int rep_in){


  burnin = burnin_in;
  total = burnin_in + rep_in;

  init();
  init_mcmc();

  int counter = 1;

  while(counter<= total){
    
    // backup current position
    map<int,int> cfg_map_cp = cfg_map;
    vector<vector<int> > config_cp = config;
    vector<int> select_snp_vec_cp = select_snp_vec;
    
    int p_type = propose();

    double prop_lik = mvlr.compute_log10_ABF(config);
    double prop_prior = compute_log10_prior(cfg_map);
    double prop_pos = prop_lik + prop_prior;
    
    
 


    double f = gsl_rng_uniform(r);
    double ratio = pow(10, prop_pos-curr_pos+log10_q);

    
    //printf("DIAG: %d    %f  %f  %f      %f %f\n",p_type, prop_lik - curr_lik, prop_prior-compute_log10_prior(cfg_map_cp), log10_q, log10(ratio), f);

   if(f<ratio){
      curr_lik = prop_lik;
      curr_pos = prop_pos;
    }else{
      config = config_cp;
      cfg_map = cfg_map_cp;
      select_snp_vec = select_snp_vec_cp;
    }
    
    
    fprintf(stderr,"%6d   %7.3f   ",counter, curr_pos);
    //fprintf(stdout,"%6d   %7.3f   ",counter, curr_pos);
    
    
    for(int i=0;i<p;i++){  
      if(cfg_map[i]>0)
	fprintf(stderr,"[%d:(%d)] ",i,cfg_map[i]);
      //fprintf(stdout,"[%d:(%d)] ",i,cfg_map[i]);
    }
    fprintf(stderr,"\n");              
    //fprintf(stdout,"\n");              
    
    
    if(counter>burnin){
    
      stringstream ss;
      for(int i=0;i<p;i++){
	//config_count[i][cfg_map[i]]++;
	if(cfg_map[i]!=0)
	  ss<<"["<<pars.geno_map[i]<<":"<<cfg_map[i]<<"] ";
      }
      
      string mid = ss.str();
      
      if(mtable.find(mid)==mtable.end()){
	mtable[mid] = 0;
	model_vec.push_back(model(mid,0,curr_pos,curr_lik,cfg_map));
	//printf("new model %s: %f   %f\n",mid.c_str(),curr_lik, mvlr.compute_log10_ABF(config));
      }
      
      mtable[mid]++;
      
      
    }
    
    counter++;
  }
  
  summarize_posterior();

}


double* controller::build_marginal_proposal(map<int,int> &control_map){
  
  MVLR smvlr;
  vector<vector<double> > Xc;
  vector<vector<double> > Xg;


  vector<int> vec0;
  for(int j=0;j<s;j++)
    vec0.push_back(0);
  
  vector<vector<int> > scfg;
  

  for(int i=0;i<p;i++){    
    if(control_map.find(i) != control_map.end())
      Xc.push_back(pars.geno_vec[i]);
    Xg.push_back(pars.geno_vec[i]);
    scfg.push_back(vec0);
  }

  smvlr.set_sigma_option(0);
  smvlr.init(pars.pheno_vec,Xg,Xc);
  smvlr.set_effect_vec(phi2_vec,omg2_vec);

  int totalc = (1<<s)-1;

  int max_index = -1;
  double max_abf = -999;

  vector<double> wv;
  
  for(int i=0;i<p;i++){
    
    vector<vector<int> > mcfg = scfg;
    vector<double> abf_vec;
    vector<double> wts_vec;
    
    for(int j=1;j<=totalc;j++){
      mcfg[i]=get_config(j);
      abf_vec.push_back(smvlr.compute_log10_ABF(mcfg));
      wts_vec.push_back(pi_vec[j]/(1-pi_vec[0]));    
    }
    
    double snp_abf = smvlr.log10_weighted_sum(abf_vec,wts_vec);
    
    if(snp_abf > max_abf){
      max_abf = snp_abf;
      max_index = i;
      
    }
    
    wv.push_back(snp_abf);
    
  }

  if(max_abf <= 1e-8)
    return 0;
  
  control_map[max_index] = 1;
  return get_weights(wv);
  
}



void controller::init_mcmc(){
 
  if(p<20)
    change_prob = 1.0;
  else
    change_prob = 0.85;
  
  remove_prob = 0.0;

  int totalc = (1<<s)-1;
  config_prop_prob = new double*[p];


  double gene_max_abf = -9999;
  int gene_max_snp = -1;
  int gene_max_config = -1;

  for(int i=0;i<p;i++){
    
    SNP snp(pars.geno_map[i],i);
    vector<vector<int> > mcfg = config;
    
    double maxv = -999;
    int maxc = -1;

    vector<double> wv;
    wv.push_back(0); 
 
    for(int j=1;j<=totalc;j++){
      mcfg[i]=get_config(j);
      double rst = mvlr.compute_log10_ABF(mcfg);
      wv.push_back(rst);
      snp.abfv.push_back(rst);
      
      if(rst>maxv){
        maxv = rst;
        maxc = j;
      }
      
    }
    config_prop_prob[i] = get_weights(wv);      
    snp.max_abf = maxv;
    snp.max_config = maxc;

    if(snp.max_abf > gene_max_abf){
      gene_max_abf = snp.max_abf;
      gene_max_snp = i;
      gene_max_config = maxc;
    }

    snp_vec.push_back(snp);
    
  
  }
  
  
  prop_size = 3;
  map<int,int> control_map;
    
  for(int i=0;i<prop_size;i++){
    double *prop_prob = build_marginal_proposal(control_map);
    if(prop_prob != 0)
      prop_prob_vec.push_back(prop_prob);
    else
      break;
  }
  
  vector<double> uniform_wts;
  for(int i=0;i<p;i++){
    uniform_wts.push_back(1);
  }
  prop_prob_vec.push_back(get_weights(uniform_wts));
  
  prop_size = prop_prob_vec.size();

  prop_wts= new double[prop_size];
  double sum = 0;
  for(int i=0;i<prop_size;i++){
    prop_wts[i] = pow(0.5,i+1);
    sum += prop_wts[i];
  }
  
  prop_wts[0] += 1 - sum;
  
  if(prop_wts[prop_size-1]> 0.01){
    prop_wts[0]+= prop_wts[prop_size-1]-0.01;
    prop_wts[prop_size-1] = 0.01;
  }

    
  /*
  prop_prob = new double[p];
  for(int i=0;i<p;i++){
    prop_prob[i] = 0;
    for(int j=0;j<prop_size;j++)
      prop_prob[i] += prop_wts[j]*prop_prob_vec[j][i];
  }
  */
  
  cfg_map[gene_max_snp] = gene_max_config;
  config[gene_max_snp] = get_config(gene_max_config);


  curr_lik = mvlr.compute_log10_ABF(config);
  curr_pos = curr_lik+compute_log10_prior(cfg_map);
  
  fprintf(stderr,"%6d   %7.3f \n",0, curr_pos);
  //fprintf(stdout,"%6d   %7.3f \n",0, curr_pos);
  string mid("");
  mtable[mid] = 1;
  model_vec.push_back(model(mid,0,curr_pos,0,cfg_map));
  
}

  





int controller::propose_snp(){
  
  double draw = gsl_rng_uniform(r);

  if(draw<0.05 && select_snp_vec.size()>0){
    int si = gsl_rng_uniform_int (r, select_snp_vec.size());
    return select_snp_vec[si];
  }

  int index = 0;
  unsigned int *in = new unsigned int[prop_size];
  gsl_ran_multinomial(r,prop_size,1,prop_wts,in);
  for(int j=0;j<prop_size;j++){
    if(in[j]==1){
      index = j;
      break;
    }
  }

  delete[] in;
  
  unsigned int *sn = new unsigned int[p];
  gsl_ran_multinomial(r,p,1,prop_prob_vec[index],sn);
  //gsl_ran_multinomial(r,p,1,prop_prob,sn);
  for(int j=0;j<p;j++){
    if(sn[j]==1){
      delete[] sn;
      return j;
    }
  }  
}

int controller::propose_config(int index){
  int totalc = (1<<s);
  unsigned int *cn = new unsigned int[totalc];
  gsl_ran_multinomial(r,totalc,1,config_prop_prob[index],cn);
  for(int j=0;j<totalc;j++){
    if(cn[j]==1){
      delete[] cn;
      return j;
    }
  }

}


double controller::compute_log10_q(int index, int curr, int target){

  return (log10(config_prop_prob[index][curr])-log10(config_prop_prob[index][target])); 
  
}


int controller::propose(){
  
  
  // use three types of proposals, 
  // 1. change a snp
  // 2. remove a snp in the existing model
  // 3. swap 2 snps
  log10_q = 0;
  double draw = gsl_rng_uniform(r);
  int propose_type = 0;
  if(draw < change_prob){
    propose_type = 1;
  }else if(draw<change_prob + remove_prob){
    propose_type = 2;
  }else{
    propose_type = 3;
  }

  if(propose_type==1){
    int index = propose_snp();      
    // change a SNP
    // pick a new configuration
    int cf_orig = cfg_map[index];
    int cf_new = propose_config(index);
    
    log10_q = compute_log10_q(index,cf_orig, cf_new);
    // bookkeeping
    cfg_map[index] = cf_new;
    config[index] = get_config(cf_new);
    
    if(cf_orig == 0&&cf_new!=0)
      select_snp_vec.push_back(index);
    
    if(cf_new == 0 && cf_orig !=0){
      vector<int> nss_vec;
      for(int i=0;i<select_snp_vec.size();i++){
	if(select_snp_vec[i]!=index){
	  nss_vec.push_back(select_snp_vec[i]);
	}
      }
      select_snp_vec = nss_vec;
    }
  
    return 1;      
  }

  if(propose_type == 2){
    if(select_snp_vec.size()==0)
      return 2;
    // remove a exiting SNP
    int si = gsl_rng_uniform_int (r, select_snp_vec.size());
    vector<int> new_snp_vec;
    for(int i=0;i<select_snp_vec.size();i++){
      if(i==si){
	int index = select_snp_vec[si];
	log10_q = compute_log10_q(index,cfg_map[index],0);
	cfg_map[index] = 0;
	config[index] = get_config(0);
	continue;
      }
      new_snp_vec.push_back(select_snp_vec[i]);
    }
       
    select_snp_vec = new_snp_vec;
    return 2;
  }

  if(propose_type == 3){
    // pick two SNPs and swap 
    
    int index = propose_snp();      
    int target;
    while(1){
      target = propose_snp();
      if(index != target)
	break;
    }
    
    int removal=-1;
    int addition = -1;
    if(cfg_map[index] == 0 &&  cfg_map[target] != 0){
      removal = target;
      addition = index;
    }

    if(cfg_map[index] != 0 &&  cfg_map[target] == 0){
      removal = index;
      addition = target;
    }

    if(removal>=0){
      vector<int> nss_vec;
      for(int i=0;i<select_snp_vec.size();i++){
	if(select_snp_vec[i]!=removal){ 
	  nss_vec.push_back(select_snp_vec[i]);
	}
      }
      nss_vec.push_back(addition);
      select_snp_vec = nss_vec;
    }


    int rcd = cfg_map[target];
    cfg_map[target] = cfg_map[index];
    config[target] = get_config(cfg_map[index]);
    cfg_map[index] = rcd;
    config[index]= get_config(rcd);
    
    log10_q = 0;
      
    return 3;
  }
    
}
  

void controller::set_prior(double pes, double lambda){
  
  int size = 1<<s;
  double pi1 = pes/double(p);
  if(pi1>=1)
    pi1 = 0.5;
  for(int i=0;i<size;i++){
    if(i==0){
      pi_vec.push_back(1-pi1); 
    }else if(i==size-1){
      pi_vec.push_back(pi1*lambda);
    }else{
      pi_vec.push_back(pi1*(1-lambda)/double(size-2));
    }
  }

  
  return;
}

void controller::set_prior(char *prior_file){

  if(prior_file==0)
    return;
  
  int size = 1<<s;

  ifstream ifile(prior_file);
  string line;
  istringstream ins;
  
  vector<double> pvec(size);
  double pi1 = -1;
  while(getline(ifile,line)){
    
    ins.clear();
      ins.str(line);
      double value;
      string type;
      ins>>type;
      
      if(strcmp(type.c_str(),"pi1")==0||strcmp(type.c_str(),"Pi1")==0){
	ins>>value;
        pvec[0] = 1-value;
	pi1 = value;
	continue;
      }
      if(strcmp(type.c_str(),"config")==0){
	for(int i=0;i<size-1;i++){
	  ins>>value;
	  pvec[i+1] = value;
	}
	continue;
      }
  }
  
  if(pi1<=0){
    pi1 = 1.0/p;
    pvec[0] = 1 - pi1;
  }


  for(int i=1;i<size;i++)
    pvec[i] *= pi1;

  pi_vec = pvec;
    
  return;    
  
}

 
double controller::compute_log10_prior(map<int,int> & mcfg_map){
  
  double lp=0;
  for(int i=0;i<p;i++){
    lp += log(pi_vec[mcfg_map[i]]);
  }
  
  return lp/log(10);

}







void controller::summarize_posterior(){

  
  map<string,int> snp_map;
  vector<SNP> post_snp_vec;
 
  int sample_size = total - burnin ;
 
  std::sort(model_vec.begin(), model_vec.end(),sort_model_dec);
  
  for(int i=0;i<model_vec.size();i++){
    string mid = model_vec[i].id;
    model_vec[i].prob = double(mtable[mid])/double(sample_size);
    vector<double> pv = get_rb_incl_prob(model_vec[i].msum);
    
    int counter = 0;
    for(int j=0;j<p;j++){
      if(model_vec[i].msum[j] !=0){
	stringstream ss;
	ss<<pars.geno_map[j]<<":"<<model_vec[i].msum[j];
	string id = ss.str();
	if(snp_map.find(id)==snp_map.end()){
	  snp_map[id]= post_snp_vec.size();
	  SNP psnp(pars.geno_map[j],j);
	  psnp.max_config = model_vec[i].msum[j];
	  psnp.incl_prob = 0;
	  post_snp_vec.push_back(psnp);
	}
	
	post_snp_vec[snp_map[id]].incl_prob += model_vec[i].prob*pv[counter];
	counter++;
      }
    }
    
    if(mid=="")
      mid = string("[NULL]");
    fprintf(outfd, "%5d   %7.4e        %7.3f %7.3f     %s\n",i+1, model_vec[i].prob, model_vec[i].post, model_vec[i].lik, mid.c_str());

  }


  fprintf(outfd,"\n\n");
  
  std::sort(post_snp_vec.begin(),post_snp_vec.end(),sort_snp_dec_by_ip);
  for(int i=0;i<post_snp_vec.size();i++){
    fprintf(outfd,"%5d %10s  %3d  %8.5e          ",i+1,pars.geno_map[post_snp_vec[i].index].c_str(), post_snp_vec[i].max_config, post_snp_vec[i].incl_prob);
    int index = post_snp_vec[i].index;
    for(int j=0; j<snp_vec[index].abfv.size();j++)
      fprintf(outfd,"(%d) %6.3f   ",j+1, snp_vec[index].abfv[j]);    
    fprintf(outfd,"\n");
  }
}




vector<int> controller::get_config(int c){

  vector<int> cfg;

  for(int i=0;i<s;i++){
    int r = c%2;
    cfg.push_back(r);
    c = c/2;
  }

  return cfg;
}



double *controller::get_weights(vector<double>& vec){

  double max = vec[0];
  for(int i=0;i<vec.size();i++){
    if(vec[i]>max)
      max = vec[i];
  }
  double sum = 0;
  for(int i=0;i<vec.size();i++){
    sum += pow(10, (vec[i]-max));
  }
  double *pp = new double[vec.size()];
  for(int i=0;i<vec.size();i++){
    pp[i] = pow(10, (vec[i]-max))/sum;
  }
  return pp;
  
}


vector<double> controller::get_rb_incl_prob(map<int,int>& mcfg_map){
  
  vector<double> rsv;
  vector<int> indv;
  vector<vector<int> > mconfig;

  for(int i=0;i<p;i++){
    mconfig.push_back(get_config(mcfg_map[i]));
    if(mcfg_map[i]!=0){
      indv.push_back(i);
    }
  }
  
  for(int i=0;i<indv.size();i++){
    //printf("%d: %d\n",i,indv[i]);
    rsv.push_back(get_conditional_prob(indv[i],mconfig,mcfg_map));
  }
  return rsv;

}





double controller::get_conditional_prob(int index,vector<vector<int> >& mconfig, map<int,int>& mcfg_map){
  
  vector<vector<int> > wconfig = mconfig;
  map<int,int> wcfg_map = mcfg_map;
  
  int totalc = (1<<s);
  vector<double> postw;
  for(int i=0;i<totalc;i++){
    wconfig[index] = get_config(i);
    wcfg_map[index] = i;
    postw.push_back(mvlr.compute_log10_ABF(wconfig) + compute_log10_prior(wcfg_map));
  }
 
  double *postv = get_weights(postw);
  double rst = postv[mcfg_map[index]];
  delete[] postv;
  return rst;
  
}






bool sort_snp_dec_by_ip(const SNP &lhs, const SNP &rhs){
  return lhs.incl_prob > rhs.incl_prob;
}


bool sort_snp_dec_by_abf(const SNP &lhs, const SNP &rhs){
  return lhs.max_abf > rhs.max_abf;
}

bool sort_model_dec(const model &lhs, const model &rhs){
  return lhs.post > rhs.post;
}

