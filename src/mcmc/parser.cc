using namespace std;

#include "parser.h"
#include "stdlib.h"


void parser::process_data(char *filename){
  
  
  // reading data file
  ifstream infile(filename);
  string line;

  while(getline(infile,line))
    process_line(line);

    
}


void parser::process_line(string line){
    
  // parsing
  char *content = new char[strlen(line.c_str())+1];
  strcpy(content, line.c_str());
 
  char *header = strtok(content, " ");
  if(header==0)
    return;

  char *name = strtok(0, " ");

  

  vector<double> vecv;
  
  while(1){
    char *val = strtok(0, " ");
    if(val==0)
      break;
    vecv.push_back(atof(val));
  }
  

  if(strcmp(header,"pheno") == 0||strcmp(header,"response")==0){
    pheno_vec.push_back(vecv);  
    pheno_map[pheno_vec.size()-1] = string(name);
  }
  if(strcmp(header,"geno") == 0 || strcmp(header,"covariate")==0){
    geno_vec.push_back(vecv);
    geno_map[geno_vec.size()-1] = string(name);
  }

  if(strcmp(header, "controlled") == 0){
    covar_vec.push_back(vecv);
  }
  
  
  delete[] content;


}

  
void parser::output(){
  
  for(int i=0;i<pheno_vec.size();i++)
    printf("%10s  ", pheno_map[i].c_str());
    
  
  for(int i=0;i<geno_vec.size();i++)
    printf("%10s  ", geno_map[i].c_str());

  printf("\n");

  int n = pheno_vec[0].size();
  
  for(int i=0;i<n;i++){
    for(int j=0;j<pheno_vec.size();j++)
      printf("%9.3f  ",pheno_vec[j][i]);
    
    for(int j=0;j<geno_vec.size();j++)
      printf("%.1f  ", geno_vec[j][i]);
    printf("\n");
  }
  
}
