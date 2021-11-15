#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <random>
#include "Eigen/Core"

#include "REML.hpp"
using namespace std;
using namespace Eigen;

// Main funtion used on cov_files.txt with a sequential access
int main(int argc, char *argv[]){
  // Input arguments
  string grmPrefix   = "";
  string mgrmFiles   = "";
  string outREML     = "";
  string phenoFile   = "";
  string covarFile   = "";
  
  bool verbose    =  true;
  bool scaleGRM   = false;
  
  // Indices
  string sw;
  int i,j;
  
  int maxit    = 100;
  float tol    = 1e-6;
  int nrand    = 100;
  int nThreads = 1;
  int mpheno   = 1;
  int nvarComp = 0;
  int nGRMs    = 0;
  int seed     = -1;
  
  VectorXf vars0;
  
  if(argc==1){
    cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
    exit(1);
  }
  
  // Read arguments
  sw = argv[1];
  if (sw == "--help" or sw == "-h"){
    cerr<<"\t--grm/--mgrm : Prefix for output GRM files (as input or output)."<<endl;
    cerr<<"\t--pheno      : Phenotype file just a single column - no missing values. ** IND ordered as in *.fam file. **"<<endl;
    cerr<<"\t--mpheno     : Index of the phenotype to analyse. Default is 1."<<endl;
    cerr<<"\t--covar      : Covariate file."<<endl;
    cerr<<"\t--init       : Initial values (1-component analysis) for 2 variance components.[e.g., --init 0.5 0.5]."<<endl;
    cerr<<"\t--maxit      : Maximum number of AI-REML iterations.[default is "<<maxit<<"]."<<endl;
    cerr<<"\t--nrand      : Number of random vector to simulate. Default is 100."<<endl;
    cerr<<"\t--seed       : Seed for random number generator. Default is -1, i.e. based on current time."<<endl;
    cerr<<"\t--nthreads   : Number of threads used in analysis. Default uses 1."<<endl;
    cerr<<"\t--scale-grm  : Divide all GRMs by their mean trace."<<endl;
    cerr<<"\t--tol        : Convergence tolerance parameter.[default is "<<tol<<"]."<<endl;
    cerr<<"\t--out        : Prefix for output REML results. Default is [none]."<<endl;
    exit(1);
  }else{
    if (argc == 1) {
      cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
      exit(1);
    }
  }
  
  for(i = 1; i<argc;i++){
    sw = argv[i];
    if (sw == "--grm"){
      grmPrefix = argv[i + 1];
    }
    if (sw == "--mgrm"){
      mgrmFiles = argv[i + 1];
    }
    if (sw == "--out"){
      outREML = argv[i + 1];
    }
    if (sw == "--pheno"){
      phenoFile = argv[i + 1];
    }
    if (sw == "--mpheno"){
      mpheno = atoi(argv[i + 1]);
    }
    if (sw == "--seed"){
      seed = atoi(argv[i + 1]);
    }
    if (sw == "--covar"){
      covarFile = argv[i + 1];
    }
    if (sw == "--maxit"){
      maxit = atoi(argv[i + 1]);
    }
    if (sw == "--nrand"){
      nrand = atoi(argv[i + 1]);
    }
    if (sw == "--nthreads"){
      nThreads = atoi(argv[i + 1]);
    }
    if (sw == "--tol"){
      tol = atof(argv[i + 1]);
    }
    if (sw == "--silent"){
      verbose = false;
    }
    if (sw == "--scale-grm"){
      scaleGRM = true;
    }
  }
  
  // Check how many variance components to estimate
  if(phenoFile==""){
    cerr<<"\tA phenotype file must be provided. [Use --pheno][Check help: --help]."<<endl;
    exit(1);
  }
  if(outREML==""){
    cerr<<"\tAn output prefix must be specified. [Use --out][Check help: --help]."<<endl;
    exit(1);
  }
  if(grmPrefix=="" and mgrmFiles==""){
    cerr<<"\tA GRM (or list of GRM) must be provided. [Use --grm or --mgrm][Check help: --help]."<<endl;
    exit(1);
  }else{
    if(grmPrefix!="" and mgrmFiles!=""){
      cerr<<"\tOnly use [--grm] or [--mgrm]. Using both is not allowed. [Check help: --help]."<<endl;
      exit(1);
    }else{
      if(grmPrefix!=""){
        nvarComp = 2;
      }
      if(mgrmFiles!=""){
        ifstream mgrmStream;
        string line;
        mgrmStream.open(mgrmFiles.c_str());
        nvarComp = 1;
        while(mgrmStream){
          getline(mgrmStream,line);
          if(line!="") nvarComp++;
        }
        mgrmStream.close();
        //cout<<"Detected "<<nvarComp<<" variance components.\n";
      }
      if(nvarComp==1){
        cerr<<"\tGRM input file has not been (correctly) specified (#variance component = 1). [Use --grm/--mgrm][Check help: --help]."<<endl;
        exit(1);
      }
      vars0.resize(nvarComp);
      vars0.setZero();
      // Read input values
      for(i = 1; i<argc;i++){
        sw = argv[i];
        if (sw == "--init"){
          for(j=0;j<nvarComp;j++){
            vars0(j) = atof(argv[i + j + 1]);
          }
        }
      }
      //cout<<vars<<endl;
    }
  }
  
  // Checkin threads
  initParallel();
  int nThreadsDetected = nbThreads();
  
  if(nThreads > nThreadsDetected){
    nThreads = nThreadsDetected;
    cerr<<"*** Number of threads constrained to "<<nThreadsDetected<<" available.\n";
  }
  
  setNbThreads(nThreads);
  
  clock_t tic = clock();
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  string logFile = outREML+".greml.log";
  ofstream fileLog(logFile.c_str());
  
  if(verbose){
    cout <<"\033[1;31m>>>> Random Projection based GREML <<<<\033[0m"<<endl;
    //cout <<"\033[1;31m>>>> (c) L. Yengo / License: TBD <<<<\033[0m"<<endl;
    cout <<"#\033[1;34m Analysis starts : \033[0m";
    cout << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday << " at "
         <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
         <<  ".\n";
    if(grmPrefix!="") cout<<"# GRM input: "<<grmPrefix<<".\n";
    if(mgrmFiles!="") cout<<"# GRM input: "<<mgrmFiles<<".\n";
    cout<<"# Phenotype file: "<<phenoFile<<" -- pheno #" <<mpheno<<".\n";
    if(covarFile!="") cout<<"# Covariate file: "<<covarFile<<".\n";
    cout<<"# Output file names: "<<outREML<<".greml\n";
    cout<<"# Number of variance components: "<<nvarComp<<".\n";
    cout<<"# Initial variance components: ";
    for(i=1;i<nvarComp;i++){
      cout<<"| V("<<i<<") = "<<vars0(i-1)<<" ";
    }
    cout<<"| V(e) = "<<vars0(nvarComp-1)<<" |\n";
    cout<<"# Max number of AI-REML iteration(s): "<<maxit<<".\n";
    cout<<"# Number of random vector vectors: "<<nrand<<" (seed = "<<seed<<").\n";
    cout<<"# Convergence tolerance (on parameters value) between AI-REML iterations: "<<tol<<".\n";
    cout<<"# Using "<<nThreads<<"/"<<nThreadsDetected<<" thread(s) detected.\n";
  }

  fileLog <<"\033[1;31m>>>> Random Projection based GREML <<<<\033[0m"<<endl;
  //fileLog <<"\033[1;31m>>>> (c) L. Yengo / License: TBD <<<<\033[0m"<<endl;
  fileLog <<"#\033[1;34m Analysis starts : \033[0m";
  fileLog << (now->tm_year + 1900) << '-'
          << (now->tm_mon + 1) << '-'
          <<  now->tm_mday << " at "
          <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
          <<  ".\n";
  if(grmPrefix!="") fileLog<<"# GRM files: "<<grmPrefix<<".\n";
  if(mgrmFiles!="") fileLog<<"# GRM files: "<<mgrmFiles<<".\n";
  fileLog<<"# Phenotype file: "<<phenoFile<<" -- pheno #" <<mpheno<<".\n";
  if(covarFile!="") fileLog<<"# Covariate file: "<<covarFile<<".\n";
  fileLog<<"# Output file names: "<<outREML<<".greml\n";
  fileLog<<"# Initial variance components: ";
  for(i=1;i<nvarComp;i++){
    fileLog<<"| V("<<i<<") = "<<vars0(i-1)<<" ";
  }
  fileLog<<"| V(e) = "<<vars0(nvarComp-1)<<" |\n";
  fileLog<<"# Max number of AI-REML iteration(s): "<<maxit<<".\n";
  fileLog<<"# Number of random vector vectors: "<<nrand<<" (seed = "<<seed<<").\n";
  fileLog<<"# Convergence tolerance (on parameters value) between AI-REML iterations: "<<tol<<".\n";
  fileLog<<"# Using "<<nThreads<<"/"<<nThreadsDetected<<" thread(s) detected.\n";
  
  nGRMs = nvarComp - 1;
  string grmInput;
  if(nGRMs==1){
    grmInput = grmPrefix;
  }else{
    grmInput = mgrmFiles;
  }
  
  REML myGREML = REML(grmInput,phenoFile,covarFile,
                      mpheno,nGRMs,verbose,fileLog);
  
  //myGREML.summaryGRMs(verbose,fileLog);
  
  if(scaleGRM){
    if(verbose) cout<<"# Scaling GRMs by their mean trace.\n";
    fileLog<<"# Scaling GRMs by their mean trace.\n";
    myGREML.scaleGRM();
  }
  
  myGREML.ai_reml(vars0,maxit,tol,nrand,seed,verbose,fileLog);
  
  myGREML.writeResults(outREML,verbose,fileLog);
  
  time_t t2 = time(0);   // get time now
  struct tm * now2 = localtime( & t2 );
  if(verbose){
    cout <<"#\033[1;34m Analysis ends: \033[0m";
    cout << (now2->tm_year + 1900) << '-'
         << (now2->tm_mon + 1) << '-'
         <<  now2->tm_mday << " at "
         <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
         << ".\n";
  }
  clock_t toc = clock();
  float time_elapsed = (float)(toc - tic) / CLOCKS_PER_SEC;
  
  if(verbose) printf("# Total time elapsed: %f seconds.\n", time_elapsed);
  if(verbose) cout <<"\033[1;31m<<<< Random Projection based GREML >>>\033[0m"<<endl;
  
  fileLog <<"# Analysis ends: ";
  fileLog << (now2->tm_year + 1900) << '-'
          << (now2->tm_mon + 1) << '-'
          <<  now2->tm_mday << " at "
          <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
          << endl;
  fileLog<<"# Time elapsed: "<<time_elapsed<<" seconds\n";
  fileLog <<"\033[1;31m<<<< Random Projection based GREML >>>\033[0m"<<endl;
  fileLog.close();
  return EXIT_SUCCESS;
}


