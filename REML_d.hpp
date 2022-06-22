#ifndef REML_d_H
#define REML_d_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <math.h>
#include <unordered_map>
#include "Eigen/Core"
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;


class REML_d{
public:
  // Data
  VectorXf y;  // vector of phenotype
  MatrixXf X;  // matrix of covariates
  
  // Maps for GRMs
  int Ngrm_max;
  VectorXi Ngrms;
  MatrixXi mapG_1toG_k;
  MatrixXi mapG_ktoG_1;
  VectorXi map_i_to_k;
  string *grmPrefix;
  
  //int maxit;
  //float tol;
  //int nrand;
  int n;
  int nvarComp;
  int nGRMs;
  int ncovar;
  
  // Parameters to estimate
  VectorXf vars;
  MatrixXf AI;
  MatrixXf AI_;
  VectorXf dL;
  VectorXf dv;
  
  // Multiple GRMs -- specifies nvarComp
  REML_d(string grmInput, string phenoFile, string covarFile, 
       int mpheno, int nGRMs,bool verbose, ofstream &fileLog);
  ~REML_d(){};
  
  void readGRM_d(string grmBinfile, MatrixXf &G,int k, bool verbose, ofstream &fileLog);
  
  void ai_reml(VectorXf &vars0,int maxit, float tol, int nrand,
               int seed, bool verbose, ofstream &fileLog);

  void writeResults(string outREML, bool verbose, ofstream &fileLog);
};
#endif
