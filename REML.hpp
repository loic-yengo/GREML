#ifndef REML_H
#define REML_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <unordered_map>
#include "Eigen/Core"
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;


class REML{
public:
  // Data
  VectorXf y;  // vector of phenotype
  MatrixXf X;  // matrix of covariantes
  MatrixXf *G; // list of GRMs
  
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
  REML(string grmInput, string phenoFile, string covarFile, 
       int mpheno, int nGRMs,bool verbose, ofstream &fileLog);
  ~REML(){};
  
  void summaryGRMs(bool verbose, ofstream &fileLog);
  
  void ai_reml(VectorXf &vars0,int maxit, float tol, int nrand,
               int seed, bool verbose, ofstream &fileLog);
  
  void scaleGRM(); // Divide all GRM by their mean trace
  
  void writeResults(string outREML, bool verbose, ofstream &fileLog);
};
#endif
