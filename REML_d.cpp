#include "REML_d.hpp"

void REML_d::readGRM_d(string grmBinfile, MatrixXf &G,int k, bool verbose, ofstream &fileLog){
  ifstream binData(grmBinfile.c_str(), ios::in | ios::binary);
  if(!binData){
    cerr << "\t[readGRM] Error reading file "<<grmBinfile<<endl;
    exit(1);
  }
  if(verbose) cout << "# Reading the GRM from [" + grmBinfile + "]." << endl;
  fileLog << "# Reading the GRM from [" + grmBinfile + "]." << endl;
  
  float f_buf = 0.0;
  int size = sizeof (float);
  
  double sGd  = 0.;
  double sG2d = 0.;
  
  double sGo  = 0.;
  double sG2o = 0.;
  
  double npd = 0.;
  double npo = 0.;
  
  for (int i = 0; i < this->Ngrms(k); i++) {
    for (int j = 0; j <= i; j++) {
      if (!(binData.read((char*) &f_buf, size))){
        throw (("\tError: the size of the [" + grmBinfile + "] file is incomplete?").c_str());
        cerr<<"Not good at all.\n";
      }
      //binData.read((char*) &f_buf, size);
      //GRM(i,j) = GRM(j,i) = f_buf;
      if(k==0){
        int k1 = this->map_i_to_k(i);
        int k2 = this->map_i_to_k(j);
        if(k1>=0 and k2>=0){
          G(k1,k2) = G(k2,k1) = f_buf;
        }
      }else{
        int i_in_G1 = this->mapG_ktoG_1(i,k);
        int j_in_G1 = this->mapG_ktoG_1(j,k);
        if(i_in_G1>=0 and j_in_G1>=0){
          int k1 = this->map_i_to_k(i_in_G1);
          int k2 = this->map_i_to_k(j_in_G1);
          if(k1>=0 and k2>=0){
            G(k1,k2) = G(k2,k1) = f_buf;
          }
        }
      }

      if(i==j){
        sGd  += f_buf;
        sG2d += f_buf * f_buf;
        npd  += 1.;
      }else{
        sGo  += f_buf;
        sG2o += f_buf * f_buf;
        npo  += 1.;
      }
    }
  }
  
  if(npd != Ngrms(k)) cout<<"# **** we might have a problem here!\n";
  
  // cout<<"npd = "<<npd<<" - npo = "<<npo<<".\n";
  
  double meanGd = sGd / npd;
  double meanGo = sGo / npo;
  double varGd  = sG2d / (npd-1.) - meanGd * meanGd * (1. + 1. / (npd - 1.));
  double varGo  = sG2o / (npo-1.) - meanGo * meanGo * (1. + 1. / (npo - 1.));
  
  //varGo = sG2o / npo - meanGo * meanGo;
  
  if(verbose) cout<<"# ---- Mean of diagonal elements: "<<meanGd<<".\n";
  if(verbose) cout<<"# ---- Mean of off-diagonal elements: "<<meanGo<<".\n";
  if(verbose) cout<<"# ---- Variance of diagonal elements: "<<varGd<<".\n";
  if(verbose) cout<<"# ---- Variance of off-diagonal elements: "<<varGo<<".\n";

  fileLog<<"# ---- Mean of diagonal elements: "<<meanGd<<".\n";
  fileLog<<"# ---- Mean of off-diagonal elements: "<<meanGo<<".\n";
  fileLog<<"# ---- Variance of diagonal elements: "<<varGd<<".\n";
  fileLog<<"# ---- Variance of off-diagonal elements: "<<varGo<<".\n";

  binData.close();
};


REML_d::REML_d(string grmInput, string phenoFile, string covarFile,int mpheno, int nGRMs,bool verbose, ofstream &fileLog){
  string line_grm, line_id;
  string token;
  string fid, iid, key;
  string connector = "_";
  string missing_value = "NA";

  // We start
  ifstream tmpStream;
  string grmIdFile;
  
  // Read GRM first
  // Now read single GRM -- should be resize!
  this->nGRMs    = nGRMs;
  this->nvarComp = nGRMs + 1;
  this->vars     = VectorXf::Zero(this->nvarComp);
  this->AI       = MatrixXf::Zero(this->nvarComp,this->nvarComp);
  this->AI_      = MatrixXf::Zero(this->nvarComp,this->nvarComp);
  this->dL       = VectorXf::Zero(this->nvarComp);
  this->dv       = VectorXf::Zero(this->nvarComp);
  this->Ngrms.resize(nGRMs);
  this->grmPrefix = new string[nGRMs];
  
  unordered_map<string, int> *indsInGRMs = new unordered_map<string, int>[nGRMs];
  if(nGRMs==1){
    this->grmPrefix[0] = grmInput;
  }else{
    ifstream mgrmStream;
    mgrmStream.open(grmInput.c_str());
    for(int k=0;k<nGRMs;k++){
      getline(mgrmStream,line_grm);
      if(line_grm!=""){
        stringstream ss_grm;
        ss_grm << line_grm;
        ss_grm >> this->grmPrefix[k];
      }
    }
    mgrmStream.close();
  }
  
  for(int k=0;k<nGRMs;k++){
    grmIdFile = this->grmPrefix[k] + ".grm.id";
    this->Ngrms(k)  = 0;
    tmpStream.open(grmIdFile.c_str());
    while(tmpStream){
      getline(tmpStream,line_id);
      if(line_id!=""){
        stringstream ss_id;
        ss_id << line_id;
        ss_id >> fid;
        ss_id >> iid;
        key = fid + connector + iid;
        indsInGRMs[k].insert({key,this->Ngrms(k)});
        this->Ngrms(k)++;
      }
    }
    tmpStream.close();
    if(verbose){
      cout<<"# Found "<<this->Ngrms(k)<<" Individuals in : "<<this->grmPrefix[k]<<".grm.id."<<endl;
    }
    fileLog<<"# Found "<<this->Ngrms(k)<<" Individuals in : "<<this->grmPrefix[k]<<".grm.id."<<endl;
  }
  

  // nGRMs >= 2 here otherwise boom!
  this->Ngrm_max     = this->Ngrms.maxCoeff();
  this->mapG_1toG_k.resize(this->Ngrms(0),this->nGRMs); //  = MatrixXi::Constant(Ngrms(0),nGRMs,-1);
  this->mapG_ktoG_1.resize(this->Ngrm_max,this->nGRMs); //  = MatrixXi::Constant(Ngrm_max,nGRMs,-1);
  
  // Smart initizialization here...
  
  // First GRM used as reference
  for(int i=0;i<this->Ngrms(0);i++){
    this->mapG_1toG_k(i,0) = i;
    this->mapG_ktoG_1(i,0) = i;
  }
  
  // Loop over GRMs
  for(int k=1;k<this->nGRMs;k++){
    for (auto& x: indsInGRMs[k]){ // Loop over IDs in other GRMs
      key   = x.first;
      int j = x.second; // goes from 0 to Ngrms[k]-1
      unordered_map<string,int>::const_iterator got = indsInGRMs[0].find(key);
      if(got != indsInGRMs[0].end()){
        int i = got->second;
        this->mapG_1toG_k(i,k) = j;
        this->mapG_ktoG_1(j,k) = i;
      }
    }
  }

  // We use the first GRM as reference
  bool *has_pheno = new bool[this->Ngrms(0)]; 
  bool *has_covar = new bool[this->Ngrms(0)];
  
  for(int i=0;i<this->Ngrms(0);i++){
    has_pheno[i] = false;
    has_covar[i] = false;
  }
  
  // Read phenotype
  int Ny = 0;
  string line_y;
  tmpStream.open(phenoFile.c_str());
  while(tmpStream){
    getline(tmpStream,line_y);
    if(line_y!=""){
      Ny++;
    }
  }
  tmpStream.close();
  
  if(verbose){
    cout<<"# Found "<<Ny<<" Individuals in : "<<phenoFile<<endl;
  }
  fileLog<<"# Found "<<Ny<<" Individuals in : "<<phenoFile<<endl;

  tmpStream.open(phenoFile.c_str());
  int nPhenos = -2;
  getline(tmpStream,line_y);
  if(line_y!=""){
    stringstream ss;
    ss << line_y;
    while( ss >> token ){
      ++nPhenos;
    } 
  }
  tmpStream.close();

  if(verbose){
    cout<<"# Found "<<nPhenos<<" phenotype(s) in "<<phenoFile<<".\n";
  }
  fileLog<<"# Found "<<nPhenos<<" phenotype(s) in "<<phenoFile<<".\n";
  if(nPhenos<mpheno){
    cerr<<"\tPhenotype Index [--mpheno] is larger than the number of phenotypes."<<endl;
    exit(1);
  }

  string y_str;
  double y_i;
  VectorXf Y = VectorXf::Zero(this->Ngrms(0));
  tmpStream.open(phenoFile.c_str());

  for(int i=0;i<Ny;i++){
    getline(tmpStream,line_y);
    if(line_y!=""){
      stringstream ss;
      ss << line_y;
      ss >> fid;
      ss >> iid;
      key = fid + connector + iid;
      
      unordered_map<string,int>::const_iterator got = indsInGRMs[0].find(key);
      
      if(got != indsInGRMs[0].end()){
        has_pheno[got->second] = true;
        for(int j=0;j<mpheno;j++){
          ss >> y_str;
        }
        if(y_str==missing_value){
          has_pheno[got->second] = false;
        }else{
          y_i = atof(y_str.c_str());
          Y(got->second) = y_i;
        }
      }
    }
  }
  tmpStream.close();

  // Read covariates
  string line_x;
  MatrixXf X;  //unordered_map<string, VectorXf> _X;
  if(covarFile!=""){
    int Nx = 0;
    tmpStream.open(covarFile.c_str());
    while(tmpStream){
      getline(tmpStream,line_x);
      if(line_x!=""){
        Nx++;
      }
    }
    tmpStream.close();

    if(verbose){
      cout<<"# Found "<<Nx<<" Individuals in : "<<covarFile<<endl;
    }
    fileLog<<"# Found "<<Nx<<" Individuals in : "<<covarFile<<endl;

    tmpStream.open(covarFile.c_str());
    int nCovars = -2;
    getline(tmpStream,line_x);
    if(line_x!=""){
      stringstream ss;
      ss << line_x;
      while( ss >> token ){
        ++nCovars;
      }
    }
    tmpStream.close();

    if(verbose){
      cout<<"# Found "<<nCovars<<" covariate(s) in "<<covarFile<<".\n";
    }
    fileLog<<"# Found "<<nCovars<<" covariate(s) in "<<covarFile<<".\n";

    string X_str;
    int nCovars_plus_1 = nCovars + 1;
    X = MatrixXf::Zero(this->Ngrms(0),nCovars_plus_1);
    for(int i=0;i<this->Ngrms(0);i++){
      X(i,0) = 1.0;
    }
    if(nCovars>=1){
      double X_ij;
      tmpStream.open(covarFile.c_str());
      for(int i=0;i<Nx;i++){
        getline(tmpStream,line_x);
        if(line_x!=""){
          stringstream ss;
          ss << line_x;
          ss >> fid;
          ss >> iid;
          key = fid + connector + iid;
          
          unordered_map<string,int>::const_iterator got = indsInGRMs[0].find(key);
          
          if(got != indsInGRMs[0].end()){
            has_covar[got->second] = true;
            for(int j=0;j<nCovars;j++){
              ss >> X_str;
              if(X_str==missing_value){
                has_covar[got->second] = false;
              }else{
                X_ij = atof(X_str.c_str());
                X(got->second,j+1) = X_ij;
              }
            }
          }
        }
      }
      tmpStream.close();
    }else{ // Mean that covarFile has not any covariate
      cerr<<"\tCovariate file has not covariates in it. Continuing analysis just with the intercept.\n"<<endl;
      cerr<<"\tUse [--covar] to specify a proper file.\n"<<endl;
      for(int i=0;i<this->Ngrms(0);i++){
        has_covar[i] = true; // This is necessary to continue the analysis
      }
    }
  }else{
    X = MatrixXf::Zero(this->Ngrms(0),1);
    for(int i=0;i<this->Ngrms(0);i++){
      has_covar[i] = true;
      X(i,0) = 1.0;
    }
  }
  
  int n = 0;
  this->map_i_to_k.resize(this->Ngrms(0)); // = VectorXi::Constant(Ngrms(0),-1);
  // Smart initialization here
  
  for(int i=0;i<this->Ngrms(0);i++){
    bool have_all_data = (has_covar[i] and has_pheno[i]);
    for(int k=0;k<this->nGRMs;k++){
      have_all_data = have_all_data and (this->mapG_1toG_k(i,k)>=0);
    }
    if(have_all_data){
      this->map_i_to_k(i) = n;
      n++;
    }else{
      this->map_i_to_k(i) = -1;
    }
  }
  if(verbose){
    cout<<"# Found "<<n<<" individuals included in analysis.\n";
  }
  fileLog<<"# Found "<<n<<" individuals included in analysis.\n";

  if(n==0){
    cerr<<"\tStopped because no individuals left to analyse. [Use --pheno][Check help: --help]."<<endl;
    exit(1);
  }

  // Allocate object now
  this->ncovar = X.cols();
  this->y.resize(n);
  this->X.resize(n,this->ncovar);

  if(verbose){
    cout<<"# "<<this->ncovar<<" fixed effect(s) included in the analysis.\n";
  }
  fileLog<<"# "<<this->ncovar<<" fixed effect(s) included in the analysis.\n";

  for(int i=0;i<this->Ngrms(0);i++){
    int k = this->map_i_to_k(i);
    if(k>=0 and k<n){
      this->y(k) = Y(i);
      for(int j=0;j<this->ncovar;j++){
        this->X(k,j) = X(i,j);
      }
      k++;
    }
  }

  // Other stuff
  this->n = n;

  delete [] has_pheno;
  delete [] has_covar;
  //delete [] grmPrefix;
  delete [] indsInGRMs;
};

void REML_d::ai_reml(VectorXf &vars0,int maxit, float tol, int nrand,int seed, bool verbose, ofstream &fileLog){
  string grmBinfile;
  float time_elapsed;
  float traceV_   = 0.;
  float traceP    = 0.;
  float *traceV_G = new float [nGRMs];
  float *tracePG  = new float [nGRMs];
  for(int k=0;k<nGRMs;k++){
    traceV_G[k] = tracePG[k] = 0.;
  }
  LLT<MatrixXf> llt_V; // size n * n
  LLT<MatrixXf> llt_xTV_x; // size ncovar * ncovar
  
  MatrixXf V               = MatrixXf::Zero(n,n);
  MatrixXf G               = MatrixXf::Zero(n,n);
  //MatrixXf iG              = MatrixXf::Zero(Ngrm_max,Ngrm_max);
  
  MatrixXf V_X             = MatrixXf::Zero(n,ncovar);
  MatrixXf xTV_x           = MatrixXf::Zero(ncovar,ncovar);
  MatrixXf xTV_2x          = MatrixXf::Zero(ncovar,ncovar);
  MatrixXf xTV_x_xTV_2x    = MatrixXf::Zero(ncovar,ncovar);

  VectorXf V_y     = VectorXf::Zero(n);
  VectorXf V_Py    = VectorXf::Zero(n);
  VectorXf xTV_y   = VectorXf::Zero(ncovar);
  VectorXf xTV_Py  = VectorXf::Zero(ncovar);
  VectorXf Py      = VectorXf::Zero(n);
  VectorXf PPy     = VectorXf::Zero(n);
  
  MatrixXf *GV_X = new MatrixXf[nGRMs];
  MatrixXf *xTV_GV_x = new MatrixXf[nGRMs];
  MatrixXf *xTV_x_xTV_GV_x = new MatrixXf[nGRMs];
  VectorXf *V_GPy = new VectorXf[nGRMs];
  VectorXf *xTV_GPy = new VectorXf[nGRMs];
  VectorXf *GPy = new VectorXf[nGRMs];
  VectorXf *PGPy = new VectorXf[nGRMs];
  
  for(int k=0;k<nGRMs;k++){
    GV_X[k].resize(n,ncovar);
    xTV_GV_x[k].resize(ncovar,ncovar);
    xTV_x_xTV_GV_x[k].resize(ncovar,ncovar);
    V_GPy[k].resize(n);
    xTV_GPy[k].resize(ncovar);
    GPy[k].resize(n);
    PGPy[k].resize(n);
  }
  
  float *yPGPy = new float[nGRMs];
  
  if(seed==-1){
    srand(time(NULL));
  }else{
    srand(seed); 
  }
  
  MatrixXf x     = MatrixXf::Zero(n,nrand);
  MatrixXf L_x   = MatrixXf::Zero(n,nrand);
  MatrixXf *GL_x = new MatrixXf[nGRMs];
  for(int k=0;k<nGRMs;k++){
    GL_x[k].resize(n,nrand);
  }
  
  // Hutchinson estimator
  clock_t tic_sim_rand = clock();
  for(int i=0;i<n;i++){
    for(int k=0;k<nrand;k++){
      if(rand()%2==0){
        x(i,k) = +1.;
      }else{
        x(i,k) = -1.;
      }
    }
  }
  clock_t toc_sim_rand = clock();
  time_elapsed = (float)(toc_sim_rand - tic_sim_rand) / CLOCKS_PER_SEC;
  if(verbose) cout<<"# Simulating "<<nrand<<" random vectors ["<<time_elapsed<<" sec.].\n";
  fileLog<<"# Simulating "<<nrand<<" random vectors ["<<time_elapsed<<" sec.].\n";
  
  int it = 0;
  float eps = 1.;
  // Initialize
  for(int k=0;k<nvarComp;k++){
    this->vars(k) = vars0(k);
  }

  double logLikOld = 0.;
  double logLikCurrent;
  double deltaLL = 1.;
 
  cout.precision(4); 

  if(verbose){
    cout<<"#\n#\tIter.";
    cout<<"\tlogLik";
    for(int i=1;i<nvarComp;i++){
      cout<<"\tV(G"<<i<<")";
    }
    cout<<"\tV(e)";
    cout<<"\teps";
    //cout<<"\teps(logLik)";
    cout<<"\n";
  }
  
  fileLog<<"#\n#\tIt";
  fileLog<<"\tlogLik";
  for(int i=1;i<nvarComp;i++){
    fileLog<<"\tV("<<i<<")";
  }
  fileLog<<"\tV(e)";
  fileLog<<"\teps";
  //fileLog<<"\teps(logLik)";
  fileLog<<"\n";

  while(it<maxit and eps>tol){
    clock_t tic_reml_iteration = clock();
    
    // Read GRMs here...
    clock_t tic_readGRM = clock();
    grmBinfile   = grmPrefix[0]+".grm.bin";
    readGRM_d(grmBinfile,G,0,false,fileLog);
    clock_t toc_readGRM = clock();
    float time_elapsed = (float)(toc_readGRM - tic_readGRM) / CLOCKS_PER_SEC;
    //if(verbose) cout<<"# Reading GRM #1 took "<<time_elapsed<<" second.\n";
    //fileLog<<"# Reading GRM #1 took "<<time_elapsed<<" second.\n";
    
    V.noalias() = vars(0) * G;
    
    for(int k=1;k<nGRMs;k++){
      clock_t tic_readGRM = clock();
      grmBinfile   = grmPrefix[k]+".grm.bin";
      readGRM_d(grmBinfile,G,k,false,fileLog);
      clock_t toc_readGRM = clock();
      float time_elapsed  = (float)(toc_readGRM - tic_readGRM) / CLOCKS_PER_SEC;
      int indexGRM = k + 1;
      //if(verbose) cout<<"# Reading GRM #"<<indexGRM<<" took "<<time_elapsed<<" second.\n";
      //fileLog<<"# Reading GRM #"<<indexGRM<<" took "<<time_elapsed<<" second.\n";
      
      V.noalias() += vars(k) * G;
    }

    for(int i=0;i<n;i++){
      V(i,i) += vars(nGRMs);
    }

    clock_t tic_cholesky = clock();
    llt_V.compute(V);
    clock_t toc_cholesky = clock();
    float time_elapsed_cholesky  = (float)(toc_cholesky - tic_cholesky) / CLOCKS_PER_SEC;
    //if(verbose) cout<<"# Cholesky Decomposition took "<<time_elapsed_cholesky<<" second.\n";
    //fileLog<<"# Cholesky Decomposition took "<<time_elapsed_cholesky<<" second.\n";

    int success_llt = llt_V.info();
    
    if(success_llt==0){
      // NOW using stochastic approximation
      V_X.noalias()   = llt_V.solve(X);
      V_y.noalias()   = llt_V.solve(y);
      xTV_x.noalias() = X.transpose() * V_X;
      
      llt_xTV_x.compute(xTV_x);
      
      xTV_y.noalias()  = V_X.transpose() * y;
      Py.noalias()     = V_y - V_X * llt_xTV_x.solve(xTV_y);
  
      // Calculate log-likelihood
      double logDetV = 0.;
      for (int i=0;i<n;i++){
        logDetV += 2. * log( llt_V.matrixL()(i,i) );
      }
      double logDetxTV_x = 0.;
      for(int j=0;j<ncovar;j++){
        logDetxTV_x += 2. * log( llt_xTV_x.matrixL()(j,j) );
      }   
      logLikCurrent = -0.5 * ( logDetV + logDetxTV_x + y.dot(Py) );
      deltaLL       = logLikCurrent - logLikOld;
      logLikOld     = logLikCurrent;
    
      // Simulate Xrand matrix -- n x nrand
      L_x.noalias()  = llt_V.matrixL().transpose().solve(x);
      
      for(int k=0;k<nGRMs;k++){
        grmBinfile = grmPrefix[k]+".grm.bin";
        readGRM_d(grmBinfile,G,k,false,fileLog);
        GPy[k].noalias()  = G * Py;
        GV_X[k].noalias() = G * V_X;
        GL_x[k].noalias() = G.transpose() * L_x; //G * L_x;
      }
      
      // Calculate PPy
      V_Py.noalias()   = llt_V.solve(Py);
      xTV_Py.noalias() = V_X.transpose() * Py;
      PPy.noalias()    = V_Py - V_X * llt_xTV_x.solve(xTV_Py);
      
      // Calculate PGPy
      for(int k=0;k<nGRMs;k++){
        V_GPy[k].noalias()   = llt_V.solve(GPy[k]);
        xTV_GPy[k].noalias() = V_X.transpose() * GPy[k];
        PGPy[k].noalias()    = V_GPy[k] - V_X * llt_xTV_x.solve(xTV_GPy[k]);
      }
      
      // Calculate tools for Traces
      xTV_2x.noalias()       = V_X.transpose() * V_X;   // X'V^-2X
      xTV_x_xTV_2x.noalias() = llt_xTV_x.solve(xTV_2x); // (x'V^-1x)^-1 (X'V^-2X)
      
      for(int k=0;k<nGRMs;k++){
        xTV_GV_x[k].noalias() = V_X.transpose() * GV_X[k];
        xTV_x_xTV_GV_x[k].noalias() = llt_xTV_x.solve(xTV_GV_x[k]);
      }
      
      // Trace V_
      traceV_  = (L_x.transpose() * L_x).trace() / nrand;
      //cout<<"Trace(V_) = "<<traceV_<<endl;
      
      for(int k=0;k<nGRMs;k++){
        traceV_G[k] = (GL_x[k].transpose() * L_x).trace() / nrand;
        //cout<<"Trace(V_G["<<k<<"]) = "<<traceV_G[k]<<endl;
      }
      
      // Trace PAi = Trace()
      traceP  =  traceV_  - xTV_x_xTV_2x.trace();
      for(int k=0;k<nGRMs;k++){
        tracePG[k] =  traceV_G[k] - xTV_x_xTV_GV_x[k].trace();
      }
      
      AI(nGRMs,nGRMs) = Py.transpose()  * PPy;
      for(int i=0;i<nGRMs;i++){
        AI(nGRMs,i)  = AI(i,nGRMs) = GPy[i].transpose() * PPy;
        for(int j=0;j<=i;j++){
          AI(i,j) = AI(j,i) = GPy[i].transpose() * PGPy[j];
        }
      }
      AI = 0.5 * AI;
      //cout<<"AI = \n"<<AI<<endl;
      
      for(int k=0;k<nGRMs;k++){
        yPGPy[k] = y.transpose() * PGPy[k];
      }
      double yPPy  = y.transpose() * PPy;
      //cout<<"yPPy = "<<yPPy<<endl;
      
      for(int k=0;k<nGRMs;k++){
        dL(k)   = -0.5 * ( tracePG[k] - yPGPy[k] );
      }
      dL(nGRMs)   = -0.5 * ( traceP  - yPPy );
      //cout<<"dL = \n"<<dL<<endl;

      // Display previous iteration
      clock_t toc_reml_iteration = clock();
      time_elapsed = (float)(toc_reml_iteration - tic_reml_iteration) / CLOCKS_PER_SEC;

      if(verbose){
        cout<<"#\t"<<it;
        cout<<"\t"<<logLikCurrent;
        for(int k=0;k<nGRMs;k++){
          cout<<"\t"<<vars(k);
        }
        cout<<"\t"<<vars(nGRMs);
        if(it==0){
          cout<<"\t---";
        }else{
          cout<<"\t"<<eps;
        }
        //cout<<"\t"<<deltaLL;
        cout<<"\t["<<time_elapsed<<" sec.]\n";
      }
     
      // Log file
      fileLog<<"#\t"<<it;
      fileLog<<"\t"<<logLikCurrent;
      for(int k=0;k<nGRMs;k++){
        fileLog<<"\t"<<vars(k);
      }
      fileLog<<"\t"<<vars(nGRMs);
      if(it==0){
        fileLog<<"\t---";
      }else{
        fileLog<<"\t"<<eps;
      }
      //fileLog<<"\t"<<deltaLL;
      fileLog<<"\t["<<time_elapsed<<" sec.]\n"; 
        
      // Update parameters
      AI_.noalias()   = AI.inverse();
      dv.noalias()    = AI_ * dL;
      vars += dv;
      eps   = dv.transpose() * AI * dv;
      
      if(eps<0.){
        cerr<<"# **** AI matrix is not positive definite!\n";
        cerr<<"# **** Choosing different starting values may help!\n";
        exit(1);
        //break;
      }
      it++;
    }else{
      cerr<<"# **** V matrix is not positive definite!\n";
      cerr<<"# **** Choosing different starting values may help!\n";
      exit(1);
      //break;
    }// success LLT
  }
  
  delete [] traceV_G;
  delete [] tracePG;
  delete [] GV_X;
  delete [] xTV_GV_x;
  delete [] xTV_x_xTV_GV_x;
  delete [] V_GPy;
  delete [] xTV_GPy;
  delete [] GPy;
  delete [] PGPy;
  delete [] yPGPy;
  delete [] GL_x;
  
}

void REML_d::writeResults(string outREML, bool verbose, ofstream &fileLog){
  string fileREML = outREML+".greml.res";
  ofstream fileREMLoutput(fileREML.c_str());
  
  // Header
  fileREMLoutput<<"Parameter";
  for(int i=1;i<=this->nGRMs;i++){
    fileREMLoutput<<"\tVg("<<i<<")";
  }
  fileREMLoutput<<"\tVe\n";
  
  // Variance components
  fileREMLoutput<<"Estimates";
  for(int i=0;i<this->nGRMs;i++){
    fileREMLoutput<<"\t"<<this->vars(i);
  }
  fileREMLoutput<<"\t"<<this->vars(this->nGRMs)<<"\n";
  
  // Standard Errors
  fileREMLoutput<<"Stand.Err";
  for(int i=0;i<this->nGRMs;i++){
    fileREMLoutput<<"\t"<<sqrt(this->AI_(i,i));
  }
  fileREMLoutput<<"\t"<<sqrt(this->AI_(this->nGRMs,this->nGRMs))<<"\n";
  
  // Variance - Covariance Matrix
  for(int i=0;i<this->nvarComp;i++){
    if(i<this->nGRMs){
      fileREMLoutput<<"CV"<<(i+1)<<"_";
    }else{
      fileREMLoutput<<"CVe_";
    }
    for(int j=0;j<this->nvarComp;j++){
      fileREMLoutput<<"\t"<<this->AI_(i,j);
    }
    fileREMLoutput<<"\n";
  }
  fileREMLoutput.close();
  if(verbose) cout<<"# Results written in "<<fileREML<<".\n";
  fileLog<<"# Results written in "<<fileREML<<".\n";
}
