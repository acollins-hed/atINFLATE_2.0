struct Current_Genotype
{
 public:
  tRNAs trnas;
  aaRSs aarss;
  Eigen::MatrixXd code;
  Eigen::MatrixXd kd;
  void Get_Code(){
    kd=Get_kd();
    Eigen::MatrixXd c(T,A);
    double hx;
    for(int i=0;i<Tcap;i++){
      hx = inverted_hsum(kd,i);
      //hx = hmean(kd,i)/double(Tcap);
      for(int j=0;j<Acap;j++){
	c(i,j) = hx/kd(i,j);
      }
    }
    code = c;
  }
  
  Current_Genotype(){
    tRNAs t;
    aaRSs a;
    trnas=t;
    aarss=a;
    Get_Code();
  }
  Current_Genotype(tRNAs t, aaRSs a)
  :trnas(t),aarss(a){
    Get_Code();
  }
  Eigen::MatrixXd Get_kd(){
    Eigen::MatrixXd krate(T,A);
    krate.setZero();

    for(int i=0;i<Tcap;i++){
      for(int j=0;j<Acap;j++){
	krate(i,j) = kdnc*exp(-1*epsilon*(__builtin_popcount((trnas.iis(1,i)&aarss.iis(1,j))&(~(trnas.iis(0,i)^aarss.iis(0,j))))));
      }
    }    
    
    
    return krate;
  }
};
