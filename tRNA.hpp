struct tRNAs
{
  Eigen::VectorXi anticodons;
  //The iis of the tRNAs are 2xT matrices.
  Eigen::MatrixXi iis;
  tRNAs(){
    Eigen::MatrixXi t(2,T);
    Eigen::VectorXi a(T);
    t.setZero();
    iis = t;
    for(int i=0;i<T;i++)
      a(i) = i;
    anticodons=a;
  }
  tRNAs(Eigen::MatrixXi x, Eigen::VectorXi z){
    iis = x;
    anticodons = z;
  }
  void  print(){
    for(int i = 0;i<2;i++){
      for(int j=0;j<T;j++){
	for(int ii=n-1;ii>=0;ii--){
	  if(iis(i,j)&(1<<ii))
	    cout<<"1";
	  else
	    cout<<"0";
	}
	cout<<" ";
      }
      cout<<endl;
    }
  }

  string print_trna(int i, int j){
    string s = "";
    for(int ii=n-1;ii>=0;ii--){
      if(iis(i,j)&(1<<ii))
	s += "1";
      else
	s += "0";
    }
    return s;
  }
  
};
