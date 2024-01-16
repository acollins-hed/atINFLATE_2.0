struct aaRSs
{
  Eigen::VectorXf aas;
  //The iis of the aaRSs are 2xA matrices.
  Eigen::MatrixXi iis;
  aaRSs(){
    Eigen::MatrixXi a(2,A);
    Eigen::VectorXf aa(A);
    a.setZero();
    iis = a;
    for(int i=0;i<A;i++)
      aa(i) = ((float)i)/(A-1);
    aas=aa;
  }
  aaRSs(Eigen::MatrixXi x, Eigen::VectorXf y){
    iis = x;
    aas = y;
  }
  void  print(){
    for(int i = 0;i<2;i++){
      for(int j=0;j<A;j++){
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

  string print_aars(int i, int j){
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

