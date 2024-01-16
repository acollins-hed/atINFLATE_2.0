struct Population
{
  Current_Genotype genotype;
  double fitness;
  Eigen::MatrixXd codon_frequency;
  Population(){
    Current_Genotype genotype;
    Get_Codon_Freq();
    fitness = fitness_func(&codon_frequency, &genotype.code, &genotype.kd);
  }
  Population(Current_Genotype g)
    :genotype(g){
    Get_Codon_Freq();
    fitness = fitness_func(&codon_frequency, &genotype.code, &genotype.kd);
  }
  void Get_Codon_Freq(){
    Eigen::MatrixXd Q(T,T),wfit(T,T),diag(S,T);
    diag = Diagonal();
    codon_frequency.resize(S,T);
    for(int i=0;i<S;i++){
      wfit.setZero();
      wfit.diagonal() = diag.row(i).array();
      Q = mutation*wfit;
      Eigen::EigenSolver<Eigen::MatrixXd> ev(Q);
      //I have seen claims on the internet that Eigen automatically puts the
      //largest eigenvalue into cell 0 (I've not seen this on Eigen's own
      //website) and I've found this to be untrue, hence the following line.
      codon_frequency.row(i) = ev.eigenvectors().col(max_eigenvalue(&ev.eigenvalues())).real();
      //The following makes sure the components sum to one. The eigenvectors are
      //already normalized but of course that just makes their norm 1, not their
      //components sum to one.
      Sum_to_one(&codon_frequency);
    }
  }

private:
  void Sum_to_one(Eigen::MatrixXd * emat){
    double sum=0;
    for(int i=0;i<emat->rows();i++){
      for(int j=0;j<emat->cols();j++){
	sum = emat->coeff(i,j) + sum;
      }
      emat->row(i) *= (1/sum);
	sum = 0;
    }
  }

  Eigen::MatrixXd Diagonal(){
    Eigen::MatrixXd diag(S,T);
    for(int stype = 0;stype<S;stype++){
      for(int codon=0;codon<Tcap;codon++){
	diag(stype,codon) = 0;
	for(int alpha=0;alpha<Acap;alpha++){
	  diag(stype,codon) += W(stype,aa_to_st[amino_acids(alpha)])*genotype.code(codon,alpha);
	}
      }
    }
    return diag;
  }

  int max_eigenvalue(const Eigen::EigenSolver<Eigen::MatrixXcd>::EigenvalueType * evals){
    int max = 0,i;
    for(i=1;i<evals->size();i++){
      if(evals->coeff(max).real() < evals->coeff(i).real()){
	max = i;
      }
    }
    return max;
  }
  
};
