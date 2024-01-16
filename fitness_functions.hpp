double A_fit(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, int stype){
  double f=0;

  for(int aa=0;aa<Acap;aa++){
    for(int codon=0;codon<Tcap;codon++){
      f += W(stype,aa_to_st[amino_acids(aa)])*eff_mutant_code->coeff(codon,aa)*codon_freq->coeff(stype,codon);
    }
  }
  return f;
}

double Fitness_rate_indep(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, Eigen::MatrixXd* kd){
  double f=1;
  
  for(int stype=0;stype<S;stype++){
    f *= pow(A_fit(codon_freq,eff_mutant_code,stype),l(stype));
  }
  
  return f;
}

double Fitness_rate_dep(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, Eigen::MatrixXd* kd){
  
  double kdbar=0;
  
  //As described in Collins-Hed, Ardell 2019
  for(int i=0; i<Tcap;i++){
    double sum = 0;
    for(int j=0; j<Acap; j++){
      sum += kd->coeff(i,j);
    }
    kdbar+=1/sum;
  }
  kdbar = Tcap/kdbar;

  //////////////////
  
  return Fitness_rate_indep(codon_freq,eff_mutant_code,kd)*(1 - exp(-rexponent*kdbar));  
}
