double inverted_hsum(Eigen::MatrixXd k, int i){
  double sum=0;
  for(int j=0;j<Acap;j++){
	sum += 1/k.coeff(i,j);
  }
  return (double(1)/(sum));
}

double Transition_Probability(double current_fitness, double mutant_fitness,int hamming){
  double Fixation_Probability;
  if(current_fitness == mutant_fitness)
    Fixation_Probability = 1/((double) N);
  else
    Fixation_Probability = (1 - (current_fitness)/(mutant_fitness))/(1 - pow((current_fitness)/(mutant_fitness),N));
  return (N*pow(mu,hamming)*Fixation_Probability); 
}
