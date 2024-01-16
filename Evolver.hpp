struct Evolver
{
  Eigen::MatrixXd single_mutant, double_mutant_wout, double_mutant_win;

  Population population;
  double current_fitness;


  Evolver(Population p)
    :population(p){
    current_fitness = population.fitness;
    single_mutant = Eigen::MatrixXd(L,5);
    double_mutant_win = Eigen::MatrixXd(n_d_win,5);
    double_mutant_wout = Eigen::MatrixXd(n_d_wout,3);
    Get_Mutant_Web();
  }

  void Get_Mutant_Web(){
    double sum = 0;
    int index = 0;
    single_mutant.setZero();
    double_mutant_wout.setZero();
    double_mutant_win.setZero();
    
    //Beginning of single_mutant loops
    //
    //column 0 = 0 for tRNA, 1 for aaRS
    //column 1 = row of mutation
    //column 2 = column of mutation
    //column 3 = value in (row,column)
    //column 4 = transition probability

    for(int row=0;row<2;row++){
      for(int col=0;col<Tcap;col++){
	for(int i=0;i<n;i++){
	  Population mutant_pop = population;
	  mutant_pop.genotype.trnas.iis(row,col) ^= 1<<i;
	  mutant_pop.genotype.Get_Code();
	  mutant_pop.fitness = fitness_func(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);

	  single_mutant(index,0) = 0;
	  single_mutant(index,1) = row;
	  single_mutant(index,2) = col;
	  single_mutant(index,3) = mutant_pop.genotype.trnas.iis(row,col);
	  single_mutant(index,4) = Transition_Probability(current_fitness,mutant_pop.fitness,1);
	  sum += single_mutant(index,4);
	  index++;
	}
      }
    }

    for(int row=0;row<2;row++){
      for(int col=0;col<Acap;col++){
	for(int i=0;i<n;i++){
	  Population mutant_pop = population;
	  mutant_pop.genotype.aarss.iis(row,col) ^= 1<<i;
	  mutant_pop.genotype.Get_Code();
	  mutant_pop.fitness = fitness_func(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);

	  single_mutant(index,0) = 1;
	  single_mutant(index,1) = row;
	  single_mutant(index,2) = col;
	  single_mutant(index,3) = mutant_pop.genotype.aarss.iis(row,col);
	  single_mutant(index,4) = Transition_Probability(current_fitness,mutant_pop.fitness,1);
	  sum += single_mutant(index,4);
	  index++;
	}
      }
    }
    

    //Beginning of double_mutant_win loops
    //
    //column 0 = 0 for tRNA, 1 for aaRS
    //column 1 = row of mutation
    //column 2 = column of mutation
    //column 3 = value in (row,column)
    //column 4 = transition probability
    index = 0;
    for(int row=0;row<2;row++){
      for(int col=0;col<Tcap;col++){
	for(int i=0;i<n-1;i++){
	  for(int j=n-1;j>i;j--){
	    Population mutant_pop = population;
	    mutant_pop.genotype.trnas.iis(row,col) ^= 1<<i;
	    mutant_pop.genotype.trnas.iis(row,col) ^= 1<<j;
	    mutant_pop.genotype.Get_Code();
	    mutant_pop.fitness = fitness_func(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);

	    double_mutant_win(index,0) = 0;
	    double_mutant_win(index,1) = row;
	    double_mutant_win(index,2) = col;
	    double_mutant_win(index,3) = mutant_pop.genotype.trnas.iis(row,col);
	    double_mutant_win(index,4) = Transition_Probability(current_fitness,mutant_pop.fitness,2);
	    sum += double_mutant_win(index,4);
	    index++;
	  }
	}
      }
    }

    for(int row=0;row<2;row++){
      for(int col=0;col<Acap;col++){
	for(int i=0;i<n-1;i++){
	  for(int j=n-1;j>i;j--){
	    Population mutant_pop = population;
	    mutant_pop.genotype.aarss.iis(row,col) ^= 1<<i;
	    mutant_pop.genotype.aarss.iis(row,col) ^= 1<<j;
	    mutant_pop.genotype.Get_Code();
	    mutant_pop.fitness = fitness_func(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);

	    double_mutant_win(index,0) = 1;
	    double_mutant_win(index,1) = row;
	    double_mutant_win(index,2) = col;
	    double_mutant_win(index,3) = mutant_pop.genotype.aarss.iis(row,col);
	    double_mutant_win(index,4) = Transition_Probability(current_fitness,mutant_pop.fitness,2);
	    sum += double_mutant_win(index,4);
	    index++;
	  }
	}
      }
    }


    //Beginning of double_mutant_wout loops
    //
    //column 0 = row of first mutation in single_mutant
    //column 1 = row of second mutation in single_mutant
    //column 2 = transition probability

    //An explanation for the indexing of the following set of nested
    //loops. i and j are used to pull mutation information from
    //single_mutant. i begins at 0 and ticks up to L-n-1. j begins
    //at q*n+n and ticks up to L-1 where  q is the quotient i with n, i.e.
    //i = qn + r. This way, j begins at n bit multiples and is at least
    //one interaction interface ahead of i.
    index = 0;

    for(int i=0;i<L-n;i++){
      for(int j = n*(i/n)+n;j<L;j++){
	double_mutant_wout(index,0) = i;
	double_mutant_wout(index,1) = j;
	Population mutant_pop = population;
	if(single_mutant(i,0))
	  mutant_pop.genotype.aarss.iis(single_mutant(i,1),single_mutant(i,2)) = single_mutant(i,3);
	else
	  mutant_pop.genotype.trnas.iis(single_mutant(i,1),single_mutant(i,2)) = single_mutant(i,3);
	if(single_mutant(j,0))
	  mutant_pop.genotype.aarss.iis(single_mutant(j,1),single_mutant(j,2)) = single_mutant(j,3);
	else
	  mutant_pop.genotype.trnas.iis(single_mutant(j,1),single_mutant(j,2)) = single_mutant(j,3);
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = fitness_func(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);

	double_mutant_wout(index,2) = Transition_Probability(current_fitness,mutant_pop.fitness,2);
	sum += double_mutant_wout(index,2);
	index++;
      }
    }

    single_mutant.col(4) = (1/sum)*single_mutant.col(4);
    double_mutant_win.col(4) = (1/sum)*double_mutant_win.col(4);
    double_mutant_wout.col(2) = (1/sum)*double_mutant_wout.col(2);
  }

  void Evolve(){
    double rand_n = dist(mt), sum=0;
    bool transition = 0;

    Get_Mutant_Web();
    //Checking if the single mutants transition probabilities
    //exceed the random number
    
    for(int i=0;i<L;i++){
      sum+=single_mutant(i,4);
      if(sum > rand_n){
	if(single_mutant(i,0))
	  population.genotype.aarss.iis(single_mutant(i,1),single_mutant(i,2)) = single_mutant(i,3);
	else
	  population.genotype.trnas.iis(single_mutant(i,1),single_mutant(i,2)) = single_mutant(i,3);
	population.genotype.Get_Code();
	population.Get_Codon_Freq();
	population.fitness = fitness_func(&population.codon_frequency,&population.genotype.code,&population.genotype.kd);

	transition = 1;
	mutation_type="Single";
	break;
      }
    }

    //Checking if the double mutants_win transition probabilities
    //exceed the random number
    if(!transition){
      for(int i=0;i<n_d_win;i++){
	sum+=double_mutant_win(i,4);
	if(sum > rand_n){
	  if(double_mutant_win(i,0))
	    population.genotype.aarss.iis(double_mutant_win(i,1),double_mutant_win(i,2)) = double_mutant_win(i,3);
	  else
	    population.genotype.trnas.iis(double_mutant_win(i,1),double_mutant_win(i,2)) = double_mutant_win(i,3);
	  population.genotype.Get_Code();
	  population.Get_Codon_Freq();
	  population.fitness = fitness_func(&population.codon_frequency,&population.genotype.code,&population.genotype.kd);

	  transition = 1;
	  mutation_type="Double_win";
	  break;
	}
      }
    }

    //Checking if the double mutants_wout transition probabilities
    //exceed the random number
    if(!transition){
      for(int i = 0;i<n_d_wout;i++){
	sum+=double_mutant_wout(i,2);
	if(sum>rand_n){
	  if(single_mutant(double_mutant_wout(i,0),0))
	    population.genotype.aarss.iis(single_mutant(double_mutant_wout(i,0),1),single_mutant(double_mutant_wout(i,0),2)) = single_mutant(double_mutant_wout(i,0),3);
	  else
	    population.genotype.trnas.iis(single_mutant(double_mutant_wout(i,0),1),single_mutant(double_mutant_wout(i,0),2)) = single_mutant(double_mutant_wout(i,0),3);
	  if(single_mutant(double_mutant_wout(i,1),0))
	    population.genotype.aarss.iis(single_mutant(double_mutant_wout(i,1),1),single_mutant(double_mutant_wout(i,1),2)) = single_mutant(double_mutant_wout(i,1),3);
	  else
	    population.genotype.trnas.iis(single_mutant(double_mutant_wout(i,1),1),single_mutant(double_mutant_wout(i,1),2)) = single_mutant(double_mutant_wout(i,1),3);
	  population.genotype.Get_Code();
	  population.Get_Codon_Freq();
	  population.fitness = fitness_func(&population.codon_frequency,&population.genotype.code,&population.genotype.kd);

	  transition = 1;
	  mutation_type="Double_out";
	  break;
	}
      }
    }
    population.genotype.Get_Code();
    population.Get_Codon_Freq();
    population.fitness = fitness_func(&population.codon_frequency,&population.genotype.code,&population.genotype.kd);

    current_fitness = population.fitness;
  }
  
};
