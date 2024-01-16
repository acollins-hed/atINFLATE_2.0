void Initialize_Environment(){
  if(codon_ring_space){
    for(int i = 0;i<T;i++){
      for(int j = 0; j<T;j++){
	if(i == j)
	  mutation(i,j) = 1 - 2*Mu;
	else{
	  if(j == i+1 || j == i-1 || (j == 0 && i == T - 1) || (j == T - 1 && i == 0))
	    mutation(i,j) = Mu;
	  else
	    mutation(i,j) = 0;
	}
      }
    }
  }
  else{
    if(nbase == 2){
      Eigen::MatrixXd Codon_Mutation(2,2);
      for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){
	  if(i == j)
	    Codon_Mutation(i,j) = 1-Mu;
	  else
	    Codon_Mutation(i,j) = Mu;
	}
      }

      if(T == 2)
	mutation = Codon_Mutation;
      if(T == 4){
	mutation = kroneckerProduct(Codon_Mutation,Codon_Mutation);
      }
      if(T == 8){
	mutation = kroneckerProduct(Codon_Mutation,kroneckerProduct(Codon_Mutation,Codon_Mutation));
      }
    }
    if(nbase == 4){
      Eigen::MatrixXd Codon_Mutation(4,4);
      
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  if(i == j)
	    Codon_Mutation(i,j) = 1-Mu;
	  else{
	    if((i == 0 && j == 1)||(i == 1 && j==0) ||(i == 2 && j == 3)||(i == 3 && j == 2))
	      Codon_Mutation(i,j) = ((double) kappa) * Mu/((double)kappa+2);
	    else
	      Codon_Mutation(i,j) = Mu/((double) kappa + 2);
	  }
	}
      }
      
      if(T == 4)
	mutation = Codon_Mutation;
      else{
	
	if(T == 16){
	  
	  mutation = kroneckerProduct(Codon_Mutation,Codon_Mutation);
	  
	}
	if(T == 64)
	  mutation = kroneckerProduct(Codon_Mutation,kroneckerProduct(Codon_Mutation,Codon_Mutation));
      }      
    }
  }
  cout<<"The Codon Mutations are\n"<<mutation<<endl;
}

