#include "atinflate_options.hpp"
#include "initialization.hpp"
#include "Stats.hpp"
#include "fitness_functions.hpp"
#include "tRNA.hpp"
#include "aaRS.hpp"
#include "Current_Genotype.hpp"
#include "Population.hpp"
#include "Evolver.hpp"

int main(int argc, char* argv[]){
  Eigen::initParallel();
  if(!initialize_variables(argc, argv))    
    return 0;

  cout<<"(a)aRS-(t)RNA (I)nteraction (N)etwork (F)itness (LA)ndscape (T)opography (E)xpress (atINFLATE)\n\n";
  clock_t t1,t2;
  t1 = clock();
  Eigen::MatrixXi t(2,T);
  Eigen::MatrixXi aiis(2,A);
  Eigen::VectorXi a(T);
  Eigen::MatrixXd Mod(S,A);
  double frac_on=0;
  
  uniform_int_distribution<int> dist_int(0,(1<<n)-1);
  bernoulli_distribution add_trans_machinery(0.001),dist_bern(p);
  ofstream traj_file(filename+"_traj.dat",std::ios_base::app);
  ofstream code_file(filename+"_code.dat",std::ios_base::app);
  ofstream prob_file(filename+"_prob.dat",std::ios_base::app);
  ofstream int_file(filename+"_int.dat",std::ios_base::app);
  ofstream codon_file(filename+"_codon.dat",std::ios_base::app);
  ofstream ocheckpoint_file(filename+"_checkpoint.log");
  /*
    Making a checkpoint file
    1  File name
    2  n
    3  k
    4  T
    5  Tcap
    6  A
    7  Acap
    8  S
    9  phi
    10 kappa
    11 N
    12 Site-type frequencies
    13 Rate
    14 proofreading
    15 mu
    16 Mu
    17 kd noncognate
    18 kd cognate
    19 Physicochemical values
    20 Codon ring space
    21 Number of bases
    22 Binomial parameter p
    23 Trajectory
    24 Halting fixation
    25 Halting fitness
  */

  ocheckpoint_file<<filename<<endl<<n<<endl<<k<<endl<<T<<endl<<Tcap<<endl<<A<<endl<<Acap<<endl<<S<<endl<<phi<<endl<<kappa<<endl<<N<<endl;
  for(int i=0;i<S;i++)
    ocheckpoint_file<<l(i)<<" ";
  ocheckpoint_file<<endl;
  ocheckpoint_file<<rate<<endl<<proofreading<<endl<<rconstant<<endl<<mu<<endl<<Mu<<endl<<kdnc<<endl<<kdc<<endl;
  for(int i = 0; i<S;i++)
    ocheckpoint_file<<sitetypes(i)<<" ";
  ocheckpoint_file<<endl;
  ocheckpoint_file<<codon_ring_space<<endl<<nbase<<endl<<p<<endl<<trial<<endl<<fxhalt<<endl<<fthalt<<endl;
  

  if(!bif){
    traj_file.close();
    traj_file.open(filename+"_traj.dat");
    code_file.close();
    code_file.open(filename+"_code.dat");
    prob_file.close();
    prob_file.open(filename+"_prob.dat");
    int_file.close();
    int_file.open(filename+"_int.dat");
    codon_file.close();
    codon_file.open(filename+"_codon.dat");
    traj_file<<"Trajectory Fixation Fitness Percent_On Mutation_Type\n";
    prob_file<<"Trajectory Fixation Site_Type Amino_Acid Probability\n";
    code_file<<"Trajectory Fixation tRNA aaRS Prob_Interaction Match kd\n";
    codon_file<<"Trajectory Fixation Site_Type Codon Codon_Frequency\n";
    int_file<<"Trajectory Fixation Molecule Number Type Value\n";
    
    tRNAStateinit.setZero();
    aaRSStateinit.setZero();
    tRNAMaskinit.setZero();
    aaRSMaskinit.setZero();
    if(b1){
      tRNAStateinit.setConstant((1<<n)-1);
      tRNAMaskinit.setConstant((1<<n)-1);
      aaRSStateinit.setConstant((1<<n)-1);
      aaRSMaskinit.setConstant((1<<n)-1);
    }
    else{
      if(!b0){
	for(int i=0;i<T;i++)
	  for(int j=0;j<trial;j++)
	    for(int intn = 0;intn<n;intn++){
	      tRNAStateinit(j,i) += (dist_bern(mt)<<intn);
	      tRNAMaskinit(j,i) += (dist_bern(mt)<<intn);
	    }
	
	for(int i=0;i<A;i++)
	  for(int j=0;j<trial;j++)
	    for(int intn = 0;intn<n;intn++){
	      aaRSStateinit(j,i) += (dist_bern(mt)<<intn);
	      aaRSMaskinit(j,i) += (dist_bern(mt)<<intn);
	    }
      }
    }
  }

  if(rate)
    fitness_func = Fitness_rate_dep;
  else
    fitness_func = Fitness_rate_indep;
  if(proofreading){
    pread=2;
  }
  else{
    pread=1;
  }
  rexponent = rconstant;
  cout<<"tRNAStateinit is\n"<<tRNAStateinit<<endl<<tRNAMaskinit<<endl<<endl<<"aaRSStateinit is\n"<<aaRSStateinit<<endl<<aaRSMaskinit<<endl;
  
  for(int trajectory = 0;trajectory<trial;trajectory++){
    for(int i=0;i<T;i++)
	  a(i) = i;
    
    //Creating tRNAs and aaRSs
    t.row(0) = tRNAStateinit.row(trajectory);
    t.row(1) = tRNAMaskinit.row(trajectory);
    tRNAs trna(t,a);
    
    //cout<<"The aaRSs are\n"<<atest.aas<<endl<<endl;
    aiis.row(0) = aaRSStateinit.row(trajectory);
    aiis.row(1) = aaRSMaskinit.row(trajectory);
    aaRSs aars(aiis,amino_acids);

    
    //Creating starting genotype with a code
    Current_Genotype current_genotype(trna, aars);
    
    //Creating starting population with a codon frequency
    Population population(current_genotype);
    
    
    Evolver evolver(population);
    
    cout<<"The tRNAs are\n";
    evolver.population.genotype.trnas.print();
    cout<<endl<<endl;
    cout<<"The aaRSs are\n";
    evolver.population.genotype.aarss.print();
    cout<<endl<<endl;
    
    cout<<"The fitness is "<<evolver.population.fitness<<endl<<endl;
    
    int fixation = endfix;


    mutation_type = "None";

    if(!bif){
      Mod = evolver.population.codon_frequency*evolver.population.genotype.code;
            traj_file<<trajectory<<" "<<fixation<<" "<<evolver.population.fitness<<" ";
      frac_on = 0;
      for(int i=0; i<T;i++)
	frac_on += __builtin_popcount(evolver.population.genotype.trnas.iis(1,i));
      for(int i=0;i<A;i++)
	frac_on += __builtin_popcount(evolver.population.genotype.aarss.iis(1,i));
      frac_on = frac_on/(n*(T+A));
      traj_file<<frac_on<<" "<<mutation_type<<endl;

      for(int trna_no=0;trna_no<T;trna_no++){
	for(int aars_no=0;aars_no<A;aars_no++){
	  
	  code_file<<trajectory<<" "<<fixation<<" "<<trna_no<<" "<<amino_acids(aars_no)<<" "<<evolver.population.genotype.code(trna_no,aars_no)<<" "<<__builtin_popcount((evolver.population.genotype.trnas.iis(1,trna_no)&evolver.population.genotype.aarss.iis(1,aars_no))&(~(evolver.population.genotype.trnas.iis(0,trna_no)^evolver.population.genotype.aarss.iis(0,aars_no))))<<" "<<evolver.population.genotype.kd(trna_no,aars_no)<<endl;
	}
      }

      for(int stype_no=0;stype_no<S;stype_no++){
	for(int codon_no=0;codon_no<T;codon_no++){
	  codon_file<<trajectory<<" "<<fixation<<" "<<sitetypes(stype_no)<<" "<<codon_no<<" "<<evolver.population.codon_frequency(stype_no,codon_no)<<endl;
	}
      }

      for(int i = 0; i<T;i++){
	int_file<<trajectory<<" "<<fixation<<" "<<"tRNA "<<i<<" State "<<evolver.population.genotype.trnas.print_trna(0,i)<<endl;
      }
      for(int i = 0; i<T;i++){
	  int_file<<trajectory<<" "<<fixation<<" "<<"tRNA "<<i<<" Mask "<<evolver.population.genotype.trnas.print_trna(1,i)<<endl;
      }
      for(int i = 0; i<A;i++){
	int_file<<trajectory<<" "<<fixation<<" "<<"aaRS "<<amino_acids(i)<<" State "<<evolver.population.genotype.aarss.print_aars(0,i)<<endl;
      }
      for(int i = 0; i<A;i++){
	int_file<<trajectory<<" "<<fixation<<" "<<"aaRS "<<amino_acids(i)<<" Mask "<<evolver.population.genotype.aarss.print_aars(1,i)<<endl;
      }


      for(int s=0;s<S;s++){
	for(int alpha=0;alpha<A;alpha++){
	  prob_file<<trajectory<<" "<<fixation<<" "<<sitetypes(s)<<" "<<amino_acids(alpha)<<" "<<Mod(s,alpha)<<"\n";
	}
      }

    }
    else
      fixation--;

    while(evolver.population.fitness < fthalt && fixation < fxhalt-1){

      evolver.Evolve();
      fixation++;

      if(add_trans_machinery(mt)){
	if(Acap<A){
	  Acap++;
	  cout<<"Fixation: "<<fixation<<" Acap: "<<Acap<<endl;
	}
	if(Tcap<T){
	  Tcap++;
	  cout<<"Fixation: "<<fixation<<" Tcap: "<<Tcap<<endl;
	}
      }

      
      Mod = evolver.population.codon_frequency*evolver.population.genotype.code;
      traj_file<<trajectory<<" "<<fixation<<" "<<evolver.population.fitness<<" ";
      frac_on = 0;
      for(int i=0; i<T;i++)
	frac_on += __builtin_popcount(evolver.population.genotype.trnas.iis(1,i));
      for(int i=0;i<A;i++)
	frac_on += __builtin_popcount(evolver.population.genotype.aarss.iis(1,i));
      frac_on = frac_on/(n*(T+A));
      traj_file<<frac_on<<" "<<mutation_type<<endl;

      for(int trna_no=0;trna_no<T;trna_no++){
	for(int aars_no=0;aars_no<A;aars_no++){
	  
	  code_file<<trajectory<<" "<<fixation<<" "<<trna_no<<" "<<amino_acids(aars_no)<<" "<<evolver.population.genotype.code(trna_no,aars_no)<<" "<<__builtin_popcount((evolver.population.genotype.trnas.iis(1,trna_no)&evolver.population.genotype.aarss.iis(1,aars_no))&(~(evolver.population.genotype.trnas.iis(0,trna_no)^evolver.population.genotype.aarss.iis(0,aars_no))))<<" "<<evolver.population.genotype.kd(trna_no,aars_no)<<endl;
	}
      }

      for(int stype_no=0;stype_no<S;stype_no++){
	for(int codon_no=0;codon_no<T;codon_no++){
	  codon_file<<trajectory<<" "<<fixation<<" "<<sitetypes(stype_no)<<" "<<codon_no<<" "<<evolver.population.codon_frequency(stype_no,codon_no)<<endl;
	}
      }

      for(int i = 0; i<T;i++){
	int_file<<trajectory<<" "<<fixation<<" "<<"tRNA "<<i<<" State "<<evolver.population.genotype.trnas.print_trna(0,i)<<endl;
      }
      for(int i = 0; i<T;i++){
	  int_file<<trajectory<<" "<<fixation<<" "<<"tRNA "<<i<<" Mask "<<evolver.population.genotype.trnas.print_trna(1,i)<<endl;
      }
      for(int i = 0; i<A;i++){
	int_file<<trajectory<<" "<<fixation<<" "<<"aaRS "<<amino_acids(i)<<" State "<<evolver.population.genotype.aarss.print_aars(0,i)<<endl;
      }
      for(int i = 0; i<A;i++){
	int_file<<trajectory<<" "<<fixation<<" "<<"aaRS "<<amino_acids(i)<<" Mask "<<evolver.population.genotype.aarss.print_aars(1,i)<<endl;
      }


      for(int s=0;s<S;s++){
	for(int alpha=0;alpha<A;alpha++){
	  prob_file<<trajectory<<" "<<fixation<<" "<<sitetypes(s)<<" "<<amino_acids(alpha)<<" "<<Mod(s,alpha)<<"\n";
	}
      }

    }
    
    
    cout<<"The Mod is\n"<<Mod<<endl;

    cout<<"The tRNAs are\n";
    evolver.population.genotype.trnas.print();
    cout<<"The aaRSs are\n";
    evolver.population.genotype.aarss.print();
    cout<<endl<<endl;
   
    for(int i=0; i<T;i++)
      ocheckpoint_file<<evolver.population.genotype.trnas.iis(0,i)<<" ";
    ocheckpoint_file<<endl;
    for(int i=0; i<T;i++)
      ocheckpoint_file<<evolver.population.genotype.trnas.iis(1,i)<<" ";
    ocheckpoint_file<<endl;
    for(int i=0; i<A;i++)
      ocheckpoint_file<<evolver.population.genotype.aarss.iis(0,i)<<" ";
    ocheckpoint_file<<endl;
    for(int i=0; i<A;i++)
      ocheckpoint_file<<evolver.population.genotype.aarss.iis(1,i)<<" ";
    ocheckpoint_file<<endl;
    
    cout<<"The fitness is "<<evolver.population.fitness<<endl<<endl;
  }
  traj_file.close();
  prob_file.close();
  code_file.close();
  codon_file.close();
  int_file.close();
  ocheckpoint_file.close();

  cout<<endl;

  t2 = clock();
  cout<<(double)(t2-t1)/CLOCKS_PER_SEC<<" seconds\n";
  cout<<(double)(t2-t1)/CLOCKS_PER_SEC/60<<" minutes\n";
  cout<<(double)(t2-t1)/CLOCKS_PER_SEC/3600<<" hours\n";

  return EXIT_SUCCESS;
}

//g++ atinflate.cpp -std=c++11 -O2 --I /PATH/TO/Eigen -o atinflate
