#include<iostream>
#include<fstream>
#include<string>
#include<cfloat>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<unsupported/Eigen/KroneckerProduct>
#include<cmath>
#include<time.h>
#include<random>
#include<map>

//compile with g++ -std=c++11 -O2 -I /PATH/TO/Eigen/ atinflate.cpp -o atinflate

using namespace std;

//Let
//T be the number of tRNAs
//A be the number of aaRSs
//n be the interaction interface length

unsigned long int seed;
int S,A,Acap,T,Tcap,n,k,N,L,trial,fxhalt,nbase=0,endfix=0,pread=1;
int n_d_oall;
int n_d_wout;
int n_d_win;
string filename, mutation_type="None",inputfile;
Eigen::VectorXi l;
Eigen::VectorXf amino_acids,sitetypes;
Eigen::MatrixXd W, mutation;
Eigen::MatrixXi tRNAStateinit,tRNAMaskinit, aaRSStateinit,aaRSMaskinit;
double kdc, kdnc, epsilon, mu, Mu,phi,p,fthalt,kappa,(*fitness_func)(Eigen::MatrixXd*,Eigen::MatrixXd*,Eigen::MatrixXd*),rexponent,rconstant;
random_device rd;
map<double,int> aa_to_st;
mt19937_64 mt;
uniform_real_distribution<double> dist(0,1);
bool rate=true,b0=false,b1=false,codon_ring_space=true,bif=false,proofreading=false,brconstant=false;

class BadConversion : public std::runtime_error {
public:
  BadConversion(const std::string& s)
    : std::runtime_error(s)
    { }
};

inline double convertToDouble(const std::string& s)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return x;
}


int initialize_variables(int argc, char* argv[])
{
  bool bw=false,bN=false,bu=false,bf=false,bp=false,btrna=false,baars=false,bkmax=false,bkmin=false,bns=false,bbp=false,btrial=false,bMu=false,bst=false, bfh=false,bgh=false,bcodonspace2=false,bcodonspace4=false,bkappa=false,bS=false,baarsstart=false,btrnasstart=false,bseed=false,bpcv=false;//bnu=false;
  string str_arg,sitestring,codonstring;
  uniform_real_distribution<double> aadist(0,std::nextafter(1,DBL_MAX));

  
  if(argc == 1)
    cout<<"\n\nEnter \"atinflate -h\" or \"atinflate --help\" for help\n\n";
  
  for(int i=1;i<argc;i++){
    str_arg = string(argv[i]);
    if(str_arg == "-h" || str_arg == "--help")
      {
	int help = system("echo -e 'ATINFLATE(1)			User Commands			ATINFLATE(1) \n \nNAME \n	atinflate - simulate aaRS-tRNA interaction network evolution. \n \nSYNOPSIS \n	atinflate [OPTION] \n \nDESCRIPTION \n	Initializes aaRS and tRNA interaction interfaces and evolve them through \n	a Moran Process \n	 \n	-0\n\t\tInitial genotype is all 0 genotype. Equivalent to\n\t\tatinflate --bp=0\n\t\tRandom by default.\n\n\t-1 \n\t\tInitial genotype is all 1 genotype. Equivalent to\n\t\tatinflate --bp=1\n\t\tRandom by default.\n \n\t-A=x, --aaRSs=x\n\t\twhere x is a positive integer greater than 1 representing the\n\t\tnumber of aaRSs. Default is 4.\n\n\t--bp=x\n\t\twhere x is in [0,1] representing the parameter p for a\n\t\tbinomial. Default is 0.5.\n\n\t--codon-space-x\n\t\twhere x is in {2,4} representing the number of bases.\n\t\tThis also switches from the default space which is a ring space\n\t\twhere each codon can only mutate to one of two neighbors\n\t\ton the ring.\n\t\tDefault is ring space.\n\n\t-f=x\n\t\tsets the rate selection parameter to x. Default: x = 1/44 for no\n\t\tproofreading and x = 1/9680 for proofreading.\n\n\t--fthalt=x \n\t\twhere x is in (0,1] representing the halting fitness. Default \n\t\tis 0.001.\n\n\t--fxhalt=x\n\t\twhere x is in [1,infinity) representing the halting fixation.\n\t\tDefault is 1.\n\n\t-i=file_name, --ifile=file_name\n\t\twhere file_name is the name of the input checkpoint file\n\t\tDefault: none.\n\n\t--kappa=x\n\t\twhere x is in [1,infinity) representing the transition bias.\n\t\tOnly used with --codon-space-4\n\t\tDefault is 1.\n\n\t--kmax=x\n\t\twhere x is a positive float representing the maximum\n\t\tdissociation rate. Default is 10000\n\n\t--kmin=x \n\t\twhere x is a positive float representing the minimum\n\t\tdissociation rate. Default is 220\n\n\t--mu=x \n\t\twhere x is in (0,1) representing substitution rate in the\n\t\ttranslation machinery. Default is 1e-6.\n\n\t--Mu=x\n\t\twhere x is in (0,1) representing the mutation rate between\n\t\tcodons in the entire genome. Default is 1e-4.\n\n\t-n=x\n\t\twhere x is a positive integer representing interface size.\n\t\tDefault is 4.\n\n\t-N=x, --popsize=x\n\t\twhere x is a positive integer representing population size.\n\t\tDefault is 100\n\n\t--no-rate\n\t\tMakes fitness rate independent. Fitness is rate dependent by\n\t\tdefault.\n\n\t-o=file_name, --ofile=file_name\n\t\twhere file_name is the name of the output file without \n\t\textension. Default is \"run.\"\n\n\t--phi=x\n\t\twhere x is in (0,1) representing missense tolerance, phi.\n\t\tDefault is 0.99.\n\n\t--proofreading\n\t\tapplies proofreading. Default: no proofreading.\n\n\t--seed=x\n\t\twhere x is the seed for the PRNG. Default is a random device.\n\n\t--starting-A=x \n\t\twhere x is a positive integer greater than 0 but less than or\n\t\tequal to A representing the number of beginning aaRSs.\n\t\tDefault is A.\n \n\t--starting-T=x \n\t\twhere x is a positive integer greater than 0 but less than or\n\t\tequal to T representing the number of beginning tRNAs.\n\t\tDefault is T.\n \n\t-S=x\n\t\t where x is the number of site-types. Default is A\n\n\t--Site-Types=x,y,z,...\n\t\tto customize site-type frequencies. Example, if a system has 5\n\t\tsite-types and their frequencies are 3, 4, 1, 7, 9 \n\t\trespectively, enter\n\t\t\t\tatinflate -S=5 --Site-Types=3,4,1,7,9\n\t\tDefault: all site-type frequencies are set to 1 \n\n\t--trial=x\n\t\twhere x is a positive integer representing the number of\n\t\tgenotypes to sample from the binomial for the run. Default is 1.\n \n\t-T=x, --tRNAs=x \n\t\twhere x is a positive integer greater than 1 representing the\n\t\tnumber of tRNAs. Default is 4.\n\n\t--uniform-amino-acids\n\t\tSets amino acid physicochemical values to uniform across the\n\t\tunit interval, e.g. if A=5, then the amino acid\n\t\tphysicochemical values are 0, 0.25, 0.5, 0.75, and 1.\n\t\tDefault: random across unit interval.\n \n \n \nAUTHOR \n	Written by Andrea Collins-Hed \n \nREPORTING BUGS \n	Report atinflate bugs to <acollins-hed@ucmerced.edu> \n \nCOPYLEFT \n	Copyleft (2023) A.I. Collins-Hed All Wrongs Reversed.Please cite \n	Collins-Hed et al. 2023 in published works using this software. \n \nUCM Ardell Lab	  	      	December 2023	      	 	ATINFLATE(1) \n'| less");
	if(!help)
	return 0;
      }

    if(str_arg.substr(0,7) == "--seed=")
      {
	seed = stoul(str_arg.substr(7));
	mt.seed(seed);
	bseed=true;
      }
    
    if(str_arg.substr(0,3) == "-n=")
      {
	n = stoi(str_arg.substr(3));
	if(n <= 0 ){
	  cout<<"\nn must be a positive integer. Use -h or --help for more.\n\n";
	  return 0;
	}
	bw=true;
      }

    if(str_arg.substr(0,3) == "-S=")
      {
	S = stoi(str_arg.substr(3));
	if(S < 2 ){
	  cout<<"\nS must be a positive integer > 2. Use -h or --help for more.\n\n";
	  return 0;
	}
	bS=true;
      }

    
    if(str_arg.substr(0,8) == "--trial=")
      {
	trial = stoi(str_arg.substr(8));
	if(trial < 1 ){
	  cout<<"\ntrial must be a positive integer. Use -h or --help for more.\n\n";
	  return 0;
	}
	btrial=true;
      }
    

    if(str_arg.substr(0,8) == "--kappa=")
      {
	kappa = convertToDouble(str_arg.substr(8,str_arg.length()-1));
	if(kappa < 1){
	  cout<<"\nkappa must be a positive number, at least 1. Use -h or --help for more.\n\n";
	  return 0;
	}
	bkappa=true;
      }

    
    if(str_arg.substr(0,3) == "-N=" || str_arg.substr(0,10) == "--popsize=")
      {
	if(str_arg.substr(0,3) == "-N=")
	  N = stoi(str_arg.substr(3));
	else
	  N = stoi(str_arg.substr(10));
	if(N <= 0){
	  cout<<"\nN must be a positve integer. Use -h or --help for more.\n\n";
	  return 0;
	}
	bN=true;
      }

    if(str_arg.substr(0,3) == "-A=" || str_arg.substr(0,8) == "--aaRSs=")
      {
	if(str_arg.substr(0,3) == "-A=")
	  A = stoi(str_arg.substr(3));
	else
	  A = stoi(str_arg.substr(8));
	if(A < 2)
	  {
	    cout<<"\nThe number of aaRSs must be an integer greater than 1. Use -h or --help for more.\n\n";
	    return 0;
	  }
	else
	  baars=true;
      }

    if(str_arg.substr(0,3) == "-T=" || str_arg.substr(0,8) == "--tRNAs=")
      {
	if(str_arg.substr(0,3) == "-T=")
	  T = stoi(str_arg.substr(3));
	else
	  T = stoi(str_arg.substr(8));
	if(T < 1)
	  {
	    cout<<"\nThe number of tRNAs must be an integer greater than 1. Use -h or --help for more.\n\n";
	    return 0;
	  }
	else
	  btrna=true;
	  }

    if(str_arg.substr(0,13) == "--starting-A=")
      {
	Acap = stoi(str_arg.substr(13));
	if(Acap < 1)
	  {
	    cout<<"\nThe number of beginning aaRSs must be an integer greater than 0. Use -h or --help for more.\n\n";
	    return 0;
	  }
	else
	  baarsstart=true;
      }
    
    if(str_arg.substr(0,13) == "--starting-T=")
      {
	Tcap = stoi(str_arg.substr(13));
	if(Tcap < 1)
	  {
	    cout<<"\nThe number of beginning tRNAs must be an integer greater than 0. Use -h or --help for more.\n\n";
	    return 0;
	  }
	else
	  btrnasstart=true;
      }
	
    if(str_arg.substr(0,3) == "-o=" || str_arg.substr(0,8) == "--ofile=")
      {
	if(str_arg.substr(0,3) == "-o=")
	  filename = str_arg.substr(3,str_arg.length()-1);
	else
	  filename = str_arg.substr(8,str_arg.length()-1);
	if(filename.empty()){
	  cout<<"\nDid you forget to give your output file a name? Use -h or --help for more.\n\n";
	  return 0;
	}
	bf=true;
      }

    if(str_arg.substr(0,3) == "-i=" || str_arg.substr(0,8) == "--ifile=")
      {
	if(str_arg.substr(0,3) == "-i=")
	  inputfile = str_arg.substr(3,str_arg.length()-1);
	else
	  inputfile = str_arg.substr(8,str_arg.length()-1);
	if(inputfile.empty()){
	  cout<<"\nDid you forget to give your input file a name? Use -h or --help for more.\n\n";
	  return 0;
	}
	bif=true;
      }

    if(str_arg.substr(0,21) == "--uniform-amino-acids")
      {
	bpcv=true;
      }    

    if(str_arg.substr(0,9) == "--no-rate")
      {
	rate=false;
      }

    if(str_arg.substr(0,14) == "--proofreading")
      {
	proofreading=true;
      }

    if(str_arg.substr(0,2) == "-1")
      {
	b1=true;
	p = 1;
	bbp = true;
      }

    if(str_arg.substr(0,2) == "-0")
      {
	b0=true;
	p = 0;
	bbp = true;
      }

    
    if(str_arg.substr(0,13) == "--Site-Types=")
      {
	sitestring = str_arg.substr(13,str_arg.length()-1)+",";
	bst=true;
      }

    if(str_arg.substr(0,15) == "--codon-space-2")
      {
	bcodonspace2=true;
	if(bcodonspace4){
	  cout<<"Can only have either 2 bases for codons or 4 but not both. Use --help or -h for more.\n\n";
	  return 0;
	}
      }

        if(str_arg.substr(0,15) == "--codon-space-4")
      {
	bcodonspace4=true;
	if(bcodonspace2){
	  cout<<"Can only have either 2 bases for codons or 4 but not both. Use --help or -h for more.\n\n";
	  return 0;
	}
      }
    
    if(str_arg.substr(0,5) == "--mu=")
      {
	mu=convertToDouble(str_arg.substr(5,str_arg.length()-1));
	if(mu <= 0 || mu >= 1){
	  cout<<"\nmu must be in (0,1). Use -h or --help for more.\n\n";
	  return 0;
	}
	bu=true;
      }

    if(str_arg.substr(0,5) == "--bp=")
      {
	p=convertToDouble(str_arg.substr(5,str_arg.length()-1));
	if(p < 0 || p > 1){
	  cout<<"\np must be in [0,1]. Use -h or --help for more.\n\n";
	  return 0;
	}
	bbp=true;
      }

    
    if(str_arg.substr(0,5) == "--Mu=")
      {
	Mu=convertToDouble(str_arg.substr(5,str_arg.length()-1));
	if(Mu <= 0 || Mu >= 1){
	  cout<<"\nMu must be in (0,1). Use -h or --help for more.\n\n";
	  return 0;
	}
	bMu=true;
      }
    
    if(str_arg.substr(0,9) == "--fthalt=")
      {
	fthalt=convertToDouble(str_arg.substr(9,str_arg.length()-1));
	if(fthalt <= 0 || fthalt > 1){
	  cout<<"\nfthalt must be in (0,1]. Use -h or --help for more.\n\n";
	  return 0;
	}
	bfh=true;
      }

    if(str_arg.substr(0,9) == "--fxhalt=")
      {
	fxhalt=stoi(str_arg.substr(9));
	if(fxhalt < 1){
	  cout<<"\nfxhalt must be in [1, infinity). Use -h or --help for more.\n\n";
	  return 0;
	}
	bgh=true;
      }


    if(str_arg.substr(0,7) == "--kmax=")
      {
	kdnc=float(convertToDouble(str_arg.substr(7,str_arg.length()-1)));
	if(kdnc <= 0){
	  cout<<"\nkmax must be in (0,infinity). Use -h or --help for more.\n\n";
	  return 0;
	}
	bkmax=true;
      }

        if(str_arg.substr(0,7) == "--kmin=")
      {
	kdc=float(convertToDouble(str_arg.substr(7,str_arg.length()-1)));
	if(kdc <= 0){
	  cout<<"\nkmin must be in (0,infinity). Use -h or --help for more.\n\n";
	  return 0;
	}
	bkmin=true;
      }
	if(str_arg.substr(0,3) == "-f=")
	  {
	    rconstant=convertToDouble(str_arg.substr(3,str_arg.length()-1));
	    if(rconstant <= 0 ){
	      cout<<"\nf must be in (0,infinity). Use -h or --help for more.\n\n";
	      return 0;
	    }
	    brconstant=true;
	  }
	
        if(str_arg.substr(0,6) == "--phi=")
      {
	phi=convertToDouble(str_arg.substr(6,str_arg.length()-1));
	if(phi <= 0 || phi >= 1){
	  cout<<"\nphi must be in (0,1). Use -h or --help for more.\n\n";
	  return 0;
	}
	bp=true;
      }
	if(str_arg.substr(0,7) != "--seed=" && str_arg.substr(0,6) != "--phi=" && str_arg.substr(0,7) != "--kmin=" && str_arg.substr(0,7) != "--kmax=" && str_arg.substr(0,5) != "--mu=" && str_arg.substr(0,9) != "--no-rate" && str_arg.substr(0,21) != "--uniform-amino-acids" && str_arg.substr(0,2) != "-m" && str_arg.substr(0,3) != "-o=" && str_arg.substr(0,8) != "--ofile=" && str_arg.substr(0,3) != "-T=" && str_arg.substr(0,8) != "--tRNAs=" && str_arg.substr(0,3) != "-A=" &&  str_arg.substr(0,13) != "--starting-A=" && str_arg.substr(0,13) != "--starting-T="&& str_arg.substr(0,8) != "--aaRSs=" && str_arg.substr(0,3) != "-N=" && str_arg.substr(0,10) != "--popsize=" && str_arg.substr(0,3) != "-n=" && str_arg.substr(0,2) != "-1" && str_arg.substr(0,2) != "-0" && str_arg.substr(0,8) != "--trial=" && str_arg.substr(0,5) != "--bp=" && str_arg.substr(0,5) != "--Mu=" && str_arg.substr(0,13) != "--Site-Types=" && str_arg.substr(0,9) != "--fthalt="&& str_arg.substr(0,9) != "--fxhalt="&& str_arg.substr(0,3) != "-S="&&str_arg.substr(0,15) != "--codon-space-4"&& str_arg.substr(0,15) != "--codon-space-2" && str_arg.substr(0,8) != "--kappa=" && str_arg.substr(0,3) != "-i=" && str_arg.substr(0,8) != "--ifile=" && str_arg.substr(0,14) != "--proofreading" && str_arg.substr(0,3) != "-f="){
	  cout<<endl<<str_arg<<" is not a recognized parameter. Use -h or --help for more.\n\n";
	  return 0;
	}
  }
  unsigned int random_number_from_random_device = rd();    
  if(bif){
    ifstream icheckpoint_file;
    icheckpoint_file.open(inputfile);
    if(icheckpoint_file){
      string cp;
      icheckpoint_file>>cp;
      filename = cp;
      icheckpoint_file>>cp;
      n = stoi(cp);
      icheckpoint_file>>cp;
      k = stoi(cp);
      icheckpoint_file>>cp;
      T = stoi(cp);
      icheckpoint_file>>cp;
      Tcap = stoi(cp);
      icheckpoint_file>>cp;
      A = stoi(cp);
      icheckpoint_file>>cp;
      Acap = stoi(cp);
      icheckpoint_file>>cp;
      S = stoi(cp);
      icheckpoint_file>>cp;
      phi = convertToDouble(cp);
      icheckpoint_file>>cp;
      kappa = convertToDouble(cp);
      icheckpoint_file>>cp;
      N = stoi(cp);
      l.resize(S);
      for(int i=0;i<S;i++){
	icheckpoint_file>>cp;
	l(i) = stoi(cp);
      }
      icheckpoint_file>>cp;
      rate = stoi(cp);
      icheckpoint_file>>cp;
      proofreading = stoi(cp);
      icheckpoint_file>>cp;
      rconstant = convertToDouble(cp);
      icheckpoint_file>>cp;
      mu = convertToDouble(cp);
      icheckpoint_file>>cp;
      Mu = convertToDouble(cp);
      icheckpoint_file>>cp;
      kdnc = convertToDouble(cp);
      icheckpoint_file>>cp;
      kdc = convertToDouble(cp);
      sitetypes.resize(S);
      for(int i = 0;i<S;i++){
	icheckpoint_file>>cp;
	sitetypes(i) = convertToDouble(cp);
      }
      icheckpoint_file>>cp;
      codon_ring_space = stoi(cp);
      icheckpoint_file>>cp;
      nbase = stoi(cp);
      if(!codon_ring_space && nbase == 2)
	bcodonspace2=true;
      else{
	if(nbase == 4)
	  bcodonspace4=true;
      }
      icheckpoint_file>>cp;
      p = convertToDouble(cp);
      icheckpoint_file>>cp;
      trial = stoi(cp);
      tRNAStateinit.resize(trial,T);
      aaRSStateinit.resize(trial,A);
      tRNAMaskinit.resize(trial,T);
      aaRSMaskinit.resize(trial,A);
      icheckpoint_file>>cp;
      endfix = stoi(cp);
      icheckpoint_file>>cp;
      fthalt = convertToDouble(cp);
      if(!bgh){
	cout<<"You need to use --fxhalt with -i or --ifile. See --help or -h for more.\n\n";
	return 0;
      }
      else{
	if(fxhalt < endfix){
	  cout<<"The new halting fixation, --fxhalt, must be larger than the ending fixation in the checkpoint.log file. See --help or -h for more.\n\n";
	  return 0;
	}
      }
      amino_acids.resize(A);
      for(int i=0;i<A;i++){
	amino_acids(i) = sitetypes(i);
	aa_to_st[amino_acids(i)] = i;
      }

      int index=0;
      for(int i = 0;i<4*trial;i++){
	if(i != 0 && i%4 == 0)
	  index++;
	if(i%4 == 0){
	  for(int j = 0;j<T;j++){
	    icheckpoint_file>>cp;
	    tRNAStateinit(index,j) = stoi(cp);
	  }
	}
	else{
	  if(i%4 == 1){
	    for(int j = 0;j<T;j++){
	      icheckpoint_file>>cp;
	      tRNAMaskinit(index,j) = stoi(cp);
	    }
	  }
	  else{
	    if(i%4 == 2){
	      for(int j = 0;j<A;j++){
		icheckpoint_file>>cp;
		aaRSStateinit(index,j) = stoi(cp);
	      }
	    }
	    else{
	      for(int j = 0;j<A;j++){
		icheckpoint_file>>cp;
		aaRSMaskinit(index,j) = stoi(cp);
	      }
	    }
	  }
	}
      }
    }
      
    icheckpoint_file.close();
   
  }
  else{
    if(!bseed)
      mt.seed(random_number_from_random_device);
    if(!bw)
      n=4;
    if(!bns)
      k=n;
    if(!btrna)
      T = 4;
    if(!btrnasstart)
      Tcap = T;
    if(Tcap > T){
      cout<<"The beginning number tRNAs cannot exceed the ending number of tRNAs. Use -h or --help for more.\n\n";
      return 0;
    }
    if(!baars)
      A = 4;
    if(!baarsstart)
      Acap = A;
    if(Acap > A){
      cout<<"The beginning number aaRSs cannot exceed the ending number of aaRSs. Use -h or --help for more.\n\n";
      return 0;
    }
    if(!bS)
      S = A;
    if(S < A){
      cout<<"There cannot have more aaRSs (A) than site-types (S). See --help or -h for more.\n\n";
      return 0;
    }
    sitetypes.resize(S);
    cout<<endl;
    
    if(!bN)
      N=100;
    if(!bu)
      mu=1e-6;
    if(!bMu)
      Mu=1e-4;
    if(!bkappa)
      kappa=1;
    if(!bp)
      phi=0.99;
    if(!brconstant){
      if(proofreading)
	rconstant = ((double) 1/9680);
      else
	rconstant = ((double) 1/44);
    }
    if(!bfh)
      fthalt = 0.001;
    if(!bgh)
      fxhalt = 1;
    if(!bbp)
      p=0.5;
    if(!btrial)
      trial = 1;
    tRNAStateinit.resize(trial,T);
    aaRSStateinit.resize(trial,A);
    tRNAMaskinit.resize(trial,T);
    aaRSMaskinit.resize(trial,A);
    if(!bf)
      filename="run";
    if(!bkmax)
      kdnc = 10000;
    if(!bkmin)
      kdc = 220;
    if(kdnc <= kdc){
      cout<<"\n\nkmax must be larger than kmin. Use -h or --help for more.\n\n";
      return 0;
    }
    if(k > n){
      cout<<"\nk must be at most n.\n\n";
      return 0;
    }
    if(bcodonspace2){
      nbase = 2;
      codon_ring_space=false;
      if(T != 2 && T!= 4 && T!=8){
	cout<<"\nThe number of tRNAs must be 2, 4, or 8 for a codon space with 2 bases. Use --help or -h for more.\n\n";
	return 0;
      }
    }
    if(bcodonspace4){
      nbase = 4;
      codon_ring_space=false;
      if(T != 4 && T!= 16 && T!=64){
	cout<<"\nThe number of tRNAs must be 4, 16, or 64 for a codon space with 4 bases. Use --help or -h for more.\n\n";
	return 0;
      }
      
    }
    
    l.resize(S);
    for(int i=0;i<S;i++){
      l(i) = 1;
    }
    if(bst){
      for(int i=0, j=0, stp=0;i<((int)sitestring.size());i++){
	if(stp >= S){
	  cout<<"There are more fequencies given than site-types. Use --help or -h for more.\n\n";
	  return 0;
	}
	if(sitestring.at(i)==','){
	  l(stp) = stoi(sitestring.substr(j,i-j));
	  stp++;
	  j = i + 1;
	  if(i == ((int)sitestring.size())-1 && stp < S){
	    cout<<"WARNING: there are fewer fequencies given than site-types! Use --help or -h for more.\n\n";
	  }
	}
      }
      
    }

    if(bpcv){
      for(int i=0;i<S;i++)
	sitetypes(i) = ((double) i)/((double) S-1);
    }
    else{
      for(int i=0;i<S;i++)
	sitetypes(i) = aadist(mt);
    }
    
    {
      Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> pm(S);
      pm.setIdentity();
      shuffle(pm.indices().data(),pm.indices().data()+pm.indices().size(),mt);
      sitetypes = pm*sitetypes;
      l = pm*l;
    }
    amino_acids.resize(A);
    for(int i=0;i<A;i++){
      amino_acids(i) = sitetypes(i);
      aa_to_st[amino_acids(i)] = i;
    }

  }

    
  
  ofstream log_file(filename+".log");
  time_t begin = time(0);
  tm *ltm = localtime(&begin);
  log_file<<"Started at "<<1900+ltm->tm_year<<"-"<<1+ltm->tm_mon<<"-"<<ltm->tm_mday<<" "<<ltm->tm_hour<<":"<<ltm->tm_min<<":"<<ltm->tm_sec<<endl;
  log_file<<"Saved in files: "<<filename<<endl;
  cout<<"Copyleft (2019) A.I. Collins-Hed\nAll Wrongs Reversed.\nPlease cite Collins-Hed et al. 2019 in published works using this software.\n";
  cout<<"\n--------------------------------------\nn (number of sites per interface) is "<<n<<endl;
  log_file<<"Window length: n = "<<n<<endl;
  cout<<"k (number of sites for full energy) is "<<k<<endl;
  log_file<<"Number of sites for full energy: k = "<<k<<endl;
  cout<<"T (number of tRNAs) is "<<T<<endl;
  log_file<<"Number of distinct tRNAs: T = "<<T<<endl;
  cout<<"T starting (beginning number of tRNAs) is "<<Tcap<<endl;
  log_file<<"Number of beginning tRNAs: Tcap = "<<Tcap<<endl;
  cout<<"A (number of aaRSs) is "<<A<<endl;
  log_file<<"Number of distinct aaRSs: A = "<<A<<endl;
  cout<<"A starting (beginning number of aaRSs) is "<<Acap<<endl;
  log_file<<"Number of beginning aaRSs: Acap = "<<Acap<<endl;
  cout<<"S (number of site-type) is "<<S<<endl;
  log_file<<"Number of distinct site-types: S = "<<S<<endl;
  cout<<"Site-Type frequencies are ";
  log_file<<"Site-Type frequencies: ";
  for(int i=0;i<S;i++){
    cout<<l(i);
    log_file<<l(i);
    if(i != S - 1)
      {
	log_file<<", ";
	cout<<",";
      }
  }
  cout<<endl;
  log_file<<endl;
  cout<<"The site-type physicochemical values are "<<sitetypes(0);
  log_file<<sitetypes(0);
  for(int i=1;i<S;i++){
    cout<<", "<<sitetypes(i);
    log_file<<", "<<sitetypes(i);
  }
  cout<<endl;
  log_file<<endl;
  cout<<"N (population size) is "<<N<<endl;
  log_file<<"Population size: N = "<<N<<endl;
  if(rate)
    {
      log_file<<"rate dependent\n";
      cout<<"rate dependent\n";
    }
  else
    {
      log_file<<"rate independent\n";
      cout<<"rate indepdent\n";
    }
  if(proofreading)
    {
      log_file<<"proofreading applied\n";
      cout<<"proofreading applied\n";
    }
  else
    {
      log_file<<"no proofreading\n";
      cout<<"no proofreading\n";
    }
  
  cout<<"The rate selection is "<<rconstant<<endl;
  log_file<<"Rate selection: "<<rconstant<<endl;
  
  if(b0){
    log_file<<"Initialized at all 0 genotype.\n";
    cout<<"Initialized at all 0 genotype.\n";
  }
  else{
    if(b1){
      log_file<<"Initialized at all 1 genotype.\n";
      cout<<"Initialized at all 1 genotype.\n";
    }
    else{
      log_file<<"Initialized at random genotype.\n";
      cout<<"Initialized at random genotype.\n";
      }
  }

  log_file<<"Translation machinery substitution rate: mu = "<<mu<<endl;
  cout<<"mu (translation machinery substitution rate) is "<<mu<<endl;
  log_file<<"Codon substitution rate: Mu = "<<Mu<<endl;
  cout<<"Mu (codon substitution rate) is "<<Mu<<endl;
  log_file<<"Transition bias: kappa = "<<kappa<<endl;
  cout<<"kappa (transition bias) "<<kappa<<endl;
  if(codon_ring_space){
    log_file<<"Codon ring space.\n";
    cout<<"Codon ring space.\n";
  }
  if(bcodonspace2){
    log_file<<"Codon space with 2 bases.\n";
    cout<<"Codon space with 2 bases.\n";
  }
  if(bcodonspace4){
    log_file<<"Codon space with 4 bases.\n";
    cout<<"Codon space with 4 bases.\n";
  }

  log_file<<"Baseline missense fitness: phi = "<<phi<<endl;
  cout<<"phi (the baseline fitness) is "<<phi<<endl;
  log_file<<"The halting fitness: "<<fthalt<<endl;
  cout<<"fthalt (the halting fitness) is "<<fthalt<<endl;
  log_file<<"The halting fixation: "<<fxhalt<<endl;
  cout<<"fxhalt (the halting fixation) is "<<fxhalt<<endl;
  log_file<<"The binomial parameter for genotype initialization: p = "<<p<<endl;
  cout<<"p for the binomial is "<<p<<endl;
  log_file<<"The number of trajectories: "<<trial<<endl;
  cout<<"Number of trajectories is "<<trial<<endl;
  log_file<<"The maximum dissociation rate: kmax = "<<kdnc<<endl;
  cout<<"kmax is "<<kdnc<<endl;
  log_file<<"The minimum dissociation rate: kmin = "<<kdc<<endl;
  cout<<"kmin is "<<kdc<<endl;
  log_file<<"Fixations from the Moran Process.\n";
  cout<<"Fixations from the Moran Process.\n";
  if(bseed){
    cout<<"PRNG seed: "<<seed<<endl;
    log_file<<"PRNG seed: "<<seed<<endl;
  }
  else{
    cout<<"PRNG seed: "<<random_number_from_random_device<<endl;
    log_file<<"PRNG seed: "<<random_number_from_random_device<<endl;
  }
  cout<<"Output file name is \""<<filename<<"\"\n--------------------------------------\n\n";
  log_file.close();
  L = (T+A)*2*n;

  W.resize(S,S);
  mutation.resize(T,T);

  n_d_wout = (T+A)*(2*(T+A)-1)*n*n;
  n_d_win = (n*(n-1))*(T+A);

  epsilon = (log(kdnc)-log(kdc))/double(n);
  if(proofreading && !bif){
    epsilon = 2*(log(kdnc)-log(kdc))/double(n);
    kdnc *= kdnc;
    kdc *= kdc;
  }

  
  for(int i=0;i<S;i++){
    for(int j=0;j<S;j++){
      W(i,j) = pow(phi,abs(sitetypes(i)-sitetypes(j)));
    }
  }
  
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
  
  return 1;

}
