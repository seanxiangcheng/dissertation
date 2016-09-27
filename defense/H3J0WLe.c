 /*
	This code is to find the density of states of the lattice gas problem 
	with at-most-1-neighbor density constraint in HN3; 
	Periodic Boundary Condition is used in this code; 
	The output is a file with number of occupied particles and their 
 	corresponding density of states
*/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/types.h>
#include<term.h>
#include<ncurses.h>
#include<gsl/gsl_rng.h>
#include<string.h>
typedef unsigned short int Natural;
#define MIN_INV_MU 0.05
#define INV_MU_STEP 0.001
#define INV_MU_L 800 

int L=0, K=0, Flat=80, Seed, Name; 	// Input parameters
float Flatness=0.8;
int StatesLength, MaxState, CurrentState, MaxStateFound, MinState=0;
double Logfnow=1.0, Logfmin=0.00000001, MinStateLogg;
double Initial_Steps=0.0, Check_Every_Steps=0.0, Total_Steps=0.0, Inner_Loop_Steps=0.0, Waste_Steps=0.0; 
double InvMu[INV_MU_L];

gsl_rng *rnd;	// RNG: random number

struct StateInfo{
  unsigned int n;
  double logg;
  float H; 
}; 


/* 	integer power of an integer (with 0^0 = 1) 	*/
int Intpow(int base, int power){
  int x=1;
  int i;
  for(i = 0; i<power; i++) {x=x*base;}
  return(x);
}


/* 	integer log2 of an integer	*/
int Intlog2(int num){
  int power=0;
  while(num>1){
    num=num/2;
    power++;
  }
  return(power);
}


/* 	Find the level i of site n and sequential number j	*/
int *Findij(int n){
  static int ij[2];
  ij[0]=1;
  ij[1]=0;
  while(n%2==0) {
    n=n/2;
    ij[0]++;
  }
  ij[1]=(n-1)/2;
  return(ij);
}



/*	Find the neighbors of site n in a network of length L with Periodic Boundary Conditions	  */
/*	The input is 1, 2, ...,L; output is a pointer to an array of 3 */
/*	e.g. Findneighbor(1, 8), returns *neighbors = [2, 3, 8]	       */
int * FindNeighbors(n){ 
  	static int neighbors[3];
  	int *ij, halfL=L/2;
  	if(n<1 || n>L){
    		fprintf(stderr,"\nError:  Cannot find the neighbor of site %d in HN3 with length %d! \n", n, L);
  	}

  	if(n==1){
    		neighbors[0]=2;
    		neighbors[1]=3;
    		neighbors[2]=L;
  	}
  	else if(n==halfL){
    		neighbors[0]=n-1;
    		neighbors[1]=n+1;
    		neighbors[2]=L;
  	}
  	else if(n==L){
    		neighbors[0]=1;
    		neighbors[1]=halfL;
    		neighbors[2]=L-1;
  	}
  	else{
    		ij=Findij(n);
    		neighbors[0]=n-1;
    		neighbors[1]=n+1;
    		if(ij[1]%2==0){
      			neighbors[2]=Intpow(2,ij[0]-1)*(2*(ij[1]+1)+1);
    		}
   		else{
      			neighbors[2]=Intpow(2,ij[0]-1)*(2*(ij[1]-1)+1);
    		}
  	}
  	return(neighbors);
}


/*	Check whether to add a particle to an empty site i	*/
int Add_or_Not(int site[], int i, int neighbors[]){
	int j, neighloc_begin, neighloc_end;
	int neig_num=0, site_1_particle=0;
	if(site[i]!=0){/*Could Delete it*/
		fprintf(stderr, "\n \n Error: The site tried to add particle has %d particle!\n", site[i]); 
        	exit(0);	
	}
	neighloc_begin=3*i;
	neighloc_end=neighloc_begin+3;
	/*	First Check its own neighbors	*/
	for(j=neighloc_begin; j<neighloc_end; j++){
		neig_num += site[neighbors[j]];
	}

	if(neig_num>=1){
		return(0);
	}
	else if(neig_num==0){
		return(1);
	}
	return(0);
}



/*	Initialize the StateInfo, neighbors and other parameters*/
void Initialize(struct StateInfo state[], int site[], int neighbors[]){
  	int i;
  	int *neigh;
  	for(i=0; i<StatesLength; i++){
		state[i].n=i;
        	state[i].logg=0.0;
		state[i].H=0.0;
  	}
  	for(i=0; i<L; i++){
		site[i]=0;
        	neigh=FindNeighbors(i+1);
        	neighbors[3*i]=neigh[0]-1;
	        neighbors[3*i+1]=neigh[1]-1;
	        neighbors[3*i+2]=neigh[2]-1;
	}
	for(i=0; i<INV_MU_L; i++){
		InvMu[i]=MIN_INV_MU + i*INV_MU_STEP;
	}
  	CurrentState=0;
  	MaxState=StatesLength-1;
  	MaxStateFound=0;
}



/*	Update the system by adding or removing particles*/
int Update(struct StateInfo state[], int site[], int neighbors[]){
	int i, n1, n2;
	i=(int)(L*gsl_rng_uniform(rnd));
	n1=CurrentState;
	if(site[i]==1){
		n2=CurrentState-1;
		if(exp(state[n1].logg-state[n2].logg)>gsl_rng_uniform(rnd) || (state[n1].logg-state[n2].logg)>200){
			site[i]=0;
			CurrentState--;
			state[CurrentState].logg += Logfnow;
			state[CurrentState].H += 1.0;		
		}
		else{
			state[CurrentState].logg += Logfnow;
			state[CurrentState].H += 1.0;
		}
	}
	else{
		n2=CurrentState+1;
		if(Add_or_Not(site, i, neighbors)){
			if(n2>MaxState){
				return(i+1);
			}
			if(exp(state[n1].logg-state[n2].logg)>gsl_rng_uniform(rnd) || (state[n1].logg-state[n2].logg)>200){
				site[i]=1;
				CurrentState++;
				state[CurrentState].logg += Logfnow;
				state[CurrentState].H += 1.0;
				if(MaxStateFound<CurrentState){
					MaxStateFound=CurrentState;
				}
			}
			else{
				state[CurrentState].logg += Logfnow;
				state[CurrentState].H += 1.0;
			}
		}
		else{
			Total_Steps -= 1.0;
			Inner_Loop_Steps -= 1.0;
			Waste_Steps +=1.0;
			return(0);
		}
	}
	return(0);
}


/*	Check the flatness of the histogram	*/
int CheckFlat(struct StateInfo state[]){
	int i;
	float meanH, totalH=0.0, minH, maxH;
	minH=state[0].H;
	maxH=state[0].H;
 	totalH=maxH;
	for(i=1; i<StatesLength; i++){
		totalH += state[i].H;
		if(state[i].H < minH){minH=state[i].H;}
		else if(state[i].H > maxH){maxH=state[i].H;}
	}
	meanH=totalH/(float)StatesLength;
	if(minH > meanH*Flatness && meanH > maxH*Flatness){
    		printf("   MC steps: %-.4e; logf=%-.2e; minH= %-.4e; maxH= %-.4e; meanH=%-.4e\n", Total_Steps/(double)L, Logfnow/2.0, minH, maxH, meanH); 
		return(1);
	}
	return(0);
}




/*	Function to read input parameters	*/
void Commandlineparse(int argc, char **argv, int *L, int *Flat,int *K, int *Seed, int *Name){
  	int i;
  	*Seed = getpid();
  	for (i = 1; i < argc; i++){  //Start at i = 1 to skip the command name.
    		if (argv[i][0] == '-'){
      			switch (argv[i][1]){
      				case 'L':       *L = atoi(argv[++i]);
        			break;
      				case 'f':       *Flat = atoi(argv[++i]);
        			break;
      				case 's':       *Seed = atoi(argv[++i]);
        			break;
      				case 'k':       *K = atoi(argv[++i]);
        			break;
      				case 'n':       *Name=atoi(argv[++i]);
        			break;
      				default:
        			fprintf(stderr,"\nError:  Incorrect option %s\n",argv[i]);
        			fprintf(stderr,"\nAvailable options: \n\
  				-L = L (HN3 length; it has to be 2^n, n=2,3,4,...)\n\
  				-k = highest level (HN3 length is 2^k))\n\
  				-f = flatness (Default = 80; or any integer in (0,99))\n\
  				-s = seed: random seed  (default=pid)\n\
    				\n");
        			exit(0);
      			}//switch
    		}//if
  	}//for
}



int main(int argc, char *argv[]){
  
  	int i, j, count=0;	// General varialbe for loops
  	unsigned int new_state_site;
  	char dos_file[15]="DOS0_";	// Name of the output file of density of states
	char pes_file[15]="PES0_";
  	char buf[8];

  	rnd =  gsl_rng_alloc(gsl_rng_mt19937);
 	gsl_rng_set(rnd, (unsigned long int)Seed);
 	for(j = 0; j <100; j++) {
      		gsl_rng_uniform(rnd); // initial cycling on the random generator
  	}

  	printf("\n    ************ HN3 Jamming Running ************ \n");
  	printf("\n      argc: %d\n      ",argc);
  	for(j = 0; j<argc; j++){
    		printf("%s ",argv[j]);
  	}

  	Commandlineparse(argc, argv, &L, &Flat, &K, &Seed, &Name);
  
  	/* Setup the simulation according to the inputs*/
  	if(L==0){
    		L=Intpow(2,K);
  	}
  	else{
    		K=Intlog2(L);
  	}
  	Flatness=(float)Flat/100.0;
  	printf("\n      L=%d; Flatness=%f;\n", L, Flatness);
  	StatesLength=(int)(L*0.8);
  	struct StateInfo *state=malloc(StatesLength*sizeof(struct StateInfo));
  	int site[L];
  	int neighbors[3*L]; 	// store the neighbors index of all sites
  	long double rho[INV_MU_L], s[INV_MU_L], shannon[INV_MU_L], rho_temp=0.0, pj;
  	long double z_temp, z[INV_MU_L];
  
  	Initial_Steps=(double)((double)K*200.0*(double)L*(double)L);
  	Check_Every_Steps=(double)((double)(K*K)*100.0*(double)L);
  
  	Initialize(state, site, neighbors);

  	Total_Steps=0.0;

  	while(Total_Steps < Initial_Steps){
		if(Update(state, site, neighbors)){
			fprintf(stderr,"\nError:  Max State greater that 0.8*L \n");    
			return(3);    
		}
        	
		Total_Steps += 1.0;
  	}

  /*printf("\n  States after Initial steps: \n n  logg        H \n");
  for(i=0; i<StatesLength; i++){
	printf(" %d  %f  %f \n", state[i].n, state[i].logg, state[i].H);
  } Could Delete it*/
  

  	printf("\n      Max State Found= %d after initial steps %.0f; \n", MaxStateFound, Initial_Steps);
  	MaxState=MaxStateFound;
  	StatesLength=MaxState+1;
  	printf("      States Length = %d; \n \n", StatesLength);
  	/* Could delete it*/
  	if(MaxState < CurrentState){
		fprintf(stderr, "\n Error: Current State %d is larger than MaxState %d!\n", CurrentState, MaxState); 
  	}/*Could delete it*/
  
  	state=(struct StateInfo *)realloc(state, StatesLength*sizeof(struct StateInfo));
  	if(state == NULL){
		fprintf(stderr, "\n Error: Memory reallocating of state(E,M,logg, H)!\n"); 
        	return(2);
  	}
  
	Inner_Loop_Steps=0.0;
  	while(Logfnow > Logfmin){
		Total_Steps += 1.0;
		Inner_Loop_Steps += 1.0;
		new_state_site=Update(state, site, neighbors);
		if(new_state_site){	// if new state found, we add another state 		
			if(CurrentState!=MaxState){
				fprintf(stderr, "\n Error: CurrentState is not the MaxState! \n"); 
        			return(3);
			}		
			MaxState++;
			CurrentState++;
			MaxStateFound=MaxState;
                	StatesLength=MaxState+1;
			state=(struct StateInfo *)realloc(state, StatesLength*sizeof(struct StateInfo));
			if(state == NULL){
        			fprintf(stderr, "\n Error: Memory reallocating of state(E,M,logg, H)!\n"); 
        			return(2);
      			}
			state[MaxState].n=MaxState;
			state[MaxState].logg=state[MaxState-1].logg-9.2;
			if(state[MaxState].logg<0){state[MaxState].logg=1.0;}
			state[MaxState].H=state[MaxState-1].H-10000;
               		if(state[MaxState].H<1){state[MaxState].H=1.0;}
                	site[new_state_site-1]=1;
			/*if(CurrentState!=MaxState){
				fprintf(stderr, "\n Error: CurrentState is not the MaxState! \n"); 
        			return(5);	
			}Could delete it*/
		}
		if(Inner_Loop_Steps > Check_Every_Steps){
			Inner_Loop_Steps = 0.0;
			if(CheckFlat(state)){
				for(i=0; i<StatesLength; i++){
					state[i].H=0.0;
				}
				Logfnow=Logfnow/2.0;
			}
		}
  	}
  
  	sprintf(buf,"%d", Name);
 	strcat(dos_file, buf);
 	strcat(pes_file, buf);

  	printf("\n      Max Occupation: %d; \n      States Length: %d; \n      Total MC steps: %-.4e; \n      Total random steps: %-.4e; \n      Wasted random steps: %-.4e", MaxState, StatesLength, Total_Steps/(double)(L), Total_Steps, Waste_Steps);
  
  	FILE *fp;
  	fp=fopen(dos_file,"w");
  	MinStateLogg=state[0].logg;
  	fprintf(fp, "n        logg \n");
  	for(i=0; i<StatesLength; i++){
		state[i].logg=state[i].logg-MinStateLogg;
		fprintf(fp, "%-8d %-18.8f\n", state[i].n, state[i].logg);
  	}
  	fclose(fp);


	for(i=0; i<INV_MU_L; i++){
		z[i]=0.0;
		rho[i]=0.0;
		s[i]=0.0;
		shannon[i]=0.0;
	}
		/*	Calculating Z, packing fraction, entropy and so on */
	for(i=0; i<INV_MU_L; i++){
		for(j=0; j<StatesLength; j++){
			z_temp = expl((long double)state[j].logg)*expl((long double)(1/InvMu[i]*(long double)j));
			z[i] += z_temp;
			rho[i] += (long double)(j*z_temp);
		}
		for(j=0; j<StatesLength; j++){
			z_temp= expl((long double)state[j].logg)*expl((long double)(1/InvMu[i]*(double)j));
			pj=(long double)(z_temp/z[i]);
			if(pj>1e-12){
				shannon[i] -= pj*log2l(pj);
			}
		}
		rho[i]=(long double)(rho[i]/z[i]/(double)L);
		s[i]=1.0/(long double)L*logl(z[i])-1/InvMu[i]*rho[i];
	}

  	fp=fopen(pes_file,"w");
	for(i=0; i<INV_MU_L; i++){
		fprintf(fp, "%-8.4f  %-18.9LG %-18.9LG  %-18.9LG \n", InvMu[i], rho[i], s[i], shannon[i]);
	}
	fclose(fp);
 	return(0);
}
