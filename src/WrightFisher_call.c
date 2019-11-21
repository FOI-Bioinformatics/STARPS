#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define sqr(x) ((x)*(x))
#define nelem(x)  (sizeof(x) / sizeof((x)[0]))

/**************************************************************/
/* simulation program for generating populations              */
/* according to a Wright-Fisher model                         */
/* test version: use R random number generators               */
/*                                                            */
/* By Petter Lindgren and Jon Ahlinder                        */
/*                                                            */
/**************************************************************/

/* compile: gcc -o WrightFisher_v3 WrightFisher_v3.c -Wall -I/usr/include -lm -lgsl -lgslcblas */
/* in command prombt: R CMD SHLIB -lm WrightFisher_call.c */

//const gsl_rng *gBaseRand;       /* global rand number generator */

/* function declaration */
/* double *concat(double *freqG,double *fitness,int nGenotype2,int nlocitot,int *pos,int *id,int *idG,int nG); */
/* void print_genetic(int nGen,int *idG,int *pos,int *id,double *fitness,FILE *fp4); */
void get_new_id(int *idG, int *id, int *idnew, int *pos, int *posnew, int nG, int nGenotype);
int get_pos_id(int *idG, int *id, int nG, int nGenotype);
int *filter_data_int(int *idG,int *index,int k);
double *filter_data_double(double *freqG,int *index,int k);
double *filter_data_double_freq(double *freqG,int *index,int k);
int *filter_data(double *freqG,int nGenotype,int k,double thresf);
int upd_ngen(double *freqG, int nGenotype, double thresf);
int *vect_cpy_int(int *x, int n);
double *vect_cpy(double *x, int n);
/* void vect_cpy2(double *x, SEXP y, int n); */
double *vect_const_mult(int *x, double a, int n);
int *get_pos(int *pos, int nG,int nGG);
int unique_pos(int *pos, int nG);
int *make_genotype_matrix(double *xab, double *fitness,int *idG,int *id,int *pos,double *freqG,int nG,int nGenotype, int size);
int *make_genotype_matrix2(int *xab, int *idG,int *id,int *pos,int nG,int nGenotype, int size);
int find_max(double *inX, int nGenotype);
int maxpos(int *pos, int nG);
double gene_div(double *freqG, int *pos, int *idG, int *id, int nGenotype, int nG,int popsize);
double *get_freq(int *numGtmp,int nGenotype);
int *reduce_loci_2(int *x,int nGenotype, int nloci, int nloci2);
int reduce_loci_1(int *x,int nGenotype,int nloci);
double *upd_geno(double *inX,int *x,double *freqtmp,double *fitness,int *idG,int nGenotype,int nloci,int nGenotype2);
int cnt_nonz(double *freqG,int nGenotype);
void read_gen_2(double *inX,int *x, double *freqG, double *fitness, int nGenotype, int nloci);
double *read_gen_3(double *inX,double *freqG, int nGenotype, int nloci,int popsize);
int upd_geno_nGen(unsigned int *numG,int nGenotype);
int upd_geno_nG(unsigned int *numG,int *id,int *idG,int nGenotype,int nG);
double nucl_div(double *freqG, int *idG,int *id,int *pos, int nGenotype, int nG, double thres);
int differ(int *id,int *pos,int *idG,int i1,int i2,int nG);
void showint(int *M,int r, int c);
void show(double *M,int r, int c);
void freq(int *x, int r, int c);
void freq_sp(int *pos, int n1,int popsize,int nloci);
double draweff(double p1, double p2, double mu);
void fix_mut(int *x, double *mut,int r, int c, int nloci);
void freq_stat(double *freqv, int nlocitot, int nGen1, int nloci,double *muteff);
int det_sp(int *x,int popsize,int nloci);
void sparse(int *x,int *pos,int *id,int popsize,int nloci);
double walltime(double *t0);
/* int count_genotypes(int *x,int popsize,int nloci); */
/* int count_mutations(int *x,int popsize,int nloci); */
int count_mut_gen(int *x,int nGenotypes,int nloci);
/* void sparse_genotypes(int *x,int *pos,int *id,double *freqG,int nG, int nGenotype,int popsize,int nloci); */
/* void sparse_gen(int *x,int *pos,int *id,double *freqG,int nG, int nGenotype,int nloci); */
void sparse_gen(int *x,int *pos,int *id,double *freqG,int *idG,int nG, int nGenotype,int nloci,int size);
void sparse_gen_pos(int *x,int *pos,int *id,double *freqG,int *idG,int nG, int nGenotype,int nloci,double *posi);
void update_mut_freq(int *pos,int *id,double *freqG,int *idG,int nG,int nGenotype, int nloci,int gen,FILE *fp);
double *calc_gen_prob(double *freqG,double *fitness,int nGenotype,int popsize);
double mean_fitness(double *fitness,double *freqG, int nGenotype);

SEXP WF2(SEXP nGenp, SEXP popsizep, SEXP expfactorp, SEXP mutratep, SEXP nlocip, SEXP verbosityp, SEXP nGenotypep, SEXP flag, SEXP inXp, SEXP posp, SEXP thresp, SEXP thresfrp, SEXP fitneffp){
  
  // nloci, nGenotype can be omitted in function call?
  // is expfactor an integer??

  // defining input variables
  int size = length(fitneffp);  
  int nGen=INTEGER_VALUE(nGenp), popsize=INTEGER_VALUE(popsizep), expfactor=INTEGER_VALUE(expfactorp), nloci=INTEGER_VALUE(nlocip), verbosity=INTEGER_VALUE(verbosityp), flagg=INTEGER_VALUE(flag), nGenotype=INTEGER_VALUE(nGenotypep);
  double mutrate = NUMERIC_VALUE(mutratep), thres=NUMERIC_VALUE(thresp), thresf=NUMERIC_VALUE(thresfrp);
  double *fitneff = REAL(fitneffp);
  double *inX = REAL(inXp);
  double *posi = REAL(posp);
  // double p1 = NUMERIC_VALUE(p1p), p2 = NUMERIC_VALUE(p2p), mu = NUMERIC_VALUE(mup),
  // defining additional variables
  int ntot=0,i,j,k,nmut=0,nlocitot,popsize_i,popsize_i0,l1,tmp1,i1;
  int *x,*pos,*id,*start,*end,*nrmut,*indvec,totpop,cnt,flag3;
  double *fitness,u1,*freqG,*freqtmp,*divers,*gdiv,*fitn;  

  unsigned int cmptmp;
  int cnt1=0,cnt2=0;
  int nG; /* number of mutations present in the genotypes */
  unsigned int *numG,*numGtmp;
  int *idtmp,*idtmp2,*postmp,*postmp2, *idG, *idGtmp,*tmpID,*xnew,*id_cand,*pos_cand, *pos_anc,*numGtmp2;
  int nG2,nGenotype2,l,tmp0,nloci2,j1,ntot2,nGtmp,maxnr;
  double *fitntmp, *pG,tmpmut,mf,pi,*freqtmp2,gd, *freqtmpid,*fitntmp2;
  double tmpu,prod;
  int ul = 1E6; /* upper limit */
  int ntmp=0,u2=0,u,mem,newgen,m;
  int *posm,*posm2,*numID,*unpos;
  int km,up;

  /* defining time variables */
  double startTime,elapsedTime;
  double clockZero = 0.0;
  
  GetRNGstate(); // initialize random number generator
  
  //PROTECT(fitneffp = AS_NUMERIC(fitneffp));

  /* error check user input */
  /* if(p1 > 1. || p1 < 0.){ */
  /*   Rprintf("probabilities for detrimental mutations cannot exceed 1.0 or be less than 0.0\n"); */
  /*   exit(-1); */
  /* } */
  /* if(p2 > 1. || p2 < 0.){ */
  /*   Rprintf("probabilities for neutral mutations cannot exceed 1.0 or be less than 0.0\n"); */
  /*   exit(-1); */
  /* } */
  /* if((p2+p1) > 1.){ */
  /*   Rprintf("probabilities for neutral and deterimental mutations cannot exceed 1.0\n"); */
  /*   exit(-1); */
  /* } */
  /* set help variables */

  start = (int *) R_alloc(nGen,sizeof(int));
  end = (int *) R_alloc(nGen,sizeof(int));
  ntot = popsize;

  for(i = 1; i <= nGen; i++){
    j = i-1;
    start[j] = ntot;
    ntot += popsize*pow(expfactor,i);
    end[j] = ntot-1;
  }
  if(verbosity==1) Rprintf("Total number of generated bacteria: %d\n",ntot);
  if(verbosity==1) Rprintf("generate mutations...\n");
  /* generate mutations */
  /* loop over all bacteria after the initial generation in the population */
  ntot2 = ntot-start[0];
 
  nrmut = (int *) calloc(nGen,sizeof(int)); // ok
  //nrmut = (int *) R_alloc(nGen,sizeof(int));
  //posm = (int *) malloc(ntot2*size*sizeof(int)); 
  posm = (int *) calloc(ntot2,sizeof(int)); // ok

  for(i = 0; i < ntot2; i++){ // loop over generated bacterias
    ntmp = rbinom(size, mutrate);//gsl_ran_binomial(gBaseRand, mutrate, size); /* generate mutations */
    if(ntmp==1){ // if a mutation has occur
      tmpu = runif(0.0,1.0);
      prod = tmpu*size;
      km = (int)prod; //gsl_rng_uniform_int (gBaseRand, size); // generate a position from a uniform distribution 
      posm[nmut] = km; // save position
      nmut++; /* count the number of mutations */   
      tmp0 = i + start[0];
      for(j = 0; j < nGen; j++)
	if(tmp0 >= start[j] && tmp0 <= end[j])
	  nrmut[j]++;
    }
  }

 if(nmut > 0) posm2 = (int *) R_alloc(nmut,sizeof(int));
 else posm2 = (int *) R_alloc(1,sizeof(int));
  for(i = 0; i < nmut; i++) // copy mutation vectors to save memory 
    posm2[i] = posm[i];
 
  free(posm); 

  if(verbosity==1) printf("Distribution of mutations over generations:\n");
  if(verbosity==1) showint(nrmut,1,nGen);
  if(verbosity==1) printf("\nNumber of mutations: %d\n",nmut);

  /* allocate memory */
  freqG = (double *) calloc(nGenotype,sizeof(double)); // ok

  if(verbosity==1) printf("nloci = %d\tnGenotype = %d\n",nloci,nGenotype);

  freqtmp = read_gen_3(inX,freqG,nGenotype,nloci,popsize); /* read input from R */ 
  /* Rprintf("Frequencies:\n"); */
  /* show(freqG,1,nGenotype); */
  maxnr = find_max(inX,nGenotype);
  Rprintf("thresf: %f\n",thresf);
 
  free(freqG);
 
  /* reduce data from non-existing clones */
  nGenotype2 = cnt_nonz(freqtmp,nGenotype);
  if(verbosity==1) printf("New number of genotypes after removing non-existant variants: %d\n",nGenotype2);
  /* allocating memory */
  xnew = (int *) calloc(nGenotype2*nloci,sizeof(int)); // ok
  idG = (int *) calloc(nGenotype2,sizeof(int)); // ok 
  fitness = (double *) calloc(nGenotype2,sizeof(double)); // ok
  freqG = upd_geno(inX,xnew,freqtmp,fitness,idG,nGenotype,nloci,nGenotype2);
  printf("xnew:\n");
  showint(xnew,nGenotype2,nloci);
  /* free(inX); */
  free(freqtmp);
 
  nGenotype = nGenotype2;
  nG = count_mut_gen(xnew,nGenotype,nloci);
  nloci2 = reduce_loci_1(xnew,nGenotype,nloci);
  x = reduce_loci_2(xnew,nGenotype,nloci,nloci2);

  printf("Genotypes:\n");
  showint(x,nGenotype,nloci2);
 
  free(xnew);
  nloci = nloci2;
  /* nG = count_mutations(x,popsize,nloci); */
  if(verbosity==1) printf("nG = %d\tnGenotype = %d\tnloci = %d\n",nG,nGenotype,nloci);
  printf("nG = %d\tnGenotype = %d\tnloci = %d\n",nG,nGenotype,nloci);
  pos = (int *) calloc(nG,sizeof(int)); // ok
  id = (int *) calloc(nG,sizeof(int)); // ok
  popsize_i0 = popsize; /* starting population size prior to expansion */

  //sparse_gen(x,pos,id,freqG,idG,nG, nGenotype,nloci,size);
  sparse_gen_pos(x,pos,id,freqG,idG,nG, nGenotype,nloci,posi);
  /* free(posi); */
//sparse_gen_pos(int *x,int *pos,int *id,double *freqG,int *idG,int nG, int nGenotype,int nloci,int *posi);
  if(verbosity==1) Rprintf("id:\n");
  if(verbosity==1) showint(id,1,nG);
  if(verbosity==1) Rprintf("pos:\n");
  if(verbosity==1) showint(pos,1,nG);
  /* Rprintf("id:\n"); */
  /* showint(id,1,nG); */
  /* Rprintf("pos:\n"); */
  /* showint(pos,1,nG); */

  free(x);
  /* simulate cell division */

  pi = nucl_div(freqG,idG,id,pos, nGenotype, nG, thres);
  Rprintf("Initial pair-wise nucleotide diversity: %f\n",pi);
  if(nG > 0){
    gd = gene_div(freqG, pos, idG, id, nGenotype, nG, popsize);
    Rprintf("Initial average gene diversity: %f\n",gd);
  }
  if(flagg==1){
    if(!(divers = (double *) calloc((nGen+1),sizeof(double))) ||
       !(fitn = (double *) calloc((nGen+1),sizeof(double))) ||
       !(gdiv = (double *) calloc((nGen+1),sizeof(double)))){ 
	printf("Error in allocating memory at line 244\n\n"); 
	return NULL;
      }
  }
  /* fprintf(fp5,"%d\t%f\t%f\n",0,1.,pi);  */
  /* if(flagg==1) divers[0] = pi; */
  /* if(flagg==1) gdiv[0] = gd; */
  /* fitn[0] = mean_fitness(fitness,freqG,nGenotype); */
  Rprintf("Initial mean fitness: %f\n",mean_fitness(fitness,freqG,nGenotype));
  totpop = maxnr+1; /* set population counter */
  startTime = walltime(&clockZero); /* start clock */
  nlocitot = nloci; /* current number of segregating loci in population */

  for(i = 1; i <= nGen; i++){
    Rprintf("********** Simulating generation %d **********\n\n",i);
    popsize_i = popsize*pow(expfactor,i); /* population size at generation i */

    i1 = i-1;
    if(verbosity==1) printf("select clones to reproduce...\n");
    /* select clones to reproduce */
    /* calculate sample probabilities for each genotype */
    pG = calc_gen_prob(freqG,fitness,nGenotype,popsize_i0);
    /* Rprintf("pG\n"); */
    /* show(pG,1,nGenotype); */
    free(freqG);
    Rprintf("Number of genotypes in population prior to selection (nGenotype): %d\nNumber of mutations in population (nG): %d\npopsize at generation %d: %d\n",nGenotype,nG,i,popsize_i);

    if(!(numG = (unsigned int *) calloc(nGenotype,sizeof(int)))){ 
	printf("Error in allocating memory at line 259\n\n"); 
	return NULL;
      } // ok
    if(verbosity==1) printf("sample from multinomial distibution to simulate new genotype frequencies\n");
    /* sample from multinomial distibution to simulate new genotype frequencies */
    
    //rmultinom(int n, double* prob, int k, int* rn)
    rmultinom(popsize_i, pG, nGenotype, numG); //gsl_ran_multinomial (gBaseRand, nGenotype, popsize_i, pG, numG);
    if(verbosity==1) Rprintf("Number of genotypes in population prior to selection: %d\n",nGenotype);
    /* Rprintf("Genotype probabilities for selection:\n"); */
    /* show(pG,1,nGenotype); */

    free(pG);

    /* Rprintf("head(idG):\n"); */
    /* if(nGenotype > 10) */
    /*   showint(idG,1,10); */
    /* else */
    /*   showint(idG,1,nGenotype); */
    /* Rprintf("head(id):\n"); */
    /* if(nG > 10) */
    /*   showint(id,1,10); */
    /* else */
    /*   showint(id,1,nG); */
    /* Rprintf("head(numG):\n"); */
    /* if(nGenotype > 10) */
    /*   showint(numG,1,10); */
    /* else */
    /*   showint(numG,1,nGenotype); */
    /* remove unselected genotypes */ // NB NB can it be postponed until updating later??
    nG2 = upd_geno_nG(numG,id,idG,nGenotype,nG); /* update number of mutations */
    nGenotype2 = upd_geno_nGen(numG,nGenotype); /* update number of genotypes */
  
    Rprintf("Number of genotypes in population after selection (nGenotype2): %d\nNumber of new mutations in population (nG2): %d\n",nGenotype2,nG2);
    /*   printf("id:\n"); */
    /* showint(id,1,nG); */
    /* allocate memory */
    /* tmp1 = nGenotype2; */
    if(verbosity==1) Rprintf("number of Genotypes after selection (%d) + new mutations (%d): tmp1 = %d\n",nGenotype2,nrmut[i1],tmp1);
    /* Rprintf("\nAllocating memory for:\nnumGtmp\t%d\nidGtmp\t%d\nfitntmp\t%d\nidtmp\t%d\npostmp\t%d\n",nGenotype2,nGenotype2,nGenotype2,nG2,nG2); */
    if(!(numGtmp = (unsigned int *) calloc(nGenotype2,sizeof(int))) || // ok
       !(idGtmp = (int *) calloc(nGenotype2,sizeof(int))) || // ok
       !(fitntmp = (double *) calloc(nGenotype2,sizeof(double))) || // ok
       !(idtmp = (int *) calloc(nG2,sizeof(int))) ||// ok 
       !(postmp = (int *) calloc(nG2,sizeof(int)))){ 
	printf("Error in allocating memory at line 282\n\n"); 
	return NULL;
      } // ok 
    if(verbosity==1) printf("numG:\n");
    if(verbosity==1) showint(numG,1,nGenotype);
    if(verbosity==1) Rprintf("id:\n");
    if(verbosity==1) showint(id,1,nG);
    /* if(verbosity==1) Rprintf("postmp:\n"); */
    /* if(verbosity==1) showint(postmp,1,nG2); */
    if(verbosity==1) printf("idG:\n");
    if(verbosity==1) showint(idG,1,nGenotype);
    if(verbosity==1) printf("numG:\n");
    if(verbosity==1) showint(numG,1,nGenotype);
    k = 0; /* reset counters */
    l1 = 0;
    /* copy data */
    for(j = 0; j < nGenotype; j++){
      if(numG[j] > 0){  // copy
	numGtmp[k] = numG[j];
	idGtmp[k] = idG[j];
	fitntmp[k] = fitness[j];
	if(fitntmp[k]<0.) Rprintf("NB! WARNING! fitntmp[%d] = %f\n",k,fitntmp[k]);
	k++;
	for(l = 0; l < nG; l++)
	  if(id[l] == idG[j]){ // copy // NB NB update!!
	    idtmp[l1] = id[l];
	    postmp[l1] = pos[l];
	    l1++;
	  }
      }
    }
    /* Rprintf("************ 1 ***************\nnumber of elements accessed in loop:\n(%d of %d) for nGenotype\n(%d of %d) for nG\n",k,nGenotype2,l1,nG2); */
    nG2 = l1;  /* update number of positions in pos and id vectors */
    Rprintf("number of copied points: %d, nGenotype2: %d, nGenotype: %d, nG2: %d\n",k,nGenotype2,nGenotype,nG2);

    newgen = maxpos(idG,nGenotype) + 1;
    Rprintf("newgen = %d, max(idG) = %d\n",newgen,maxpos(idG,nGenotype));
        
    free(fitness);
    free(id);  
    free(idG);
    free(pos);
    free(numG);

    if(verbosity==1) Rprintf("idtmp:\n");
    if(verbosity==1) showint(idtmp,1,nG2);
    /* if(verbosity==1) Rprintf("postmp:\n"); */
    /* if(verbosity==1) showint(postmp,1,nG2); */
    if(verbosity==1) printf("idGtmp:\n");
    if(verbosity==1) showint(idGtmp,1,nGenotype2);
    if(verbosity==1) printf("numGtmp:\n");
    if(verbosity==1) showint(numGtmp,1,nGenotype2);
    if(verbosity==1) printf("distribute new mutations over remaining genotypes\n");
    /* Rprintf("head(idGtmp):\n"); */
    /* if(nGenotype2 > 10) */
    /*   showint(idGtmp,1,10); */
    /* else */
    /*   showint(idGtmp,1,nGenotype2); */
    /* Rprintf("head(idtmp):\n"); */
    /* if(nG2 > 10) */
    /*   showint(idtmp,1,10); */
    /* else */
    /*   showint(idtmp,1,nG2); */
    /* Rprintf("head(numGtmp):\n"); */
    /* if(nGenotype2 > 10) */
    /*   showint(numGtmp,1,10); */
    /* else */
    /*   showint(numGtmp,1,nGenotype2); */
    /* distribute new mutations over remaining genotypes */
    if(nrmut[i1] > 0){
      /* Rprintf("\nAllocating memory for:\nfreqtmp2\t%d\nnumG\t%d\n",tmp1,tmp1); */
      freqtmp2 = get_freq(numGtmp,nGenotype2);
      if(!(numG = (unsigned int *) calloc(nGenotype2,sizeof(int)))){ 
	printf("Error in allocating memory at line 340\n\n"); 
	return NULL;
      } // ok
      if(verbosity==1) Rprintf("sample from multinomial distibution to simulate new genotype frequencies\n");
      if(verbosity==1) Rprintf("freqtmp:\n");
      if(verbosity==1) show(freqtmp2,1,nGenotype2);
      Rprintf("nrmut[%d] = %d\n",i1,nrmut[i1]);
    /* sample from multinomial distibution to get select clones to mutate */
      rmultinom(nrmut[i1], freqtmp2, nGenotype2, numG);
      free(freqtmp2);
 
      if(verbosity==1) Rprintf("step 2: assign mutations to genotypes...\n");
      /* step 2: assign mutations to genotypes */
      u = 0; /* counter for the number of mutations assigned */


      if(verbosity==1) Rprintf("Looping over selected genotypes...\n");
      if(verbosity==1) Rprintf("numG:\n");
      if(verbosity==1) showint(numG,1,nGenotype2);
      mem = 0; // find max memory size
      for(j = 0; j < nGenotype2; j++ ){
	l1 = 0;  
	for(l = 0; l < nG2; l++ ){
	  if(idtmp[l] == idGtmp[j]) l1++;	  
	}
	if(l1 > mem) mem = l1;
      }
      
      /* Rprintf("************ 2 ***************\nnumber of elements accessed in loop:\n(%d of %d) for nGenotype\n(%d of %d) for nG\n",j,nGenotype2,l,nG2); */
      //Rprintf("Max number of mutations per genotype: %d\n",mem);
      /* Rprintf("\nAllocating memory for:\nid_cand\t%d\npos_cand\t%d\npos_anc\t%d\nnumGtmp2\t%d\nfitntmp2\t%d\n",(mem+1)*nrmut[i1],(mem+1)*nrmut[i1],mem,nGenotype2,(mem+1)*nrmut[i1]); */
      if(!(id_cand = (int *) calloc((mem+1)*nrmut[i1],sizeof(int))) || // allocate memory // ok
	 !(pos_cand = (int *) calloc((mem+1)*nrmut[i1],sizeof(int))) || // ok 
	 !(numGtmp2 = (int *) calloc(nGenotype2,sizeof(int))) || // ok
	 !(fitntmp2 = (double *) calloc((mem+1)*nrmut[i1],sizeof(double)))){ 
	printf("Error in allocating memory at line 369\n\n"); 
	return NULL;
      }
	   // ok
      if(mem > 0) pos_anc = (int *) calloc(mem,sizeof(int));  // allocate memory // ok
      else pos_anc = (int *) calloc(1,sizeof(int));
      k = 0; /* reset counters */
      for(j = 0; j < nGenotype2; j++ ){
	/* Rprintf("loop over genotype %d\n",j); */
	fitntmp2[j] = fitntmp[j];
	
	numGtmp2[j] = numGtmp[j];
	//Rprintf("numGtmp:\n");
	//showint(numGtmp,1,nGenotype2);

	if(numG[j] > 0){

	  /* sample from multinomial distibution to get specific id of selected clones that mutate */
	  l1 = 0; // reset number of mutations at genotype j counter
	  //Rprintf("pos_anc:\n");
	  for(l = 0; l < nG2; l++ ){ // copy ancestral genotype
	    if(idtmp[l] == idGtmp[j]){
	      pos_anc[l1] = postmp[l]; 
	      //  Rprintf("pos_anc[%d] = %d ",l1,pos_anc[l1]);
	      //pos_ind_anc[l1] = l;
	      l1++; 
	    }
	  }
	  /* Rprintf("************ 3 ***************\nnumber of elements accessed in loop:\n(%d of %d) for mem\n",l1,mem); */
	  //printf("\n");
	  //Rprintf("\nAllocating memory for:\nfreqtmpid\t%d\nnumID\t%d\n",numGtmp[j],numGtmp[j]);
	  if(!(freqtmpid = (double *) calloc(numGtmp[j],sizeof(double))) || // ok
	     !(numID = (unsigned int *) calloc(numGtmp[j],sizeof(int)))){ 
	    printf("Error in allocating memory at line 413\n\n"); 
	    return NULL;
	  } // ok
	  
	  for(l = 0; l < numGtmp[j]; l++ ) freqtmpid[l] = (double)1./numGtmp[j];	  
	  rmultinom(numG[j], freqtmpid, numGtmp[j], numID);
	  for(l = 0; l < numGtmp[j]; l++ ){
	    
	    if(numID[l] > 0){ // genotype j mutated
	      flag3 = 0;
	      //Rprintf("numGtmp[%d] = %d, numID[%d] = %d\n",j,numGtmp[j],l,numID[l]);
	      if(numID[l] > 1) Rprintf("%d mutations on the same copy!\n",numID[l]);
	      fitntmp2[k] = fitntmp[j]; // copy ancestral fitness
	      // check if mutation is reverted!!!!!
              for(m = 0; m < l1; m++){
	      	if(pos_anc[m]==posm2[u2]){
	      	  pos_anc[m] = -9;
	      	  flag3 = 1;
	      	  Rprintf("Reverted mutation detected for clone id: %d, position: %d!\n",newgen,posm2[u2]);
	      	}
	      }

	      // save candidate
	      for(m = 0; m < l1; m++){  // ancestral contribution
		id_cand[u] = newgen;
		pos_cand[u] = pos_anc[m];
		u++;
	      }
	      for(m = 0; m < numID[l]; m++){        // new mutation
		id_cand[u] = newgen;
		if(flag3==0){
		  pos_cand[u] = posm2[u2];
		  fitntmp2[k] += fitneff[posm2[u2]];
		  if(fitntmp2[k] < 0.)
		    Rprintf("1.WARNING! Negative fitness! u = %d, u2 = %d, k = %d, posm2[%d] = %d, fitneff[%d] = %f\n",u,u2,k,u2,posm2[u2],posm2[u2],fitneff[posm2[u2]]);
		}
		else{
		  pos_cand[u] = -9;
		  fitntmp2[k] -= fitneff[posm2[u2]];
		  if(fitntmp2[k] < 0.)
		    Rprintf("2.WARNING! Negative fitness! u = %d, u2 = %d, k = %d, posm2[%d] = %d, fitneff[%d] = %f\n",u,u2,k,u2,posm2[u2],posm2[u2],fitneff[posm2[u2]]);
		}
		u++;
		u2++;
	      }
	     
	      numGtmp2[j]--; // remove copy from genotype j in temporary vector
	      /* Rprintf("fitn effect: %f at pos %d for mutation no %d\nAccessing position %d in fitntmp2\n",fitneff[posm2[u3]],posm2[u3],u3,tmp1+k); */
	      newgen++;  
	      k++;	  
	    }
	    
	  }
	  /* Rprintf("\nLooping over: 1 fitntmp2[%d] 2 fitntmp[%d] 3 numGtmp2[%d] 4 numGtmp[%d] 5 pos_ans[%d] 6 postmp[%d] 7 idGtmp[%d] 8 idtmp[%d] 9 id_cand[%d] 10 pos_cand[%d] 11 posm2[%d]\n",k,tmp1,tmp1,tmp1,l1,nG2,tmp1,nG2,u,u,u3); */
	  
	  free(freqtmpid);
	  free(numID);
	}
	
      }
      /* Rprintf("************ 4 ***************\nnumber of elements accessed in loop:\n(%d of %d) for u\n(%d of %d) for u2\n(%d of %d) for u3\n(%d of %d) for k\n",u,(mem+1)*nrmut[i1],u2,nmut,u3,nmut,k,(mem+1)*nrmut[i1]); */
      /* Rprintf("Hit1\n"); */
      free(numG);
      free(pos_anc);
      free(numGtmp);

      // add code here for the alternative scenarios of rejecting condidates
      //NBNBNB! If all candidates are accepted.
      /* allocate memory for new mutation parameters */
      nG = u + nG2; // new size of pos and id vectors
      /* Rprintf("Hit2\n"); */
      /* Rprintf("\nAllocating memory for:\nid\t%d\npos\t%d\nnumG\t%d\nidG\t%d\nfreqG\t%d\nfitness\t%d\n",nG,nG,newgen,newgen,newgen,newgen); */
      if(!(id = (int *) calloc(nG,sizeof(int))) ||
	 !(pos = (int *) calloc(nG,sizeof(int))) ||
	 !(numG = (int *) calloc(newgen,sizeof(int))) ||
	 !(idG = (int *) calloc(newgen,sizeof(int))) ||
	 !(freqG = (double *) calloc(newgen,sizeof(double))) ||
	 !(fitness = (double *) calloc(newgen,sizeof(double)))){ 
	    printf("Error in allocating memory at line 484\n\n"); 
	    return NULL;
	  } // ok
      /* merge ancestral and new genotypes and remove zero contributing ancestral genotypes (highly unlikely though) */
      m = 0;
      /* Rprintf("Hit3\n"); */
      for(j = 0; j < nG2; j++ ){ // copy old clones
	flag3 = 0;
	/* NB! reduce the number of the original clone */
	for(l = 0; l < nGenotype2; l++ )
	  if(idGtmp[l] == idtmp[j])
	    if(numGtmp2[l] > 0)
	      flag3 = 1;	  
	if(flag3 = 1){
	  id[m] = idtmp[j];
	  pos[m] = postmp[j];
	  m++;
	}
      }
      /* Rprintf("Hit4\n"); */
      free(idtmp);
      free(postmp);
      for(j = 0; j < u; j++ ){ // new candidate
        if(pos_cand[j] != -9){
	  id[m] = id_cand[j];
	  pos[m] = pos_cand[j];
	  m++;
	}
      }
      /* Rprintf("************ 5 ***************\nnumber of elements accessed in loop:\n(%d of %d) for nG\n",m,nG); */
      /* Rprintf("Hit5\n"); */
      free(pos_cand);
      l = 0;
      for(j = 0; j < nGenotype2; j++ )
	if(numGtmp2[j] > 0){
	  numG[l] = numGtmp2[j];
	  idG[l] = idGtmp[j];
	  fitness[l] = fitntmp[j];
	  freqG[l] = (double)numGtmp2[j]/popsize_i; 
	  l++;
	}

      free(fitntmp);
      l1 = id_cand[0];
      /* Rprintf("Hit6\n"); */
      for(j = 0; j < k; j++ ){
	idG[l] = l1 + j;  
	numG[l] = 1;
	fitness[l] =  fitntmp2[j];
	freqG[l] = (double)1./popsize_i; 
	l++;	
      }
      /* Rprintf("************ 6 ***************\nnumber of elements accessed in loop:\n(%d of %d) for nGenotype\n",l,newgen); */
      free(id_cand);
      free(fitntmp2);
      free(numGtmp2);
      free(idGtmp);
      nGenotype2 = l;
      /* Rprintf("Hit7\n"); */
    }
    else{ // no mutations at current generation
      // allocating memory
      // copy mutation data
      idG = vect_cpy_int(idGtmp, nGenotype2);
      numG = vect_cpy_int(numGtmp, nGenotype2);
      fitness = vect_cpy(fitntmp, nGenotype2);
      id = vect_cpy_int(idtmp, nG2);
      pos = vect_cpy_int(postmp, nG2);
      freqG = vect_const_mult(numGtmp, 1./popsize_i, nGenotype2);

      nG = nG2;
      // free memory
      free(fitntmp);
      free(idtmp);  
      free(idGtmp);
      free(postmp);   
      free(numGtmp);
    }

    if(verbosity==1) Rprintf("\nfreqG:\n");
    if(verbosity==1) show(freqG,1,nGenotype2);
    if(verbosity==1) Rprintf("\nfitness:\n");
    if(verbosity==1) show(fitness,1,nGenotype2);
    if(verbosity==1) Rprintf("idG:\n");
    if(verbosity==1) showint(idG,1,nGenotype2);
    if(verbosity==1) Rprintf("numG:\n");
    if(verbosity==1) showint(numG,1,nGenotype2);
    if(verbosity==1) Rprintf("id:\n");
    if(verbosity==1) showint(id,1,nG);
    /* if(verbosity==1) Rprintf("pos:\n"); */
    /* if(verbosity==1) showint(pos,1,nG); */

    popsize_i0 = popsize_i;
    nlocitot += nrmut[i1];
    nGenotype = nGenotype2;

    Rprintf("nGenotype: %d, nG: %d\n",nGenotype,nG);
 
    if(verbosity==1) Rprintf("update mean fitness\n");
    mf = mean_fitness(fitness,freqG,nGenotype);
    if(verbosity==1) printf("update nucleotide diversity\n");
    if(flagg==1){ 
      pi = nucl_div(freqG,idG,id,pos, nGenotype, nG, thres);
      gd = gene_div(freqG, pos, idG, id, nGenotype, nG, popsize_i);
      divers[i] = pi;
      gdiv[i] = gd;
      fitn[i] = mf;
    }
    if(verbosity==1) Rprintf("Average fitness: %f\n",mf);
    if(verbosity==1 & flagg==1) Rprintf("Average pair-wise nucleotide diversity: %f\nNeis gene diversity: %f\n",pi,gd);
  /* fprintf(fp5,"%d\t%f\t%f\n",i,mf,pi);  */
  }

  elapsedTime = walltime(&startTime)/60.0;

  /* if(verbosity==1)  */
  Rprintf("Elapsed time (in minutes): %f\n",elapsedTime);
  Rprintf("Average fitness after %d generations: %f\n",nGen,mf);
  if(flagg==1) Rprintf("Average pair-wise nucleotide diversity after %d generations: %f\n",nGen,pi);
  PutRNGstate();
  nGenotype2 = upd_ngen(freqG, nGenotype, thresf);
  Rprintf("New number of genotypes in population after filtering: %d (was %d)\n",nGenotype2,nGenotype);
  /* if(verbosity==1) Rprintf("\nold freqG:\n"); */
  /* if(verbosity==1) show(freqG,1,nGenotype); */
  /* Rprintf("\nthresf: %f\n",thresf); */
  /* Rprintf("head(id):\n"); */
  /* showint(id,1,10) */;
  indvec = filter_data(freqG,nGenotype,nGenotype2,thresf);
  if(verbosity==1) Rprintf("\nindvec:\n");
  if(verbosity==1) showint(indvec,1,nGenotype2);
  Rprintf("\nindvec:\n");
  showint(indvec,1,nGenotype2);
  idGtmp = filter_data_int(idG,indvec,nGenotype2);
  freqtmp2 = filter_data_double_freq(freqG,indvec,nGenotype2);
  fitntmp = filter_data_double(fitness,indvec,nGenotype2);
  nG2 = get_pos_id(idGtmp, id, nG, nGenotype2);
  Rprintf("New number of mutations in population after filtering: %d (was %d)\n",nG2,nG);
  // allocate memory
  if(!(idtmp = (int *) calloc(nG2,sizeof(int))) ||  
     !(postmp = (int *) calloc(nG2,sizeof(int)))){ 
    printf("Error in allocating memory at line 626\n\n"); 
    return NULL;
  } 
  get_new_id(idGtmp, id, idtmp, pos, postmp, nG, nGenotype2);
  /* Rprintf("hit1\n"); */
  if(verbosity==1) Rprintf("\nfreqtmp2:\n");
  if(verbosity==1) show(freqtmp2,1,nGenotype2);
  if(verbosity==1) Rprintf("\nfitntmp:\n");
  if(verbosity==1) show(fitntmp,1,nGenotype2);
  if(verbosity==1) Rprintf("idGtmp:\n");
  if(verbosity==1) showint(idGtmp,1,nGenotype2);
  if(verbosity==1) Rprintf("idtmp:\n");
  if(verbosity==1) showint(idtmp,1,nG2);
  if(verbosity==1) Rprintf("postmp:\n");
  if(verbosity==1) showint(postmp,1,nG2);
  /* Rprintf("\nfreqGtmp:\n"); */
  /* show(freqtmp2,1,nGenotype2); */
  /* Rprintf("\nfitntmp:\n"); */
  /* show(fitntmp,1,nGenotype2); */
  /* Rprintf("idGtmp:\n"); */
  /* showint(idGtmp,1,nGenotype2); */
  /* Rprintf("idtmp:\n"); */
  /* showint(idtmp,1,nG2); */
  /* Rprintf("postmp:\n"); */
  /* showint(postmp,1,nG2); */
  up = unique_pos(postmp, nG2);
  /* Rprintf("hit1.2\n"); */
  l1 = (int)nGenotype2 * (up + 1); // calculate size of vector
  Rprintf("l1 = %d\n",l1);
  if(l1<0){
    l1 = -l1;
    Rprintf("l1 negative!! up = %d,nGenotype2 = %d\n",up,nGenotype2);
  }  
 
  SEXP outX = PROTECT(allocVector(INTSXP,l1)); //= PROTECT(coerceVector(REALSXP, nx, ny));
  /* Rprintf("hit1\n"); */
  /* Rprintf("hit2\n"); */
  //PROTECT(outX = NEW_NUMERIC(l1)); 
  /* Rprintf("hit3\n"); */
  int *xab; // output
  xab = INTEGER(outX);

  /* Rprintf("hit4\n"); */
  /* unpos = make_genotype_matrix(xab,fitntmp,idGtmp,idtmp,postmp,freqtmp2,nG2,nGenotype2,l1); */
  unpos = make_genotype_matrix2(xab,idGtmp,idtmp,postmp,nG2,nGenotype2,l1);
  /* if(verbosity==1) Rprintf("\nUnique pos: %d\n",unique_pos(postmp, nG2)); */
  /* if(verbosity==1) Rprintf("nelements allocated : %d\n",l1); */
  /* if(verbosity==1) Rprintf("\nxab:\n"); */
  /* if(verbosity==1) show(xab,nGenotype2+1,(unique_pos(postmp, nG2) + 3)); */
  /* Rprintf("\nUnique pos: %d\n",unique_pos(postmp, nG2)); */
  /* Rprintf("nelements allocated : %d\n",l1); */
  /* Rprintf("\nxab:\n"); */
  /* showint(xab,nGenotype2,(up + 1)); */
  free(freqG);
  free(nrmut);
  free(fitness);
  free(id);
  free(pos);
  free(idG);

  /* free(idtmp); */
  /* free(postmp); */
  /* return output; */
  SEXP xdim = PROTECT(allocVector(INTSXP,1));
  SEXP ydim = PROTECT(allocVector(INTSXP,1));
  SEXP posnames = PROTECT(allocVector(INTSXP,up));
  SEXP freq2 = PROTECT(allocVector(REALSXP,nGenotype2));
  SEXP fitn2 = PROTECT(allocVector(REALSXP,nGenotype2));

  int *xn,*yn,*posna; // dimension
  double *freq3,*fitn3;

  xn = INTEGER(xdim);
  yn = INTEGER(ydim);
  posna = INTEGER(posnames);
  fitn3 = REAL(fitn2);
  freq3 = REAL(freq2);


  xn[0] = nGenotype2;
  yn[0] = up + 3;

  for(i = 0; i < up; i++)
    posna[i] = unpos[i];

  for(i = 0; i < nGenotype2; i++){
    freq3[i] = (double)freqtmp2[i];
    fitn3[i] =  (double)fitntmp[i];
  }

  if(flagg==0){
    SEXP vec = PROTECT((allocVector(VECSXP,6)));
    SET_VECTOR_ELT(vec,0,outX);
    SET_VECTOR_ELT(vec,1,posnames);
    SET_VECTOR_ELT(vec,2,freq2);
    SET_VECTOR_ELT(vec,3,fitn2);
    SET_VECTOR_ELT(vec,4,xdim);
    SET_VECTOR_ELT(vec,5,ydim);

    // free pointers
    /* free(posna); */
    /* free(unpos); */
    /* free(freqtmp2); */
    /* free(fitntmp); */
    /* free(freq3); */
    /* free(fitn3); */
    /* free(fitn); */
    /* free(gdiv); */
    /* free(divers); */
    /* free(xn); */
    /* free(yn); */
    /* Rprintf("hit5.1\n"); */
  // UNPROTECT(2);
    UNPROTECT(6);
 /* Rprintf("hit9:\n"); */
    return(vec);
  }
  else{
    SEXP fitnv = PROTECT(allocVector(REALSXP,(nGen+1)));
    SEXP gdv = PROTECT(allocVector(REALSXP,(nGen+1)));
    SEXP divv = PROTECT(allocVector(REALSXP,(nGen+1)));

    SEXP vec = PROTECT((allocVector(VECSXP,9)));

    double *fv,*gv,*dv;

    fv = REAL(fitnv);
    gv = REAL(gdv);
    dv = REAL(divv);

    // copy vectors
    for(i = 0; i < (nGen+1); i++){
      fv[i] = fitn[i];
      gv[i] = gdiv[i];
      dv[i] = divers[i];
    }

    SET_VECTOR_ELT(vec,0,outX);
    SET_VECTOR_ELT(vec,1,posnames);
    SET_VECTOR_ELT(vec,2,freq2);
    SET_VECTOR_ELT(vec,3,fitn2);
    SET_VECTOR_ELT(vec,4,fitnv);
    SET_VECTOR_ELT(vec,5,gdv);
    SET_VECTOR_ELT(vec,6,divv);
    SET_VECTOR_ELT(vec,7,xdim);
    SET_VECTOR_ELT(vec,8,ydim);


    // free pointers
    /* free(posna); */
    /* free(unpos); */
    /* free(freqtmp2); */
    /* free(fitntmp); */
    /* free(fv); */
    /* free(gv); */
    /* free(dv); */
    /* free(fitn); */
    /* free(gdiv); */
    /* free(divers); */
    /* free(freq3); */
    /* free(fitn3); */
    /* free(xn); */
    /* free(yn); */

    UNPROTECT(9);
    return vec;
  }


}
/* void vect_cpy2(double *x, SEXP y, int n){ */

/*   int i; */

/*   for(i = 0; i < n; i++) */
/*     y[i] = x[i]; */


/* } */
void get_new_id(int *idG, int *id, int *idnew, int *pos,int *posnew, int nG, int nGenotype){
  // call function after allocating memory of new id and pos vectors
  int i,j,k = 0;

  for(i = 0; i < nGenotype; i++)
    for(j = 0; j < nG; j++)
      if(idG[i] == id[j]){ 
	idnew[k] = id[j];
	posnew[k] = pos[j];
	/* Rprintf("i = %d, j = %d, k = %d, idnew[%d] = %d, posnew[%d] = %d\n",i,j,k,k,idnew[k],k,posnew[k]); */
	k++;
      }

}
int get_pos_id(int *idG, int *id, int nG, int nGenotype){
  // get size of new id and pos vectors
  int i,j,k = 0;

  for(i = 0; i < nGenotype; i++)
    for(j = 0; j < nG; j++)
      if(idG[i] == id[j]) 
	k++;

  return k;

}
int *filter_data_int(int *idG,int *index,int k){

  int i,j,*new;

  new = (int *) R_alloc(k,sizeof(int));


  for(i = 0; i < k; i++){
    new[i] = idG[index[i]];
    /* Rprintf("k = %d, i = %d, index[%d] = %d, new[%d] = %d, idG[%d] = %d\n",k,i,i,index[i],i,new[i],index[i],idG[index[i]]); */
  }

  return new;

}
double *filter_data_double(double *freqG,int *index,int k){

  int i,j;
  double *new;

  new = (double *) R_alloc(k,sizeof(double));

  for(i = 0; i < k; i++)
    new[i] = freqG[index[i]];

  return new;

}
double *filter_data_double_freq(double *freqG,int *index,int k){

  int i,j;
  double *new,sum = 0.;

  new = (double *) R_alloc(k,sizeof(double));

  for(i = 0; i < k; i++){
    new[i] = freqG[index[i]];
    sum += new[i];
  }

  for(i = 0; i < k; i++)
    new[i]=new[i]/sum;
    
  return new;

}
int *filter_data(double *freqG,int nGenotype,int k,double thresf){
  // determine which genotypes that passes the threshold set
  int i,j=0,*index;

  index = (int *) R_alloc(k,sizeof(int));

  for(i = 0; i < nGenotype; i++)
    if(freqG[i] >= thresf){
      index[j] = i;
      /* Rprintf("i = %d, j = %d, index[%d] = %d\n",i,j,j,index[j]); */
      j++;
    }

  return index;

}

int upd_ngen(double *freqG, int nGenotype, double thresf){ 
  // determine how many genotypes that passes the threshold set
  int i,k = 0;

  for(i = 0; i < nGenotype; i++)
    if(freqG[i] >= thresf)
      k++;

  return k;

}
int unique_pos(int *pos, int nG){
  
  int i,j,k=1,flag;
  int *tmp;

  tmp = (int *) calloc(nG,sizeof(int));
  tmp[0] = pos[0];
  for(i = 1; i < nG; i++){
    flag = 0;
    for(j = 0; j < i; j++){
      if(tmp[j] == pos[i])
	flag = 1;
    }
    if(flag==0){
      k++;
      tmp[k] = pos[i];
    }
  }

  free(tmp);

  return k;

}
int *get_pos(int *pos, int nG,int nGG){
  
  int i,j,k=0,flag;
  int *tmp;
  //Rprintf("entering function get_pos...\n");
  tmp = (int *) calloc(nGG,sizeof(int));
 /* Rprintf("hit1\n"); */
  tmp[0] = pos[0];
 /* Rprintf("hit2\n"); */
  for(i = 1; i < nG; i++){
    /* Rprintf("hit%d\n",2+i); */
    flag = 0;
    for(j = 0; j < i; j++){
      if(tmp[j] == pos[i])
	flag = 1;
    }
    if(flag==0){
      k++;
      tmp[k] = pos[i];
    }
  }
  //Rprintf("leaving function get_pos\n");
  return tmp;

}
int *make_genotype_matrix(double *xab, double *fitness,int *idG,int *id,int *pos,double *freqG,int nG,int nGenotype, int size){
  // set up clone id x marker matrix
  double *x;
  int flag,i=0,j,l,nGG=0,k,tmp,tmp2,nGenotype1 = nGenotype + 1,*unique;
  Rprintf("entering function make_genotype_matrix...\n");
  // get size of matrix
  nGG = unique_pos(pos, nG);
  Rprintf("unique positions: %d\nnG: %d\nnGenotype: %d\n",nGG,nG,nGenotype);
  // size = (nGenotype+1) * (nGG + 3); // rows x columns
  Rprintf("size: %d\n",size);

  if(idG[0]==0) k = 1; // determine starting row
  else k = 0;

  unique = get_pos(pos, nG,nGG);
  /* Rprintf("unique position vector\n"); */
  /* showint(unique,1,nGG); */
  for (j = 0; j < size; j++) xab[j] = 0.; // reset data 
  /* xab[0] = 0.; */
  // for (j = 0; j < nGG; j++){ 
  //  xab[(j+1)*nGenotype1] = (double)unique[j]; // set the first row to positions of mutations
    /* Rprintf("0. x[%d,%d] = %f\n",0,j+1,xab[(j+1)*nGenotype1]); */
  //}
  /* xab[(nGG+1)*nGenotype1] = 0.; */
  /* xab[(nGG+2)*nGenotype1] = 0.; */
 
  /* Rprintf("hit2\n"); */
  // first entry manual
  if(id[0] != 0) 
    xab[nGenotype+k] = 1.;
  /* Rprintf("0.5. x[%d,%d] = %f\n",k,1,xab[nGenotype1+k]); */
  for (i = 1; i < nG; i++){ // continue with the rest of the mutations
    if(id[i]!=id[i-1]){ // move to new row in genotype x marker matrix
      k++;
      /* Rprintf("Moving from %d to row %d\n",k-1,k); */
    }
    for(j = 0; j < nGG; j++){ // find mutation among unique set
      /* Rprintf("unique[%d] = %d, pos[%d] = %d\n",j,unique[j],i,pos[i]); */
      if(unique[j] == pos[i]){
	xab[(j+1)*nGenotype+k] = 1.;
	/* Rprintf("1. i = %d, j = %d, x[%d,%d] = %f\n",i,j,k,j+1,xab[(j+1)*nGenotype1+k]); */
      }
    }
      /* else */
      /* 	xab[(j+1)*nGenotype1+k] = 0.; */
  }
  /* Rprintf("hit3\n"); */
  for(i = 0; i < nGenotype; i++){
    xab[i] = idG[i];
    xab[(nGG+1)*nGenotype+i] = freqG[i];
    xab[(nGG+2)*nGenotype+i] = fitness[i];
  }
  /* Rprintf("hit4\n"); */
  return unique;

}
int *make_genotype_matrix2(int *xab,int *idG,int *id,int *pos,int nG,int nGenotype, int size){
  // set up clone id x marker matrix w.o. frequency or fitness info
  double *x;
  int flag,i=0,j,l,nGG=0,k,tmp,tmp2,nGenotype1 = nGenotype + 1,*unique;
  Rprintf("entering function make_genotype_matrix...\n");
  // get size of matrix
  nGG = unique_pos(pos, nG);
  Rprintf("unique positions: %d\nnG: %d\nnGenotype: %d\n",nGG,nG,nGenotype);
  // size = (nGenotype+1) * (nGG + 3); // rows x columns
  Rprintf("size: %d\n",size);

  if(idG[0]==0) k = 1; // determine starting row
  else k = 0;

  unique = get_pos(pos, nG,nGG);
  /* Rprintf("unique position vector\n"); */
  /* showint(unique,1,nGG); */
  for (j = 0; j < size; j++) xab[j] = 0.; // reset data 
  /* xab[0] = 0.; */
  // for (j = 0; j < nGG; j++){ 
  //  xab[(j+1)*nGenotype1] = (double)unique[j]; // set the first row to positions of mutations
    /* Rprintf("0. x[%d,%d] = %f\n",0,j+1,xab[(j+1)*nGenotype1]); */
  //}
  /* xab[(nGG+1)*nGenotype1] = 0.; */
  /* xab[(nGG+2)*nGenotype1] = 0.; */
 
  /* Rprintf("hit2\n"); */
  // first entry manual
  if(id[0] != 0) 
    xab[nGenotype+k] = 1.;
  /* Rprintf("0.5. x[%d,%d] = %f\n",k,1,xab[nGenotype1+k]); */
  for (i = 1; i < nG; i++){ // continue with the rest of the mutations
    if(id[i]!=id[i-1]){ // move to new row in genotype x marker matrix
      k++;
      /* Rprintf("Moving from %d to row %d\n",k-1,k); */
    }
    for(j = 0; j < nGG; j++){ // find mutation among unique set
      /* Rprintf("unique[%d] = %d, pos[%d] = %d\n",j,unique[j],i,pos[i]); */
      if(unique[j] == pos[i]){
	xab[(j+1)*nGenotype+k] = 1.;
	/* Rprintf("1. i = %d, j = %d, x[%d,%d] = %f\n",i,j,k,j+1,xab[(j+1)*nGenotype1+k]); */
      }
    }
      /* else */
      /* 	xab[(j+1)*nGenotype1+k] = 0.; */
  }
  /* Rprintf("hit3\n"); */
  for(i = 0; i < nGenotype; i++){
    xab[i] = idG[i];
    /* xab[(nGG+1)*nGenotype+i] = freqG[i]; */
    /* xab[(nGG+2)*nGenotype+i] = fitness[i]; */
  }
  /* Rprintf("hit4\n"); */
  return unique;

}
void showint(int *M,int r, int c){
 
 int i,j;

 for(i = 0; i < r; i++){
   for(j = 0; j < c; j++)
     printf("%d ",M[j*r+i]);
   printf("\n");
 }


}
void show(double *M,int r, int c){
 
 int i,j;

 for(i = 0; i < r; i++){
   for(j = 0; j < c; j++)
     printf("%.2f ",M[j*r+i]);
   printf("\n");
 }
 //%.2f

}
void freq(int *x, int r, int c){

 int i,j;
 int *lc;
 double freq;

 lc = (int *) calloc(c,sizeof(int));

 printf("Allele frequency at each locus:\n");
 for(i = 0; i < r; i++)
   for(j = 0; j < c; j++)
     lc[j] += x[j*r+i];
   
   

 /* for(j = 0; j < c; j++) printf("%d\t",j+1); */
 /* printf("\n"); */

 for(i = 0; i < c; i++){
   freq = (double)lc[i]/r;
   printf("%1.7f\t",freq);
 }
 printf("\n");

 free(lc);

}
double draweff(double p1, double p2, double mu){ //, const gsl_rng *r
  // draw effect of a mutation on fitness
  // p1 p2
  // mu is the mean of the exponential distribution or here:
  // the average selective advantage (or disadvantage) of mutations
  // in Rozen et al. 2002, mu = 1/35 was used
 
  double eff=0.,cdf[3],u;
  /* step 0: create cdf */
  cdf[0] = p1;
  cdf[1] = p2+p1;
  cdf[2] = 1.;
  /* step 1: decide which distribution to use */
  //runif(double a, double b) 
  u = runif(0.0,1.0);//gsl_rng_uniform (r);
 
  /* step 2: make a draw from that distribution */
  if(u <= cdf[0])
    //rexp(double sl)
    eff = - rexp(mu); //gsl_ran_exponential (r, mu);
  else if(u > cdf[1])
    eff = rexp(mu); //gsl_ran_exponential (r, mu);
  
  /* printf("Mutation effect on fitness: %f\n",eff); */
  
  return eff;
  
}
void fix_mut(int *x, double *mut,int r, int c, int nloci){

 int i,j;
 int *lc;
 double freq;

 lc = (int *) calloc(c,sizeof(int));


 for(i = 0; i < r; i++)
   for(j = 0; j < c; j++)
     lc[j] += x[j*r+i];
   
   

 /* for(j = 0; j < c; j++) printf("%d\t",j+1); */
 /* printf("\n"); */
 printf("Effect of fixed mutations:\n");
 for(i = 0; i < c; i++){
   freq = (double)lc[i]/r;
   if(freq==0. && i > nloci) printf("%f ",mut[i-nloci]);
   //printf("%1.7f\t",freq);
 }
 printf("\n");

 free(lc);

}
void freq_stat(double *freqv, int nlocitot, int nGen1, int nloci,double *muteff){

  int i,j,k=0,hit,hit2,nGen = nGen1-1,l;


  printf("Mutation\tintroduced\tfixed\teffect\n");
  for(i = nloci; i < nlocitot; i++){
    l = i-nloci;
    if(freqv[nGen*nlocitot+i] <= 0.00001){
      k++;
      hit = 0;
      hit2 = 0;
      printf("%d\t\t",k);
      for(j = 0; j < nGen1; j++){

	if(hit==0 && freqv[j*nlocitot+i] < 1.){
	  printf("%d\t\t",j);
	  hit = 1;
	}
	if(hit2==0 && freqv[j*nlocitot+i] <= 0.00001){
	  printf("%d\t",j);
	  hit2 = 1;
	}
      }
      printf("%f\n",muteff[l]);
    }
  }
  

}
int det_sp(int *x,int popsize,int nloci){

  int i,j,p=0;

  for(i = 0; i < popsize; i++)
    for(j = 0; i < nloci; j++)
      if(x[j*popsize+i]) p++;

  return p;

}
void sparse(int *x,int *pos,int *id,int popsize,int nloci){

  int i,j,p=0;

  for(i = 0; i < popsize; i++)
    for(j = 0; j < nloci; j++)
      if(x[j*popsize+i]==0){
	id[p] = i;
	pos[p] = j;
	p++;
      }



}
int count_mut_gen(int *x,int nGenotype,int nloci){

  int i,j,nG=0;

  for(i = 0; i < nGenotype; i++)
    for(j = 0; j < nloci; j++)
      if(x[j*nGenotype+i]==1) nG++;

  return nG;
}
void sparse_gen(int *x,int *pos,int *id,double *freqG,int *idG,int nG, int nGenotype,int nloci,int size){ //,const gsl_rng * r

  int i,j,k = 0;
  double u,tmp;
  unsigned long int n;
  //n = (unsigned long int) size;
  
  for(i = 0; i < nGenotype; i++)
    for(j = 0; j < nloci; j++)
      if(x[j*nGenotype+i]==0){
        u = runif(0.0,1.0);
        tmp = u*size;
	id[k] = idG[i];
	pos[k] = (int) tmp;//gsl_rng_uniform_int (r, n);
	Rprintf("u = %f, tmp = %f, pos[%d] = %d\n",u,tmp,k,pos[k]);
	k++;
      }


}
void sparse_gen_pos(int *x,int *pos,int *id,double *freqG,int *idG,int nG, int nGenotype,int nloci,double *posi){ //,const gsl_rng * r

  int i,j,k = 0;
  //unsigned long int n;
  //n = (unsigned long int) size;
  
  for(i = 0; i < nGenotype; i++)
    for(j = 0; j < nloci; j++)
      if(x[j*nGenotype+i]==1){
        //u = runif(0.0,1.0);
	//  tmp = u*size;
	id[k] = idG[i];
        if(i > 0)
	  pos[k] = (int)posi[i-1];//gsl_rng_uniform_int (r, n);
	else
	  pos[k] = (int)posi[i];
	//Rprintf("i = %d, j = %d, pos[%d] = %d\n",i,j,k,pos[k]);
	k++;
      }


}
void freq_sp(int *pos, int n1,int popsize,int nloci){

  int i,j,n;
  double f;

  printf("Allele frequency:\n");
  for(i = 0; i < nloci; i++){
    n = 0;
    for(j = 0; j < n1; j++)
      if(pos[j]==i)
	n++;
    f = (double)n/popsize;
    printf("%1.4f ",1.-f);
  }
  printf("\n");
}
double walltime( double *t0 )
{

  double mic, time;
  double mega = 0.000001;
  struct timeval tp;
  struct timezone tzp;
  static long base_sec = 0;
  static long base_usec = 0;

  (void) gettimeofday(&tp,&tzp);
  if (base_sec == 0)
    {
      base_sec = tp.tv_sec;
      base_usec = tp.tv_usec;
    }

  time = (double) (tp.tv_sec - base_sec);
  mic = (double) (tp.tv_usec - base_usec);
  time = (time + mic * mega) - *t0;
  return(time);

}
void update_mut_freq(int *pos,int *id,double *freqG,int *idG,int nG,int nGenotype, int nloci,int gen,FILE *fp){

  // nloci is the number of loci simulated
  // nG is the size of pos and id vectors
  
  int i,j,k;
  double tmp;

  /* printf("Mutation frequency:\t"); */

  
  for(i = 0; i < nloci; i++){
    tmp = 1.;
    for(j = 0; j < nG; j++)
      if(pos[j]==i){
	for(k = 0; k < nGenotype; k++)
	  if(id[j] == idG[k])
	    tmp -= freqG[k]; 
      }
    fprintf(fp,"%d\t%d\t%f\n",gen,i,tmp);
  }


}
double *calc_gen_prob(double *freqG,double *fitness,int nGenotype,int popsize){

  int i;
  double sum = 0., *dG;
  
  dG = (double *) calloc(popsize,sizeof(double));
  
  for(i = 0; i < nGenotype; i++){
    if(fitness[i] > 0.)
      sum += fitness[i]*freqG[i]*popsize;
    else{
      Rprintf("\n\n************************\nNB! WARNING! freqG[%d] = %f, fitness[%d] = %f\n\n****************************\n",i,freqG[i],i,fitness[i]);
      fitness[i] = 0.;
    }
  }
  for(i = 0; i < nGenotype; i++){ 
    dG[i] = freqG[i]*popsize*fitness[i]/sum;
    /* if(dG[i] < 0.){ */
      
    /*   dG[i] = 0.;  */
    /*   fitness[i] = 0.1; */
    /* } */
  }
  return dG;

}
double mean_fitness(double *fitness, double *freqG, int nGenotype){

  int i;
  double tmp=0.;

  for(i = 0; i < nGenotype; i++) tmp += fitness[i]*freqG[i];

  return tmp;

}
double nucl_div(double *freqG, int *idG, int *id, int *pos, int nGenotype,int nG, double thres){
  // computes average pair-wise nucleotide diversity

  int i,j,pi_ij,cnt1=0,cnt2=0;
  double pi = 0.,sum=0.;

  if(nGenotype > 1){
    for(i = 1; i < nGenotype; i++) 
      for(j = 0; j < i; j++){
	if(freqG[i] < thres || freqG[j] < thres){
	  pi_ij = 0.;
	  cnt1++;
	}
	else{
	  pi_ij = differ(id,pos,idG,i,j,nG);
	  cnt2++;
	}
	sum += freqG[i]*freqG[j]*pi_ij;

      }
    pi = (double)(nGenotype/(nGenotype-1.))*sum;
    if(pi > 4){ 
      printf("Warning!!!, unreasonably high diversity\n");
      /* show(freqG) */
    }
  }
  /* printf("number of times with threshold activated: %d\nnumber of times without activated threshold: %d\n",cnt1,cnt2); */
  return pi;

}
int differ(int *id,int *pos,int *idG,int i1,int i2,int nG){
  // count sequence difference per site between seq i1 and i2
  
  int i,j,n=0,n1=0,n2=0,n3=0;

  for(i = 0; i < nG; i++){ 
    if(id[i]==idG[i1]) n1++;  // count nr of mutations in i1 
    if(id[i]==idG[i2]) n2++; // count nr of mutations in i2
    for(j = 0; j < nG; j++){
      if((id[i]==idG[i1]) && (id[j]==idG[i2]) && (pos[i] == pos[j])) n3++; // count common mutations
    }
  }
       
  n = n1 + n2 - 2*n3;
  /* if(n > 2){ */
  /*   printf("%d vs %d: n1 = %d, n2 = %d, n3 = %d, n = %d\n",idG[i1],idG[i2],n1,n2,n3,n); */
  /*   printf("id:\n"); */
  /*   showint(id,1,nG); */
  /*   printf("pos:\n"); */
  /*   showint(pos,1,nG); */
  /*   printf("idG:\n"); */
  /*   printf("i1 = %d, i2 = %d\n",i1,i2); */
  /* } */
  return n;

}

int upd_geno_nG(unsigned int *numG,int *id,int *idG,int nGenotype,int nG){

  int j,k,nG2 = nG;

  for(j = 0; j < nGenotype; j++)
    if(numG[j] == 0){ 
      for(k = 0; k < nG; k++)
	if(id[k]==idG[j]) nG2--;// NB NB update!!
    }
 
  return nG2;

}
int upd_geno_nGen(unsigned int *numG,int nGenotype){

  int j,nGenotype2 = nGenotype;

  for(j = 0; j < nGenotype; j++)
    if(numG[j] == 0) nGenotype2--;
 
  return nGenotype2;

}
void read_gen_2(double *inX,int *x, double *freqG, double *fitness, int nGenotype, int nloci){

  int i,j,nl1 = nloci+1,nl2 = nloci+2;
  
  for(i = 0; i < nGenotype; i++){ 
    for(j = 0; j < nloci; j++)
      x[j*nGenotype+i] = (int) inX[(j+1)*nGenotype+i];
    freqG[i] = inX[nl1*nGenotype+i];
    fitness[i] = inX[nl2*nGenotype+i];
  }
}
double *upd_geno(double *inX,int *x,double *freqtmp,double *fitness,int *idG,int nGenotype,int nloci,int nGenotype2){
  
  int i,j,k=0,nl1 = nloci+1,nl2 = nloci+2;
  double *newfreq,*freqtmp2,sum = 0.;

  newfreq = (double *) calloc(nGenotype2,sizeof(double));

  for(i = 0; i < nGenotype; i++){ 
    if(freqtmp[i] > 0.){
      idG[k] = (int) inX[i];
      for(j = 0; j < nloci; j++)
	x[j*nGenotype2+k] = (int) inX[(j+1)*nGenotype+i];
      newfreq[k] = freqtmp[i];//inX[nl1*nGenotype+i];
      fitness[k] = inX[nl2*nGenotype+i];
      k++; // id counter
    }
  }
  for(i = 0; i < k; i++) sum += newfreq[i]; 
  for(i = 0; i < k; i++) newfreq[i] /= sum;
  return newfreq;

}
double *read_gen_3(double *inX,double *freqG, int nGenotype, int nloci,int popsize){
  
  int i,j,nl1 = nloci+1,nl2 = nloci+2,sum=0;
  double *newfreq,tmpfit;
  int *tmp,tmpx;

  newfreq = (double *) calloc(nGenotype,sizeof(double));
  tmp = (int *) R_alloc(nGenotype,sizeof(int));

  for(i = 0; i < nGenotype; i++){ 
    for(j = 0; j < nloci; j++)
      tmpx = (int) inX[(j+1)*nGenotype+i];
    freqG[i] = inX[nl1*nGenotype+i];
    tmpfit = inX[nl2*nGenotype+i];
  }
  /* rmultinom(nGenotype,freqG,popsize,tmp); */
  Rprintf("nGenotype = %d, popsize = %d:\n",nGenotype,popsize);
  //rmultinom(nGenotype,freqG,popsize,tmp);
  rmultinom(popsize,freqG,nGenotype,tmp); //gsl_ran_multinomial(gBaseRand, nGenotype, popsize, freqG, tmp);
  Rprintf("original frequency:\n");
  show(freqG,1,nGenotype);
  /* Rprintf("tmp:\n"); */
  /* showint(tmp,1,nGenotype); */
  for(i = 0; i < nGenotype; i++) sum += tmp[i];
  if(sum > 0.)
    for(i = 0; i < nGenotype; i++) newfreq[i] = (double)tmp[i]/sum;
  else{
    Rprintf("Error! No genotypes were selected for cultivation. Please look at your imput data\n");
    exit(-1);
  }
  Rprintf("new frequency after draws from multinomial distribution:\n");
  show(newfreq,1,nGenotype);

  /* Rprintf("tmp:\n"); */
  /* showint(tmp,1,nGenotype); */
  /* free(tmp); */

  return newfreq;

}
int cnt_nonz(double *freqG,int nGenotype){

  int i,cnt=0;

  for(i = 0; i < nGenotype; i++) 
    if(freqG[i] > 0.) cnt++;

  return cnt;

}
int reduce_loci_1(int *x,int nGenotype,int nloci){

  int nloci2=0,i,j,tmp;

  for(i = 0; i < nloci; i++){
    tmp = 0;
    for(j = 0; j < nGenotype; j++) 
      if(x[i*nGenotype+j]==0) tmp = 1;
    if(tmp==1) nloci2++;

  }

  printf("Number of segregating loci: %d\n",nloci2);

  return nloci2;

}
int *reduce_loci_2(int *x,int nGenotype, int nloci, int nloci2){

  int i,j,*x2,tmp,k=0;

  x2 = (int *) malloc(nGenotype*nloci2*sizeof(int));

  for(i = 0; i < nloci; i++){
    tmp = 0;
    for(j = 0; j < nGenotype; j++) 
      if(x[i*nGenotype+j]==0) tmp = 1;
    if(tmp==1){ // copy x to x2
      for(j = 0; j < nGenotype; j++) 
	x2[k*nGenotype+j] = x[i*nGenotype+j];
      k++;
    }
  }

  return x2;

}
double *get_freq(int *numGtmp,int nGenotype){

  int i,sum=0;
  double *tmpfreq;

  tmpfreq = (double *) calloc(nGenotype,sizeof(double));

  for(i = 0; i < nGenotype; i++) sum += numGtmp[i];

  for(i = 0; i < nGenotype; i++) tmpfreq[i] = (double)numGtmp[i]/sum;

  return tmpfreq;

}
double gene_div(double *freqG, int *pos, int *idG, int *id, int nGenotype, int nG,int popsize){
  // only for polymorphic loci

  double gd=0.;
  // save allele frequency on format: generation, locus, allele variant, frequency 
  int i,j,k,cnt1 = 0;
  double tmp1,freq,sum,sl = 0.,ps1;

  int *tmp,nloci=1;
  double *cnt;

  ps1 = (double)popsize/(popsize-1);
  tmp = (int *) R_alloc(nG,sizeof(int));

  // find number of polymorphic loci 
  tmp[0] = pos[0];
  /* Rprintf("hit3\n"); */
  for(i = 1; i < nG; i++){
    /* Rprintf("nG = %d, i = %d, nloci = %d\n",nG,i,nloci); */
    cnt1 = 0;
    for(k = 0; k < nloci; k++){
      if(tmp[k]==pos[i]) cnt1 = 1;
    }
    if(cnt1==0){ // new loci
      tmp[nloci]=pos[i];
      nloci++;
    }
  }
  /* Rprintf("total number of polymorphic loci: %d\n",nloci); */
  cnt = (double *) calloc(nloci,sizeof(double));
  //  for(i = 0; i < nloci; i++){
  for(j = 0; j < nG; j++){
    for(k = 0; k < nGenotype; k++){
      if(id[j]==idG[k]){
	for(i = 0; i < nloci; i++){
	  if(tmp[i]==pos[j])
	    cnt[i] += freqG[k];
	}
      }
    }
  }
  /* printf("allele frequencies:\n"); */
  /* show(cnt,1,nloci); */
  /* printf("allele pos:\n"); */
  /* showint(tmp,1,nloci); */
  /* free(tmp); */
  for(i = 0; i < nloci; i++){
    sum = 0.;
    sum += cnt[i]*cnt[i];
    tmp1 = 1-cnt[i];
    sum += tmp1*tmp1;
    sl += (double)(1.-sum)*ps1;
  }
  free(cnt);

  gd = (double)sl/nloci;

  /* printf("gd = %f\tsl = %f\tsum = %f\tnloci = %d\n",gd,sl,sum,nloci); */

  return gd;

}
int find_max(double *inX, int nGenotype){
  // find maximum clone ID number
  int n = 0,i;

  for(i = 0; i < nGenotype; i++){
    if(inX[i] > n)
      n = (int)inX[i];

  }
  return n;

}
int maxpos(int *pos, int nG){

  int i,m = -2;

  for(i = 0; i < nG; i++)
    if(pos[i] > m)
      m = pos[i];

  return m; 

}

int *vect_cpy_int(int *x, int n){

  int i,*res;

  res = (int *) calloc(n,sizeof(int));

  for(i = 0; i < n; i++)
    res[i] = x[i];

  return res;

}
double *vect_cpy(double *x, int n){
 
  int i;
  double *res;

  res = (double *) calloc(n,sizeof(double));

  for(i = 0; i < n; i++)
    res[i] = x[i];

  return res;

}
double *vect_const_mult(int *x, double a, int n){

  int i;
  double *res;

  res = (double *) calloc(n,sizeof(double));

  for(i = 0; i < n; i++)
    res[i] = (double)a*x[i];

  return res;

}
