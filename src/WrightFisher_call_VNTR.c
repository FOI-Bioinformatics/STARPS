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

/**************************************************************/
/* simulation program for generating populations              */
/* according to a Wright-Fisher model                         */
/* version 7: efficient memory storage model (only genotypes) */
/*                                                            */
/* By Petter Lindgren and Jon Ahlinder                        */
/*                                                            */
/* input: ./WrightFisher nGen popsize expfactor mutrate nloci */
/* nGen is number of generations                              */
/* popsize is initial population size                         */
/* expfactor is the expanding factor. If expfactor=1,         */
/* then contant population size is maintained.                */
/* If expfactor=2 population size dubbles each gen            */
/* mutrate is the mutation rate, per genome per rep           */
/* nloci is the starting number of segregating loci in pop    */
/* output                                                     */
/**************************************************************/

/* compile: gcc -o WrightFisher_v3 WrightFisher_v3.c -Wall -I/usr/include -lm -lgsl -lgslcblas */
/* in command prombt: R CMD SHLIB -lgsl -lgslcblas -lm WrightFisher_int.c */

/* function declaration */
/* double *concat(double *freqG,double *fitness,int nGenotype2,int nlocitot,int *pos,int *id,int *idG,int nG); */
/* void print_genetic(int nGen,int *idG,int *pos,int *id,double *fitness,FILE *fp4); */
int find_max2(double *inX, int nGenotype);
int find_no_of_abundant_clones(double *freqG, int nGenotype, double thresf);
double find_normalizing_constant(double *freqG, int nGenotype, double thresf);
double gene_div_VNTR(double *freqG, int *x, int *nall, int *allele, int nGenotype, int nloci,int popsize);
double nucl_div_VNTR(double *freqG, int *x, int nGenotype, int nloci, double thres);
int comp_dim(int *x2,int *xnew, int nloci,int nGenotype2,int nmut);
int concat_matr(int *x2,int *xnew, int *x,int *numG,int *numGtmp,double *fitness,double *fitntmp, int *idG,int *idGtmp,int nloci,int nGenotype2,int nmut,int dim,int maxnr);
void upd_param(int *numGtmp,int *idGtmp,double *fitntmp,unsigned int *numG,int *idG,double *fitness,int nGenotype, int nGenotype2);
int *upd_geno_matr(int *x,int *numG,int nloci,int nGenotype, int nGenotype2);
int *upd_no_alleles(int *numG,int *x,int nGenotype,int nloci);
int find_max(int *nall,int nloci);
int sum_VNTR_alleles(int *nall, int nloci);
int *get_VNTR_alleles(int *x, int nallsum, int nloci, int nGenotype);
double *get_freq(int *numGtmp,int nGenotype);
int *reduce_loci_2(int *x,int nGenotype, int nloci, int nloci2);
int reduce_loci_1(int *x,int nGenotype,int nloci);
double *upd_geno(double *inX,int *x,double *freqtmp,double *fitness,int *idG,int nGenotype,int nloci,int nGenotype2);
int cnt_nonz(double *freqG,int nGenotype);
void read_gen(FILE *fps,int *x, double *freqG, double *fitness, int nGenotype, int nloci);
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
void fix_mut(int *x, double *mut,int r, int c, int nloci);
void freq_stat(double *freqv, int nlocitot, int nGen1, int nloci,double *muteff);
int det_sp(int *x,int popsize,int nloci);
void sparse(int *x,int *pos,int *id,int popsize,int nloci);
double walltime(double *t0);
/* int count_genotypes(int *x,int popsize,int nloci); */
/* int count_mutations(int *x,int popsize,int nloci); */
int count_mut_gen(int *x,int nGenotypes,int nloci);
/* void sparse_genotypes(int *x,int *pos,int *id,double *freqG,int nG, int nGenotype,int popsize,int nloci); */
void sparse_gen(int *x,int *pos,int *id,double *freqG,int nG, int nGenotype,int nloci);
void update_mut_freq(int *pos,int *id,double *freqG,int *idG,int nG,int nGenotype, int nloci,int gen,FILE *fp);
void update_mut_freq_VNTR(int *nall,int *x, int *allele,double *freqG,int nGenotype, int nloci,int gen,FILE *fp);
double *calc_gen_prob(double *freqG,double *fitness,int nGenotype,int popsize);
double mean_fitness(double *fitness,double *freqG, int nGenotype);
int *num_all(int *x,int nloci,int nGenotype);
/* void WF(int *nGenp, int *popsizep, int *expfactorp, double *mutratep, int *sizep, int *nlocip, int *verbosityp, int *nGenotypep, double *p1p, double *mup,char **filenamep,double *fitn,double *divers); */

/* int main(int argc, char *argv[]) */
/* { */
/*   int *nGenp,*popsizep,*expfactorp,*sizep,*nlocip,*verbosityp,*nGenotypep; */
/*   double *mutratep, *p1p,*mup,*fitn,*divers; */
/*   char *filenamep[1]; */

 
/*   printf(); */
/*   *nGenp = 3; */
/*   *popsizep = 100; */
/*   *expfactorp = 1; */
/*   *sizep = 1700000; */
/*   *nlocip = 3; */
/*   *verbosityp = 1; */
/*   *nGenotypep = 4; */
/*   *mutratep = 1.0E-8; */
/*   *p1p = 0.9; */
/*   *mup = 0.1; */
/*   filenamep[0] = (char *) malloc(200*sizeof(char)); */
/*   fitn = (double *) malloc(3*sizeof(double)); */
/*   divers = (double *) malloc(3*sizeof(double));  */

/*   filenamep[0] = "input_pop.txt";  */

/*   printf("Call WF function...\n"); */

/*   WF(nGenp, popsizep, expfactorp, mutratep, sizep, nlocip, verbosityp, nGenotypep, p1p, mup,filenamep,fitn,divers); */

/*   printf("Done!\n"); */

/*   free(fitn); */
/*   free(divers); */

/*   return 1; */

/* } */
//void WF_VNTR(int *nGenp, int *popsizep, int *expfactorp, double *mutratep, int *nlocip, int *verbosityp, int *nGenotypep,double *fitn,double *divers,double *gdiv,double *inX,char **filen1,char **filen3, int *rSeedp, double *thresp, double *thresfrp,int *saveadd)
//{
SEXP WF_VNTR(SEXP nGenp, SEXP popsizep, SEXP expfactorp, SEXP mutratep, SEXP nlocip, SEXP verbosityp, SEXP nGenotypep, SEXP flagp, SEXP inXp, SEXP posp, SEXP thresp, SEXP thresfrp){
  
  /* printf("file1: %s\tfile2: %s\tfile3: %s\tfile4: %s\n",filen1[0],filen2[0],filen3[0],filen4[0]); */
    // defining input variables
  /* int size = INTEGER_VALUE(sizep);  */ 
  int nGen=INTEGER_VALUE(nGenp), popsize=INTEGER_VALUE(popsizep), expfactor=INTEGER_VALUE(expfactorp), nloci=INTEGER_VALUE(nlocip), verbosity=INTEGER_VALUE(verbosityp), flagg=INTEGER_VALUE(flagp), nGenotype=INTEGER_VALUE(nGenotypep);
  double mutrate = NUMERIC_VALUE(mutratep), thres=NUMERIC_VALUE(thresp), thresf=NUMERIC_VALUE(thresfrp);
  /* double *fitneff = REAL(fitneffp); */
  double *inX = REAL(inXp);
  int *posi = INTEGER(posp);

  double *gdiv,*divers,*fitn;
  int ntot=0,i,j,k,nmut=0,popsize_i,popsize_i0,l1,i1,tmp1;
  int *x,*x2,*pos,*id,*start,*end,*nrmut,totpop,cnt;
  double dev,dev2,*fitness,*fitness2,u1,*freqG,*freqtmp;
  double mfitness=1.,indx2;
  unsigned int cmptmp;
  unsigned long randSeed;
  int cnt1=0,cnt2=0;
  /* gsl_vector *mutID; */
  /* specifying to use Mersenne twister MT-19937 as the uniform PRNG */
  /* gBaseRand = gsl_rng_alloc(gsl_rng_mt19937); */
  /* srand(time(NULL));      */               /* initialization for rand() */
  /* if(rSeed == 1 || rSeed <= 0) randSeed = rand();                    /\* returns a non-negative integer *\/ */
  /* else randSeed = (long unsigned int)rSeed; */
  /* printf("randSeed = %lu\n",randSeed); */
  /* gsl_rng_set (gBaseRand, randSeed);    /\* seed the PRNG *\/ */
  /* gsl_rng_env_setup(); */
  /* GetRNGstate(); */ /* generate random number state */

  /* char *file,*file2,*file3,*file4; */
  /* FILE *fp1,*fp2; */
  /* FILE *fp3; */ //,*fp4;

  /* defining time variables */
  double startTime,elapsedTime;
  double clockZero = 0.0;
 
  int nG; /* number of mutations present in the genotypes */
  unsigned int *numG,*numG2,*numGtmp;
  int *idG,*idG2, *idGtmp,*tmpID,*tmpPOS,*xnew,*x3;
  int nGenotype2,l,k1,dummy=1,tmp0,nloci2,mk,idm,j1,maxnr;
  double *fitntmp, *pG,tmpmut,mf,pi,*freqtmp2;
  double tmpdev,gd;
  int ul = 1E6; /* upper limit */

  /* new variable declarations */
  int *nall,nmax,nallsum,*allele,*allele2,*nall2,*muttmp,*mut,*nperm,*posm,*posm2;
  int mcnt = 0; /* mutation counter */
  /* int gmut = 0; */ /* number of mutated clones */
  int ntmp=0,ntot2,nl,nlint,gennr=nGenotype,km,u2=0;
  double u,tmpu,prod;
  /* double *output; */
  GetRNGstate(); // initialize random number generator
  Rprintf("posi:\n");
  showint(posi,1,nloci);

  /* error check user input */

  if(mutrate > 1.0){
    fprintf(stderr,"Unrealistic high mutation rate!\nTerminating program\n");
    return NULL;
  }
  /* set help variables */
  start = (int *) calloc(nGen,sizeof(int));
  end = (int *) calloc(nGen,sizeof(int));
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
  nrmut = (int *) calloc(nGen,sizeof(int));
  muttmp = (int *) calloc(ntot2,sizeof(int));
  nperm = (int *) calloc(nloci,sizeof(int));
  posm = (int *) calloc(ntot2,sizeof(int));
  int flag = 0;
  u = (double)nloci;
  Rprintf("nloci = %d\n",nloci);
  for(i = 0; i < ntot2; i++){
   
    /* for(k = 0; k < nloci; k++){ */
      //ntmp = gsl_ran_binomial(gBaseRand, mutrate, 1); /* generate the number of mutations */
    ntmp = rbinom(u, mutrate); // has bacteria i mutated?
    if(ntmp==1){  // if a mutation has occur
      /* muttmp[nmut]++; */ /* count number of mutations at that particular genotype */
      tmpu = runif(0.0,1.0);
      prod = tmpu*nloci;
      km = (int)prod; // generate a position from a uniform distribution 
      posm[nmut] = km; // save position
      nmut++; /* count the number of mutations */    
      tmp0 = i + start[0];
      for(j = 0; j < nGen; j++)
	if(tmp0 >= start[j] && tmp0 <= end[j])
	  nrmut[j]++;
    }
    //printf("nmut = %d, gmut = %d, i = %d\n",nmut,i);
  }
  free(start);
  free(end);

 if(nmut > 0) posm2 = (int *) R_alloc(nmut,sizeof(int));
 else posm2 = (int *) R_alloc(1,sizeof(int));
  for(i = 0; i < nmut; i++) // copy mutation vectors to save memory 
    posm2[i] = posm[i];
 
  free(posm); 
  /* gsl_vector_free(mutID); */
  if(verbosity==1) printf("Distribution of mutations over generations:\n");
  if(verbosity==1) showint(nrmut,1,nGen);
 
  if(verbosity==1) printf("\nTotal number of mutations: %d\n",nmut);
 
  /* allocate memory */
  freqG = (double *) calloc(nGenotype,sizeof(double));

  if(verbosity==1) printf("nloci = %d\tnGenotype = %d\n",nloci,nGenotype);

  freqtmp = read_gen_3(inX,freqG,nGenotype,nloci,popsize); /* read input from R */
  maxnr = find_max2(inX,nGenotype);
  maxnr++; /* increase clone id nr by one */
  free(freqG);

  /* reduce data from non-existing clones */
  nGenotype2 = cnt_nonz(freqtmp,nGenotype);
  if(verbosity==1) printf("New number of genotypes after removing non-existant variants: %d\n",nGenotype2);
  /* allocating memory */
  x = (int *) calloc(nGenotype2*nloci,sizeof(int));
  fitness = (double *) calloc(nGenotype2,sizeof(double));
  idG = (int *) calloc(nGenotype2,sizeof(int));
  freqG = upd_geno(inX,x,freqtmp,fitness,idG,nGenotype,nloci,nGenotype2);
  free(freqtmp);
  /* free(inX); */
  nGenotype = nGenotype2;
 
  /* nG = count_mut_gen(x,nGenotype,nloci); */ // Needed??? NB NB NB!!!
  /* remove non existant genotypes */

  nall = num_all(x,nloci,nGenotype); /* calculate the number of alleles per loci */
  /* nmax = find_max(nall, nloci); */
  nallsum = sum_VNTR_alleles(nall, nloci);
  allele = get_VNTR_alleles(x, nallsum, nloci, nGenotype);

  if(verbosity==1) printf("nG = %d\tnGenotype = %d\tnloci = %d\tnallsum = %d\n",nG,nGenotype,nloci,nallsum);
 
  popsize_i0 = popsize; /* starting population size prior to expansion */

  /* simulate cell division */
  gd = gene_div_VNTR(freqG, x,nall, allele, nGenotype, nloci, popsize);
  pi = nucl_div_VNTR(freqG,x, nGenotype, nloci,thres);  // nb nb NB CHECK!
  if(flagg==1){
    if(!(divers = (double *) calloc((nGen+1),sizeof(double))) ||
       !(fitn = (double *) calloc((nGen+1),sizeof(double))) ||
       !(gdiv = (double *) calloc((nGen+1),sizeof(double)))){ 
	printf("Error in allocating memory at line 244\n\n"); 
	return NULL;
    }
    divers[0] = pi;
    gdiv[0] = gd;
    fitn[0] = mean_fitness(fitness,freqG,nGenotype);
  }
  if(verbosity==1) printf("Initial average gene diversity: %f\n",gd);
  if(verbosity==1) printf("Initial nucleotide diversity: %f\n",pi);
  if(verbosity==1) printf("Initial mean fitness: %f\n",mean_fitness(fitness,freqG,nGenotype));
  /* fix_mut_sp(pos, popsize, nlocitot, nloci, freqM, 0, n1); */
  /* fix_genotype_sp(pos, popsize, nlocitot, nloci, freqG, 0, n1); */
  
  totpop = maxnr+1; /* set population counter */
  startTime = walltime(&clockZero); /* start clock */

  free(nall);
  free(allele);
  /* display initial allele frequency */
  Rprintf("nloci = %d\n",nloci);
  for(i = 1; i <= nGen; i++){
    if(verbosity==1) printf("********** Simulating generation %d **********\n\n",i);
    popsize_i = popsize*pow(expfactor,i); /* population size at generation i */
    i1 = i-1;
    if(verbosity==1) printf("select clones to reproduce...\n");
    /* select clones to reproduce */
   
    /* calculate sample probabilities for each genotype */
    pG = calc_gen_prob(freqG,fitness,nGenotype,popsize_i0);
   
    free(freqG);

    numG = (unsigned int *) calloc(nGenotype,sizeof(int));
    if(verbosity==1) printf("sample from multinomial distibution to simulate new genotype frequencies\n");
    /* sample from multinomial distibution to simulate new genotype frequencies */
    rmultinom(popsize_i, pG, nGenotype, numG); /* gsl_ran_multinomial (gBaseRand, nGenotype, popsize_i, pG, numG); */
    if(verbosity==1) printf("Number of genotypes in population prior to selection: %d\n",nGenotype);

    free(pG);
    /* update the number of alleles in the population after selection */
    /* i.e. removing non-contributing alleles from the population */
    /* nG2 = upd_geno_nG(numG,id,idG,nGenotype,nG); */ /* update number of mutations */

    nGenotype2 = upd_geno_nGen(numG,nGenotype); /* update number of genotypes */
    x2 = upd_geno_matr(x,numG,nloci,nGenotype,nGenotype2);/* update genotype matrix */
    free(x);


    /* printf("Number of genotypes in population after selection (nGenotype2): %d\nNumber of new mutations in population (nG2): %d\n",nGenotype2,nG2); */

    /* allocate memory */
    if(verbosity==1) printf("number of Genotypes after selection: %d\n",nGenotype2);

    if(verbosity==1) printf("Remove non-existing genotypes from the memory...\n");
    /* upd_param(numGtmp,idGtmp,fitntmp,numG,idG,fitness,nGenotype,nGenotype2); */
    k=0;

    numGtmp = (int *) calloc(nGenotype2,sizeof(int));
    idGtmp = (int *) calloc(nGenotype2,sizeof(int));
    fitntmp = (double *) calloc(nGenotype2,sizeof(double));

    /* printf("Remove non-existing genotypes from the memory...\n"); */
    /* copy data */
    for(j = 0; j < nGenotype; j++)
      if(numG[j] > 0){  // copy
	numGtmp[k] = (int)numG[j];
	idGtmp[k] = idG[j];
	fitntmp[k] = fitness[j];
	k++;
      }
    /* printf("number of copied points: %d, nGenotype2: %d, nGenotype: %d\n",k,nGenotype2,nGenotype); */

    free(fitness);
    free(idG);
    free(numG);

    l1 = 0;
    if(verbosity==1) printf("distribute new mutations over remaining genotypes\n");
    /* distribute new mutations over remaining genotypes */
    numG = (unsigned int *) calloc(nGenotype2,sizeof(int));
    if(nrmut[i1] > 0){
    
      freqtmp2 = get_freq(numGtmp,nGenotype2);
    /* printf("Temporate frequency:\n"); */
    /* show(freqtmp2,1,nGenotype2); */

      if(verbosity==1) printf("sample from multinomial distibution to simulate new genotype frequencies\n");
      /* sample from multinomial distibution to simulate new genotype frequencies */
      rmultinom(nrmut[i1], freqtmp2, nGenotype2, numG); /* gsl_ran_multinomial (gBaseRand, nGenotype2, nrmut[i1], freqtmp2, numG); */
      /* rmultinom(popsize_i, pG, nGenotype, numG); */ /* gsl_ran_multinomial (gBaseRand, nGenotype, popsize_i, pG, numG); */
      /* printf("Draw from multinomial distribution:\n"); */
      /* for(j = 0; j < nGenotype2; j++) printf("%u ",numG[j]); */
      /* printf("\n"); */
      free(freqtmp2);
    }
    /* for(k = 0; k < nloci; k++) nperm[k] = k+1; */
    if(verbosity==1) printf("step 2: assign mutations to genotypes...\n");
    /* step 2: assign mutations to genotypes */
 
    cnt = 0; /* counter for the number of mutations assigned */
    /* if(nrmut[i1] > totpop) nrmut[i1] = totpop; */ // NB NB NB
    xnew = (int *) calloc(nrmut[i1]*nloci,sizeof(int));
    
    for(j = 0; j < nGenotype2; j++){
      //  if(numG[j] > 0){
      for(k = 0; k < numG[j]; k++){ // mutate genotype j
	numGtmp[j]--; // reduce number of ancestral genotype
  	for(l = 0; l < nloci; l++) xnew[l*nrmut[i1]+cnt] = x2[l*nGenotype2+j]; // copy original genotype
  	//	for(l = 0; l < mut[mcnt]; l++){ // select loci to mutate

  	/* gsl_ran_shuffle (gBaseRand, nperm, nloci, sizeof (int)); */ // select loci to mutate
  	//for(l = 0; l < mut[mcnt]; l++){ // mutate loci l
  	  // basic VNTR mutation model: the stepwise mutation model, see Ohta & Kimura (1973)
	/* nl = runif(0.0,1.0); /\* nl = gsl_rng_uniform (gBaseRand)*nloci; *\/ */
	/* gd = (double) (nl*nloci); */
	/* nlint = (int) gd; */
	u = runif(0.0,1.0); /* gsl_rng_uniform (gBaseRand); */
	if(u < 0.5){ // reduce allele number (i.e. number of repeats)
	  if(xnew[posm2[u2]*nrmut[i1]+cnt] > 1)
	    xnew[posm2[u2]*nrmut[i1]+cnt]--;
	  else
	    xnew[posm2[u2]*nrmut[i1]+cnt]++;
	}
	else // increase allele number
	  xnew[posm2[u2]*nrmut[i1]+cnt]++;
	/* Rprintf("Mutation no %d, mutating genotype %d, loci %d, random draw %f, new allele %d\n",cnt+1,j+1,posm2[u2],u,xnew[posm2[u2]*nrmut[i1]+cnt]); */
	  //}
  	  //}
  	cnt++;
  	mcnt++;
	u2++;
      }
    }

    if(nrmut[i1] > 0) free(numG);
    /* check if mutant genotype already exist */
    if(verbosity==1) printf("check if mutant genotype already exist...\n");
    nGenotype = comp_dim(x2, xnew, nloci, nGenotype2, nrmut[i1]);

    /* allocate new memory */
    x3 = (int *) calloc(nGenotype*nloci,sizeof(int));
    numG2 = (int *) calloc(nGenotype,sizeof(int));
    idG2 = (int *) calloc(nGenotype,sizeof(int));
    fitness2 = (double *) calloc(nGenotype,sizeof(double));

    maxnr = concat_matr(x2,xnew, x3,numG2,numGtmp,fitness2,fitntmp,idG2,idGtmp,nloci,nGenotype2,nrmut[i1],nGenotype,maxnr);
    gennr = nGenotype;
    //maxnr += nGenotype-nGenotype2;
    free(x2);
    free(xnew);
    free(numGtmp);
    free(fitntmp);
    free(idGtmp);

    if(verbosity==1) printf("copy and remove zero contributing genotypes\n");
 
    /* copy and remove zero contributing genotypes */
 
    nGenotype2 = upd_geno_nGen(numG2,nGenotype); /* update number of genotypes */
  
    if(verbosity==1) printf("Number of genotypes in population after mutatations occur, after removing zero contributing genotypes: %d\nbefore it was: %d\n",nGenotype2,nGenotype);

    /* allocate memory */
  
    idG = (int *) calloc(nGenotype2,sizeof(int));
    fitness = (double *) calloc(nGenotype2,sizeof(double));
    freqG = (double *) calloc(nGenotype2,sizeof(double));
    x = (int *) calloc(nGenotype2*nloci,sizeof(int));
    /* copy */

    k1 = 0; /* counter */
    l = 0;

    for(j = 0; j < nGenotype; j++){
      if(numG2[j] > 0){
  	freqG[k1] = (double)numG2[j]/popsize_i;
  	idG[k1] = idG2[j];
  	fitness[k1] = fitness2[j];
  	tmp0 = 1;
	for(k = 0; k < nloci; k++) x[k*nGenotype2+k1] = x3[k*nGenotype+j];
  	/* if(i==nGen && freqG[k1] > thresf){  */
	/*   fprintf(fp3,"%d\t%d\t",i,idG[k1]); */
	/*   for(k = 0; k < nloci; k++) fprintf(fp3,"%d\t",x[k*nGenotype+k1]); */
	/*   fprintf(fp3,"%f\n",freqG[k1]); */
	/* } */
  	/* if(tmp0==1 && i==nGen) */
  	/*   fprintf(fp4,"%d\t%d\t%d\t%f\n",i,idG[k1],-1,fitness[k1]); */
  	k1++;
      }
    }

    free(numG2);
    free(fitness2);
    free(idG2);
    free(x3);
    /* printf("freqG:\n"); */
    /* show(freqG,1,nGenotype2); */
  
    /* printf("\nfitness:\n"); */
    /* show(fitness,1,nGenotype2); */
    /* printf("idG:\n"); */
    /* showint(idG,1,nGenotype2); */

    popsize_i0 = popsize_i;
    nGenotype = nGenotype2;
 
  /* printf("nG = %d\tnGenotype = %d\tnlocitot = %d\n",nG,nGenotype,nlocitot); */
    if(verbosity==1) printf("update allele frequency...\n");
   
    nall = num_all(x,nloci,nGenotype); /* calculate the number of alleles per loci */
    nallsum = sum_VNTR_alleles(nall, nloci);
    allele = get_VNTR_alleles(x, nallsum, nloci, nGenotype);
    
    /* update allele frequency */
    /* print_genetic(i,idG,pos,id,fitness,fp4); */
    /* update_mut_freq_VNTR(nall,x,allele,freqG,nGenotype,nloci,i,fp1); *///update_mut_freq(pos,id,freqG,idG, nG,nGenotype, nlocitot,i,fp1);
    
    gd = gene_div_VNTR(freqG, x,nall, allele, nGenotype, nloci,popsize_i0);
    free(nall);
    free(allele);
    if(verbosity==1) printf("update mean fitness...\n");
    mf = mean_fitness(fitness,freqG,nGenotype);
    if(verbosity==1) printf("update nucleotide diversity...\n");
    pi = nucl_div_VNTR(freqG,x, nGenotype, nloci,thres);  // nb nb NB CHECK!
    if(flagg==1){
      gdiv[i] = gd;
      divers[i] = pi;
      fitn[i] = mf;
    }
    if(verbosity==1) printf("Average fitness: %f\n",mf);
    if(verbosity==1) printf("Average pair-wise nucleotide diversity: %f\n",pi);
    if(verbosity==1) printf("Average gene diversity: %f\n",gd);
   
  }
  /* free(idG); */
  free(fitness);
  /* free(x); */
  /* printf("thresf = %f, cnt1 = %d, cnt2 = %d\n",thresf,cnt1,cnt2); */
  elapsedTime = walltime(&startTime)/60.0;
  /* if(verbosity==1)  */
  printf("Elapsed time (in minutes): %f\n",elapsedTime);
  printf("Average fitness after %d generations: %f\n",nGen,mf);
  printf("Average pair-wise nucleotide diversity after %d generations: %f\n",nGen,pi);
  printf("Average gene diversity after %d generations: %f\n",nGen,gd);
  /* if(sav==1) fclose(fp1); */
  /* fclose(fp3); */
  /* /\* fclose(fp4); *\/ */
  /* /\* output = concat(freqG,fitness,nGenotype,nlocitot,pos,id,idG,nG); *\/ */
  /* /\* if(verbosity==1) printf("freqG:\n"); *\/ */
  /* /\* if(verbosity==1) show(freqG,1,nGenotype); *\/ */
  /* free(posm2); */
  /* free(mut); */
  /* Filter data */
  nGenotype2 = find_no_of_abundant_clones(freqG,nGenotype,thresf);
  Rprintf("Reducing number of clones with lower frequency than %f...\nNumber of clones prior to filtering: %d, number of clones after filtering: %d\n",thresf,nGenotype,nGenotype2);
  tmpdev = find_normalizing_constant(freqG,nGenotype,thresf);
  PutRNGstate();
  
  /* return output; */
  /* copy genotype x sample matrix */
  SEXP outX = PROTECT(allocVector(INTSXP,nGenotype2*(nloci+1)));
  SEXP xdim = PROTECT(allocVector(INTSXP,1));
  SEXP ydim = PROTECT(allocVector(INTSXP,1));
  SEXP freq2 = PROTECT(allocVector(REALSXP,nGenotype2));
  SEXP fitn2 = PROTECT(allocVector(REALSXP,nGenotype2));
  SEXP posnames = PROTECT(allocVector(INTSXP,nloci));

  int *xn,*yn,*posna; // dimension
  double *freq3,*fitn3;
  int *xab; // output matrix of allelic profiles in population
  
  xab = INTEGER(outX);
  xn = INTEGER(xdim);
  yn = INTEGER(ydim);
  fitn3 = REAL(fitn2);
  freq3 = REAL(freq2);
  posna = INTEGER(posnames);

  xn[0] = nGenotype2;
  yn[0] = nloci + 3;
  /* Rprintf("posi:\n"); */
  for(i = 0; i < nloci; i++){
    posna[i] = posi[i];
    /* Rprintf("%d ",posi[i]); */
  }
  /* Rprintf("\n"); */
  /* Rprintf("nr\tfreq\tfitness\n"); */
  k = 0;
  for(i = 0; i < nGenotype; i++){
    if(freqG[i] >= thresf){
      freq3[k] = freqG[i]/tmpdev;
      fitn3[k] = 1.;
      /* Rprintf("%d\t%f\t%f\n",i+1,freq3[i],1.); */
      for(j = 0; j < nloci; j++)
	xab[(j+1)*nGenotype2+k] = x[j*nGenotype+i];
      xab[k] = idG[i];
      k++;
    }

  }
  if(k!=nGenotype2)
    Rprintf("Warning!! k != nGenotype2....\nk = %d, nGenotype2 = %d,nGenotype = %d\n",k,nGenotype2,nGenotype);

  if(flagg==0){
    SEXP vec = PROTECT((allocVector(VECSXP,6)));
    SET_VECTOR_ELT(vec,0,outX);
    SET_VECTOR_ELT(vec,1,posnames);
    SET_VECTOR_ELT(vec,2,freq2);
    SET_VECTOR_ELT(vec,3,fitn2);
    SET_VECTOR_ELT(vec,4,xdim);
    SET_VECTOR_ELT(vec,5,ydim);
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
double find_normalizing_constant(double *freqG, int nGenotype, double thresf){

  int i;
  double g=0.;

  for(i = 0; i < nGenotype; i++)
    if(freqG[i] >= thresf)
      g += freqG[i];
  
  return g;

}
int find_no_of_abundant_clones(double *freqG, int nGenotype, double thresf){

  int i,nG = 0;
  
  for(i = 0; i < nGenotype; i++)
    if(freqG[i] >= thresf)
      nG++;
  
  return nG;

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
     printf("%f ",M[j*r+i]);
   printf("\n");
 }


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
      if(x[j*nGenotype+i]==0) nG++;

  return nG;
}
void sparse_gen(int *x,int *pos,int *id,double *freqG,int nG, int nGenotype,int nloci){

  int i,j,k = 0;
  
  for(i = 0; i < nGenotype; i++)
    for(j = 0; j < nloci; j++)
      if(x[j*nGenotype+i]==0){
	id[k] = i;
	pos[k] = j;
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
int find_max(int *nall,int nloci){

  int max=0,i;

  for(i = 0; i < nloci; i++)
    if(nall[i]>max)
      max = nall[i];
    
  return max;

}
void update_mut_freq_VNTR(int *nall,int *x, int *allele,double *freqG,int nGenotype, int nloci,int gen,FILE *fp){

  // save allele frequency on format: generation, locus, allele variant, frequency 
  
  int i,j,k,cnt = 0;
  double tmp,freq;

  for(i = 0; i < nloci; i++){
    for(j = 0; j < nall[i]; j++){
      freq = 0.;
      for(k = 0; k < nGenotype; k++)
	if(allele[cnt]==x[i*nGenotype+k]){ 
	  freq+= freqG[k];
	  /* printf("hit: allele[%d] = %d, loci = %d, allele no = %d,genotype nr = %d, freqG[%d] = %f, freq = %f\n",cnt,allele[cnt],i,j,k,k,freqG[k],freq); */
	}
      fprintf(fp,"%d\t%d\t%d\t%f\n",gen,i,allele[cnt],freq);
      cnt++;
    }
  }

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
  
  for(i = 0; i < nGenotype; i++)
    sum += fitness[i]*freqG[i]*popsize;

  for(i = 0; i < nGenotype; i++) dG[i] = freqG[i]*popsize*fitness[i]/sum;

  return dG;

}
double mean_fitness(double *fitness, double *freqG, int nGenotype){

  int i;
  double tmp=0.;

  for(i = 0; i < nGenotype; i++) tmp += fitness[i]*freqG[i];

  return tmp;

}
double nucl_div_VNTR(double *freqG, int *x, int nGenotype, int nloci, double thres){
  // computes average pair-wise nucleotide diversity for VNTRs

  int i,j,pi_ij,n1,k;
  double pi = 0.,sum=0.;

  if(nGenotype > 1){
    for(i = 1; i < nGenotype; i++) 
      for(j = 0; j < i; j++){
	n1 = 0;
	if(freqG[i] < thres || freqG[j] < thres)
	  pi_ij = 0.;	
	else{
	  for(k = 0; k < nloci; k++) 
	    if(x[k*nGenotype+i]!=x[k*nGenotype+j]) n1++;
	}
	sum += (double)freqG[i]*freqG[j]*n1;
      }
    pi = (double)(nGenotype/(nGenotype-1.))*sum;
  }
  else
    pi = 0.;
  
  return pi;

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
void read_gen(FILE *fps,int *x, double *freqG, double *fitness,int nGenotype,int nloci){

  int i,j,tmpint;
  float tmpfreq,tmpfit;

  for(i = 0; i < nGenotype; i++){
    fscanf(fps,"%d",&tmpint);
    for(j = 0; j < nloci; j++)
      fscanf(fps,"%d",&x[j*nGenotype+i]);
    fscanf(fps,"%f",&tmpfreq);
    fscanf(fps,"%f",&tmpfit);
    freqG[i] = (double)tmpfreq;
    fitness[i] = (double)tmpfit;
  }


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

  newfreq = (double *) malloc(nGenotype2*sizeof(double));

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
  tmp = (int *) calloc(nGenotype,sizeof(int));

  for(i = 0; i < nGenotype; i++){ 
    for(j = 0; j < nloci; j++)
      tmpx = (int) inX[(j+1)*nGenotype+i];
    freqG[i] = inX[nl1*nGenotype+i];
    tmpfit = inX[nl2*nGenotype+i];
  }
  rmultinom(popsize,freqG,nGenotype,tmp);
  /* gsl_ran_multinomial(gBaseRand, nGenotype, popsize, freqG, tmp); */
  /* printf("original frequency:\n"); */
  /* show(freqG,1,nGenotype); */
  for(i = 0; i < nGenotype; i++) sum += tmp[i];

  for(i = 0; i < nGenotype; i++) newfreq[i] = (double)tmp[i]/sum;

  /* printf("new frequency after draws from multinomial distribution:\n"); */
  /* show(newfreq,1,nGenotype); */
  free(tmp);

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
int sum_VNTR_alleles(int *nall, int nloci){

  int i,nallsum=0;

  for(i = 0; i < nloci; i++) nallsum += nall[i];

  return nallsum;

}
int *get_VNTR_alleles(int *x, int nallsum, int nloci, int nGenotype){

  int i,j,*allele,*tmp,k1,k,hit,cnt=0,in;

  allele = (int *) calloc(nallsum,sizeof(int));
  tmp = (int *) calloc(nGenotype,sizeof(int));

  for(i = 0; i < nloci; i++){
    k1 = 1;
    in = i*nGenotype;
    tmp[0] = x[in];
    allele[cnt] = tmp[0];
    cnt++;
    for(j = 1; j < nGenotype; j++){
      hit = 0;
      for(k = 0; k < k1; k++) 
	if(tmp[k]==x[in+j])
	  hit=1;
      if(hit==0){
	tmp[k1] = x[in+j];
	allele[cnt] = tmp[k1];
	cnt++;
	k1++;
      }
    }
    /* printf("tmp at loci %d\n",i); */
    /* showint(tmp,1,k1); */
  }
  /* printf("allocated memory = %d, indexed = %d\n",nallsum,cnt); */
  free(tmp);

  return allele;

}
int *upd_no_alleles(int *numG,int *x,int nGenotype,int nloci){

  int i,j,*nall,k1,*tmp,tmp1,k,l;

  nall = (int *) calloc(nloci,sizeof(int));
  tmp = (int *) calloc(nGenotype,sizeof(int));

  for(i = 0; i < nloci; i++){
    k1 = 1;
    for(l = 0; l < nGenotype; l++)
      if(numG[l] > 0){ 
	tmp[0] = x[i*nGenotype+l]; // Check!
	l = nGenotype; /* exit loop */
      }
    for(j = 0; j < nGenotype; j++)
      if(numG[j] > 0){
	tmp1 = 0;
	for(k = 0; k < k1; k++){
	  if(tmp[k]==x[i*nGenotype+j]) // allele already existing 
	  tmp1 = 1;	
	}
	if(tmp1==0){
	  tmp[k1]=x[i*nGenotype+j];
	  k1++;
	}
      }
    nall[i] = k1;
  }

  free(tmp);

  return nall;

}
int *num_all(int *x,int nloci,int nGenotype){
  // return the number of alleles per loci
  int *nall,*tmp,i,j,k,k1,tmp1;

  nall = (int *) calloc(nloci,sizeof(int));
  tmp = (int *) calloc(nGenotype,sizeof(int));

  for(i = 0; i < nloci; i++){
    k1 = 1;
    tmp[0] = x[i*nGenotype]; // Check!
 
    for(j = 1; j < nGenotype; j++){
      tmp1 = 0;
      for(k = 0; k < k1; k++){
	if(tmp[k]==x[i*nGenotype+j]) // allele already existing 
	  tmp1 = 1;	
      }
      if(tmp1==0){
	tmp[k1]=x[i*nGenotype+j];
	k1++;
      }
    }
    nall[i] = k1;
  }

  free(tmp);

  return nall;

}
int *upd_geno_matr(int *x,int *numG,int nloci,int nGenotype, int nGenotype2){

  int i,j, *x2,k = 0;

  x2 = (int *) malloc(nloci*nGenotype2*sizeof(int));

  for(i = 0; i < nGenotype; i++)
    if(numG[i] > 0){
      for(j = 0; j < nloci; j++)
	x2[j*nGenotype2+k] = x[j*nGenotype+i];
      k ++;
    }
  
  return x2;

}
void upd_param(int *numGtmp,int *idGtmp,double *fitntmp,unsigned int *numG,int *idG,double *fitness,int nGenotype, int nGenotype2){
    
  int j,k=0;

  numGtmp = (int *) calloc(nGenotype2,sizeof(int));
  idGtmp = (int *) calloc(nGenotype2,sizeof(int));
  fitntmp = (double *) calloc(nGenotype2,sizeof(double));

  /* printf("Remove non-existing genotypes from the memory...\n"); */
  /* copy data */
  for(j = 0; j < nGenotype; j++)
    if(numG[j] > 0){  // copy
      numGtmp[k] = (int)numG[j];
      idGtmp[k] = idG[j];
      fitntmp[k] = fitness[j];
      k++;
    }
  /* printf("number of copied points: %d, nGenotype2: %d, nGenotype: %d\n",k,nGenotype2,nGenotype); */
  /* printf("in upd_param\n"); */
  /* showint(numGtmp,1,nGenotype2); */
  free(fitness);
  free(idG);
  free(numG);

}
int comp_dim(int *x2,int *xnew,int nloci,int nGenotype2,int nmut){

  int j,k,l,flag,flag2,dim=nGenotype2,*tmp,ntot=nGenotype2+nmut;

  tmp = (int *) malloc(ntot*nloci*sizeof(int));

  for(k = 0; k < nGenotype2; k++)
    for(l = 0; l < nloci; l++) tmp[l*ntot+k] = x2[l*nGenotype2+k];

  /* step 1: determine the dimensionality of x (i.e. the number of rows) */
  for(j = 0; j < nmut; j++){
    flag2 = 0;
    for(k = 0; k < dim; k++){ 
      flag = 0;
      for(l = 0; l < nloci; l++){
	if(xnew[l*nmut+j]!=tmp[l*ntot+k]) flag = 1; /* not equal */
      }
      if(flag==0){ 
	flag2 = 1; /* equal */
      }
    }
    if(flag2==0){ /* if genotype is unique, add it to the population */
      for(l = 0; l < nloci; l++)
	tmp[l*ntot+dim] = xnew[l*nmut+j];
      dim++;
    }
  }

  free(tmp);

  return dim;

}
int concat_matr(int *x2,int *xnew, int *x,int *numG,int *numGtmp,double *fitness,double *fitntmp, int *idG,int *idGtmp,int nloci,int nGenotype2,int nmut,int dim,int maxnr){

  int j,k,l,flag,flag2,i = nGenotype2,*tmp,i1 = maxnr;

  /* printf("old genotype vector:\n"); */
  /* showint(x2,nGenotype2,nloci); */
 
  /* step 0: copy old genotypes */
  for(k = 0; k < nGenotype2; k++){
    for(l = 0; l < nloci; l++) x[l*dim+k] = x2[l*nGenotype2+k];
    numG[k] = numGtmp[k];
    fitness[k] = fitntmp[k];
    idG[k] = idGtmp[k];
  }
  //i1 = idGtmp[nGenotype2-1]+1;

  for(k = nGenotype2; k < dim; k++) numG[k] = 1;
  /* step 1: determine the dimensionality of x (i.e. the number of rows) */
  for(j = 0; j < nmut; j++){
    flag2 = 0;
    for(k = 0; k < i; k++){ 
      flag = 0;
      for(l = 0; l < nloci; l++){
	if(xnew[l*nmut+j]!=x[l*dim+k]) flag = 1; /* not equal */
      }
      if(flag==0){ 
	flag2 = 1; /* equal */
	numG[k] += 1;
      }
    }
    if(flag2==0){ /* if genotype is unique, add it to the population */
      for(l = 0; l < nloci; l++)
	x[l*dim+i] = xnew[l*nmut+j];
      fitness[i] = 1.;
      idG[i] = i1;
      i++;
      i1++;
    }
  }
 
  /* printf("new genotype vector:\n"); */
  /* showint(x,dim,nloci); */

  /* printf("new bacterial number vector:\n"); */
  /* showint(numG,1,dim); */
  return i1;
}
double gene_div_VNTR(double *freqG, int *x, int *nall, int *allele, int nGenotype, int nloci,int popsize){
  //gd = gene_div_VNTR(freqG, x,nall, allele, nGenotype, nloci);
  double gd=0.;
  // save allele frequency on format: generation, locus, allele variant, frequency 
  int i,j,k,cnt = 0;
  double tmp,freq,sum,sl = 0.,ps1;

  ps1 = (double)popsize/(popsize-1);


  if(nGenotype > 1){
    for(i = 0; i < nloci; i++){
      if(nall[i] > 1){
	sum = 0.;
	for(j = 0; j < nall[i]; j++){
	  freq = 0.;
	  for(k = 0; k < nGenotype; k++)
	    if(allele[cnt]==x[i*nGenotype+k]){ 
	      freq+= freqG[k];
	      /* printf("hit: allele[%d] = %d, loci = %d, allele no = %d,genotype nr = %d, freqG[%d] = %f, freq = %f\n",cnt,allele[cnt],i,j,k,k,freqG[k],freq); */
	    }
	  sum += freq*freq;
	  cnt++;
	}
	sl += (double)(1.-sum)*ps1; //(nall[i]/(nall[i]-1.));
      }
    }
  }    

  gd = (double)sl/nloci;

  /* printf("gd = %f\tsl = %f\tsum = %f\tnloci = %d\n",gd,sl,sum,nloci); */

  return gd;

}
int find_max2(double *inX, int nGenotype){
  // find maximum clone ID number
  int n = 0,i;

  for(i = 0; i < nGenotype; i++){
    if(inX[i] > n)
      n = (int)inX[i];

  }
  return n;

}
/* double *concat(double *freqG,double *fitness,int nGenotype2,int nlocitot,int *pos,int *id,int *idG,int nG){ */

/*   double *output; */
/*   int size = nGenotype2*(nlocitot+3),nl1 = nlocitot+1,nl2 = nlocitot+2; */
/*   int i,j,k; */

/*   output = (double *) malloc(size*sizeof(double)); */

/*   for(i = 0; i < nGenotype2; i++){ */
/*     output[i] = idG[i]; */
/*     for(j = 0; j < nlocitot; j++) output[(j+1)*nGenotype2+i] = 1.; */
/*     output[nl1*nGenotype2+i] = freqG[i]; */
/*     output[nl2*nGenotype2+i] = fitness[i]; */
/*   } */
/*   for(i = 0; i < nG; i++) output[pos[i]*nGenotype2+id[i]] = 0.; */

/*   return output; */

/* } */
/* void print_genetic(int nGen,int *idG,int *pos,int *id,double *fitness,FILE *fp4){ */




/* } */
