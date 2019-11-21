# STARPS
A simulation tool for cultivation of bacterial populations according to a Wright-Fisher model of inheritance. 


## Usage

### Compilation
To compile the C code, the [Gnu Scientific Library](https://www.gnu.org/software/gsl/) is needed. 
```
gcc -o WrightFisher_call WrightFisher_call.c -Wall -I/usr/include -lm -lgsl -lgslcblas
```
### Running a simulation 
The file [WrightFisher_call.c](src/WrightFisher_call.c) is called within the R environment: 


...
source("/mnt/powervault/jonhall/Desktop/Forensics/STARPS/R/WFwrapper_call.R")
source("/mnt/powervault/jonhall/Desktop/Forensics/STARPS/R/draweff.R")
cat(sprintf('Generate fitness effects...\n'))
fitn.eff <- draweff(size = tsize)
cat(sprintf('Generate an initial population...\n'))
initPop<-matrix(0,2,4)

initPop[1:2,4]<-1
initPop[1,3]<-0.999
initPop[2,3]<-0.001
initPop[2,2]<-1
initPop[2,1]<-1
colnames(initPop)<-c("clone",157321,"frequency","fitness")

cat(sprintf('Start a cultivation...\n'))
init.pop <- WFwrapper(nGen=nGen, 
      popsize=popsize, 
      expfactor=expfactor, 
      mutrate=mutrate, 
      size=size, 
      verbosity=verbosity, 
      p1=p1, 
      p2=p2, 
      mu=mu, 
      initPop=initPop, 
      thresf=thresf, 
      marker=marker, 
      fitn.eff=fitn.eff, 
      simple=simple)
...

The input to WFwrapper_call.R is
- nGen is number of generations per passage                                       
- popsize is the initial population size
- expfact is the expansion factor for bacterial growth
- mutrate is the mutation rate, per genome per cell division per base
- size is the bacterial or viral genome size in base pairs
- verbosity: should additional information be displayed on the screen: 0 for no additional info, 1 for additional info                         
- p1 is the probability of a harmful mutation (on fitness)
- p2 is the probability of a neutral mutation (on fitness)
- mu is the average selective advantage (or disadvantage) of mutations (mean of exp dist)
- thresn is a threshold used to avoid unneccesary extensive computational times to calculate the nucleotide diversity
- thresf is a threshold used to avoid saving unneccesary large data (i.e skip low abundant genotypes and mutations) for computationally challenging runs, set saveadd to 0
- marker is a flag for simulating biallelic (i.e. SNP or indels) data (marker = 1) or VNTR data (marker = 2)
- flag is telling whether diversity calculations should be performed within the main loop
- simple is a flag for allowing multiple mutations at the same bacteria in the same generation (i.e. at the same copy of a genotype). Default is simple = 0, i.e. allowing for multiple mutations. Warning! Might increase the required computational time substantially!
