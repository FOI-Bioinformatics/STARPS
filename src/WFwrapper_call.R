WFwrapper<-function(nGen=1000,popsize=1000,expfactor=2,mutrate=5e-10,size=1800000,verbosity=1,p1=0.4,p2=0.5,mu=0.1,initPop,thresn=0.0001,thresf=0.000001,marker=1,fitn.eff=NULL,flag=0,simple=0)
{
# Input wrapper file for the Wright-Fisher simulation program STARPS
# nGen is number of generations per passage                                       
# popsize is the initial population size
# expfact is the expansion factor for bacterial growth
# mutrate is the mutation rate, per genome per cell division per base
# size is the bacterial or viral genome size in base pairs
# verbosity: should additional information be displayed on the screen: 0 for no additional info, 1 for additional info                         
# p1 is the probability of a harmful mutation (on fitness)
# p2 is the probability of a neutral mutation (on fitness)
# mu is the average selective advantage (or disadvantage) of mutations (mean of exp dist)
# thresn is a threshold used to avoid unneccesary extensive computational times to calculate the nucleotide diversity
# thresf is a threshold used to avoid saving unneccesary large data (i.e skip low abundant genotypes and mutations)
# for computationally challenging runs, set saveadd to 0
# marker is a flag for simulating biallelic (i.e. SNP or indels) data (marker = 1) or VNTR data (marker = 2)
# flag is telling whether diversity calculations should be performed within the main loop
# simple is a flag for allowing multiple mutations at the same bacteria in the same generation (i.e. at the same copy of a genotype). Default is simple = 0, i.e. allowing for multiple mutations. Warning! Might increase the required computational time substantially!
#
# script to be used for STARPS simulation toolbox
#
# By Petter Lindgren and Jon Ahlinder
  cat(sprintf('Simulating an evolving population\n'))	
  nGen1<-nGen+1	
  nGenotype = nrow(initPop)
  nloci = ncol(initPop) - 3		
  initPop2<-as.numeric(unlist(as.vector(initPop)))	
  vsize<-length(initPop2)
  if(is.null(fitn.eff) & marker == 1){ 
    cat(sprintf('Generate fitness effects of novel variants...\n'))
    fitn.eff<-draweff(size,p1,p2,mu)
  }
 # cat(sprintf('colnames(initPop)\n'))
 # print(colnames(initPop))
  #storage.mode(x) <- storage.mode(y) <- "double"
  if(marker==1) pos<-as.numeric(colnames(initPop)[2:nGenotype])
  else pos<-seq(1,nloci)
#  cat(sprintf('initPop in function WFwrapper_call.R:\n'))
#  print(initPop)
#  cat(sprintf('pos in function WFwrapper_call.R:\n'))
#  print(pos)

  if(marker==1){
    cat(sprintf('Biallelic markers considered\n'))
    if(simple==0)
      dyn.load(file.path("/mnt/powervault/jonhall/Desktop/Forensics/STARPS/",paste("WrightFisher_call",.Platform$dynlib.ext,sep="")))		
    else
      dyn.load(file.path("/mnt/powervault/jonhall/Desktop/Forensics/STARPS/",paste("WrightFisher_call_simple",.Platform$dynlib.ext,sep="")))
    outX<-.Call("WF2",nGen,popsize,expfactor,mutrate,nloci,verbosity,nGenotype,flag,initPop2,pos,thresn,thresf,fitn.eff)
  }
  else{
    cat(sprintf('Multiallelic markers considered\n'))
    cat(sprintf('nloci = %d\n',nloci))
    dyn.load(file.path("/mnt/powervault/jonhall/Desktop/Forensics/STARPS/",paste("WrightFisher_call_VNTR",.Platform$dynlib.ext,sep="")))
    outX<-.Call("WF_VNTR",nGen,popsize,expfactor,mutrate,nloci,verbosity,nGenotype,flag,initPop2,pos,thresn,thresf)
  }
#SEXP WF_VNTR(SEXP nGenp, SEXP popsizep, SEXP expfactorp, SEXP mutratep, SEXP nlocip, SEXP verbosityp, SEXP nGenotypep, SEXP flagp, SEXP inXp, SEXP posp, SEXP thresp, SEXP thresfrp, SEXP sizep)
 # print(outX,digits=3)

 # mat=cbind(0:nGen,tmp$fit,tmp$div)
  #mat=cbind(0:nGen,tmp$fit,tmp$div,tmp$gdiv)
 # colnames(mat)=c("Generation","Fitness","Nucleotide diversity","Gene diversity")
#  cat(sprintf('outX:\n'))
 # print(outX)
  mat<-make_matrix(outX)
 # mat<-outX
  cat(sprintf('Done\n'))
  return(mat)
 # mat		
}
