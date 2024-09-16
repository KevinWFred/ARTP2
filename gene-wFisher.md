gene-wFisher
================
2024-09-14

## Gene centric gene-wFisher method (see Fu et al. for more details on the method)

An example to run the gene-wFisher procedure on data from five
populations:

``` r
#this is the output folder for intermediate result files
outfolder="./"
#number of permuation
nperm=1e4
#number of thread to use
nthread=5
print(nperm)
```

    ## [1] 10000

``` r
print(nthread)
```

    ## [1] 5

``` r

library("ARTP2")

#data from the following ancestres 
text = c("AFR","AMR","EAS","EUR","SAS")

npop = length(text)
#number of case/control
ncase = nctrl = c(4000, 6000, 6000, 10000, 4000)/2

#summary data
studys <- vector("character", npop)
extdatafolder=system.file("extdata", package = "ARTP2")
#reference data
reference=NULL
for(i in 1:npop){
  studys[i]= paste(extdatafolder,"/study_", text[i],"1.txt.gz", sep = "")
  prefix=paste0(extdatafolder,"/",text[i],"_ref")
  reference=rbind(reference,data.frame(fam=paste0(prefix,".fam"),bim=paste0(prefix,".bim"),bed=paste0(prefix,".bed")))
}
#pathway file
pathway_file = paste0(extdatafolder,"/pathway_100_new.txt.gz")


## running three sARTP methods 
lambda <- 1
lambdas = rep(lambda,npop)

family <- 'binomial'

ncases <- split(ncase, seq_along(ncase))
nctrls <- split(nctrl, seq_along(nctrl))

path_pvalue = NULL
gene_pvalue = list()

running_time = NULL

## preparing sARTP files for warm.start.multiPop.gene
res_warm <- list()
s0 = 10000
jobid ="warm_example"
outdir=paste0("./", jobid)
if (!dir.exists(outdir)) dir.create(outdir)
#jobid <- system("echo $SLURM_JOB_ID", intern=TRUE)

for(j in 1:npop){
  ref1 = reference[j,]
  study1 = studys[j]
  
  optioni_warm = list(inspect.snp.n = 2, nperm = nperm, nthread = nthread, 
                      maf = .02, HWE.p = 0, gene.R2 = .9, chr.R2=0.9,
                      min.marg.p = 1e-20,group.gap = 9*10^6,
                      id.str = paste0(text[j],"_warm"), 
                      out.dir=paste0(outdir),
                      save.setup = FALSE, delete.gene.files = FALSE,
                      seed = s0+j)
  
  # Record the starting time
  start_time <- Sys.time()
  res_warm[[j]] <- sARTP(summary.files = study1, pathway_file, family, ref1, lambdas[j],
                         ncases[j], nctrls[j], options = optioni_warm)
  # Record the ending time
  end_time <- Sys.time()
  
  # Calculate the time taken
  running_time[j] <- end_time - start_time
  
  path_pvalue[j] = res_warm[[j]]$pathway.pvalue
  gene_pvalue[[j]] = res_warm[[j]]$gene.pvalue
  
  
  cat("sARTP files for warm.gene with ",text[j],"is done!\n")
}
```

``` r
res.warm.WGTFISHER$pathway.pvalue
```

    ## [1] 0.00019998
