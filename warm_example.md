Warm
================
2024-09-14

## Warm

Add some introduction here.

The following is an example to run warm on data from five populations:

``` r
#this is the output folder for intermediate result files
outfolder="./"
#number of permuation
nperm=1e5
#number of thread to use
nthread=5
print(nperm)
```

    ## [1] 1e+05

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

    ## Loading definition of pathway: Sat Sep 14 00:28:03 2024

    ## Loading summary statistics: Sat Sep 14 00:28:04 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_AFR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 00:28:04 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 00:28:04 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 00:28:04 2024

    ## Realigning allele information of reference: Sat Sep 14 00:28:04 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 00:28:04 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 00:28:05 2024

    ## Removing constant SNPs: Sat Sep 14 00:28:05 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 00:28:05 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 00:28:13 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 00:28:37 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 00:28:42 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 00:29:05 2024

    ## Calculating P and SE if not provided: Sat Sep 14 00:29:05 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 00:29:05 2024

    ## Recovering test statistics: Sat Sep 14 00:29:05 2024

    ## Setup completed: Sat Sep 14 00:29:06 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting group 1: Sat Sep 14 00:29:06 2024

    ## Permuting group 2: Sat Sep 14 00:29:07 2024

    ## Permuting group 3: Sat Sep 14 00:29:07 2024

    ## Permuting group 4: Sat Sep 14 00:29:07 2024

    ## Permuting group 5: Sat Sep 14 00:29:08 2024

    ## Permuting group 6: Sat Sep 14 00:29:09 2024

    ## Permuting group 7: Sat Sep 14 00:29:09 2024

    ## Permuting group 8: Sat Sep 14 00:29:10 2024

    ## Permuting group 9: Sat Sep 14 00:29:10 2024

    ## Permuting group 10: Sat Sep 14 00:29:10 2024

    ## Permuting group 11: Sat Sep 14 00:29:11 2024

    ## Permuting group 12: Sat Sep 14 00:29:11 2024

    ## Permuting group 13: Sat Sep 14 00:29:12 2024

    ## Permuting group 14: Sat Sep 14 00:29:12 2024

    ## Permuting group 15: Sat Sep 14 00:29:12 2024

    ## Permuting group 16: Sat Sep 14 00:29:13 2024

    ## Permuting group 17: Sat Sep 14 00:29:13 2024

    ## Permuting group 18: Sat Sep 14 00:29:13 2024

    ## Permuting group 19: Sat Sep 14 00:29:14 2024

    ## Permuting group 20: Sat Sep 14 00:29:14 2024

    ## Permuting group 21: Sat Sep 14 00:29:15 2024

    ## Permuting group 22: Sat Sep 14 00:29:15 2024

    ## Permuting group 23: Sat Sep 14 00:29:16 2024

    ## Permuting group 24: Sat Sep 14 00:29:16 2024

    ## Permuting group 25: Sat Sep 14 00:29:17 2024

    ## Permuting group 26: Sat Sep 14 00:29:17 2024

    ## Permuting group 27: Sat Sep 14 00:29:18 2024

    ## Permuting group 28: Sat Sep 14 00:29:18 2024

    ## Permuting group 29: Sat Sep 14 00:29:18 2024

    ## Permuting group 30: Sat Sep 14 00:29:18 2024

    ## Permuting group 31: Sat Sep 14 00:29:19 2024

    ## Permuting group 32: Sat Sep 14 00:29:19 2024

    ## Permuting group 33: Sat Sep 14 00:29:19 2024

    ## Permuting group 34: Sat Sep 14 00:29:20 2024

    ## Permuting group 35: Sat Sep 14 00:29:20 2024

    ## Permuting group 36: Sat Sep 14 00:29:21 2024

    ## Permuting group 37: Sat Sep 14 00:29:21 2024

    ## Permuting group 38: Sat Sep 14 00:29:22 2024

    ## Permuting group 39: Sat Sep 14 00:29:22 2024

    ## Permuting group 40: Sat Sep 14 00:29:23 2024

    ## Permuting group 41: Sat Sep 14 00:29:23 2024

    ## Permuting group 42: Sat Sep 14 00:29:24 2024

    ## Permuting group 43: Sat Sep 14 00:29:24 2024

    ## Permuting group 44: Sat Sep 14 00:29:25 2024

    ## Permuting group 45: Sat Sep 14 00:29:26 2024

    ## Permuting group 46: Sat Sep 14 00:29:26 2024

    ## Permuting group 47: Sat Sep 14 00:29:26 2024

    ## Permuting group 48: Sat Sep 14 00:29:27 2024

    ## Permuting group 49: Sat Sep 14 00:29:27 2024

    ## Permuting group 50: Sat Sep 14 00:29:28 2024

    ## Permuting group 51: Sat Sep 14 00:29:28 2024

    ## Permuting group 52: Sat Sep 14 00:29:29 2024

    ## Permuting group 53: Sat Sep 14 00:29:29 2024

    ## Permuting group 54: Sat Sep 14 00:29:30 2024

    ## Permuting group 55: Sat Sep 14 00:29:30 2024

    ## Permuting group 56: Sat Sep 14 00:29:31 2024

    ## Permuting group 57: Sat Sep 14 00:29:33 2024

    ## Permuting group 58: Sat Sep 14 00:29:33 2024

    ## Permuting group 59: Sat Sep 14 00:29:33 2024

    ## Permuting group 60: Sat Sep 14 00:29:34 2024

    ## Permuting group 61: Sat Sep 14 00:29:35 2024

    ## Permuting group 62: Sat Sep 14 00:29:35 2024

    ## Permuting group 63: Sat Sep 14 00:29:36 2024

    ## Permuting group 64: Sat Sep 14 00:29:36 2024

    ## Permuting group 65: Sat Sep 14 00:29:37 2024

    ## Permuting group 66: Sat Sep 14 00:29:38 2024

    ## Permuting group 67: Sat Sep 14 00:29:39 2024

    ## Permuting group 68: Sat Sep 14 00:29:39 2024

    ## Permuting group 69: Sat Sep 14 00:29:40 2024

    ## Permuting group 70: Sat Sep 14 00:29:40 2024

    ## Permuting group 71: Sat Sep 14 00:29:41 2024

    ## Permuting group 72: Sat Sep 14 00:29:41 2024

    ## Permuting group 73: Sat Sep 14 00:29:42 2024

    ## Permuting group 74: Sat Sep 14 00:29:42 2024

    ## Permuting group 75: Sat Sep 14 00:29:42 2024

    ## Permuting group 76: Sat Sep 14 00:29:43 2024

    ## Permuting group 77: Sat Sep 14 00:29:43 2024

    ## Permuting group 78: Sat Sep 14 00:29:44 2024

    ## Permuting group 79: Sat Sep 14 00:29:44 2024

    ## Permuting group 80: Sat Sep 14 00:29:45 2024

    ## Permuting group 81: Sat Sep 14 00:29:45 2024

    ## Permuting group 82: Sat Sep 14 00:29:45 2024

    ## Permuting group 83: Sat Sep 14 00:29:46 2024

    ## Permuting group 84: Sat Sep 14 00:29:46 2024

    ## Permuting group 85: Sat Sep 14 00:29:47 2024

    ## Permuting group 86: Sat Sep 14 00:29:47 2024

    ## Permuting group 87: Sat Sep 14 00:29:47 2024

    ## Permuting group 88: Sat Sep 14 00:29:48 2024

    ## Permuting group 89: Sat Sep 14 00:29:48 2024

    ## Permuting group 90: Sat Sep 14 00:29:49 2024

    ## Permuting group 91: Sat Sep 14 00:29:49 2024

    ## Permuting group 92: Sat Sep 14 00:29:50 2024

    ## Permuting group 93: Sat Sep 14 00:29:50 2024

    ## Permuting group 94: Sat Sep 14 00:29:51 2024

    ## Permuting group 95: Sat Sep 14 00:29:51 2024

    ## Permuting group 96: Sat Sep 14 00:29:52 2024

    ## Permuting group 97: Sat Sep 14 00:29:52 2024

    ## Permuting group 98: Sat Sep 14 00:29:52 2024

    ## Permuting group 99: Sat Sep 14 00:29:53 2024

    ## Permuting group 100: Sat Sep 14 00:29:53 2024

    ## Permutation completed: Sat Sep 14 00:29:54 2024

    ## Computing pathway p-value: Sat Sep 14 00:29:54 2024

    ## sARTP files for warm.gene with  AFR is done!

    ## Loading definition of pathway: Sat Sep 14 00:29:54 2024

    ## Loading summary statistics: Sat Sep 14 00:29:54 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_AMR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 00:29:54 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 00:29:54 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 00:29:54 2024

    ## Realigning allele information of reference: Sat Sep 14 00:29:55 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 00:29:55 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 00:29:55 2024

    ## Removing constant SNPs: Sat Sep 14 00:29:55 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 00:29:55 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 00:30:03 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 00:30:19 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 00:30:25 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 00:30:38 2024

    ## Calculating P and SE if not provided: Sat Sep 14 00:30:38 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 00:30:39 2024

    ## Recovering test statistics: Sat Sep 14 00:30:39 2024

    ## Setup completed: Sat Sep 14 00:30:39 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting group 1: Sat Sep 14 00:30:39 2024

    ## Permuting group 2: Sat Sep 14 00:30:40 2024

    ## Permuting group 3: Sat Sep 14 00:30:40 2024

    ## Permuting group 4: Sat Sep 14 00:30:40 2024

    ## Permuting group 5: Sat Sep 14 00:30:41 2024

    ## Permuting group 6: Sat Sep 14 00:30:42 2024

    ## Permuting group 7: Sat Sep 14 00:30:43 2024

    ## Permuting group 8: Sat Sep 14 00:30:43 2024

    ## Permuting group 9: Sat Sep 14 00:30:43 2024

    ## Permuting group 10: Sat Sep 14 00:30:44 2024

    ## Permuting group 11: Sat Sep 14 00:30:44 2024

    ## Permuting group 12: Sat Sep 14 00:30:45 2024

    ## Permuting group 13: Sat Sep 14 00:30:45 2024

    ## Permuting group 14: Sat Sep 14 00:30:46 2024

    ## Permuting group 15: Sat Sep 14 00:30:46 2024

    ## Permuting group 16: Sat Sep 14 00:30:46 2024

    ## Permuting group 17: Sat Sep 14 00:30:46 2024

    ## Permuting group 18: Sat Sep 14 00:30:47 2024

    ## Permuting group 19: Sat Sep 14 00:30:47 2024

    ## Permuting group 20: Sat Sep 14 00:30:47 2024

    ## Permuting group 21: Sat Sep 14 00:30:48 2024

    ## Permuting group 22: Sat Sep 14 00:30:48 2024

    ## Permuting group 23: Sat Sep 14 00:30:48 2024

    ## Permuting group 24: Sat Sep 14 00:30:49 2024

    ## Permuting group 25: Sat Sep 14 00:30:49 2024

    ## Permuting group 26: Sat Sep 14 00:30:49 2024

    ## Permuting group 27: Sat Sep 14 00:30:50 2024

    ## Permuting group 28: Sat Sep 14 00:30:50 2024

    ## Permuting group 29: Sat Sep 14 00:30:50 2024

    ## Permuting group 30: Sat Sep 14 00:30:50 2024

    ## Permuting group 31: Sat Sep 14 00:30:51 2024

    ## Permuting group 32: Sat Sep 14 00:30:51 2024

    ## Permuting group 33: Sat Sep 14 00:30:51 2024

    ## Permuting group 34: Sat Sep 14 00:30:51 2024

    ## Permuting group 35: Sat Sep 14 00:30:52 2024

    ## Permuting group 36: Sat Sep 14 00:30:52 2024

    ## Permuting group 37: Sat Sep 14 00:30:53 2024

    ## Permuting group 38: Sat Sep 14 00:30:53 2024

    ## Permuting group 39: Sat Sep 14 00:30:53 2024

    ## Permuting group 40: Sat Sep 14 00:30:53 2024

    ## Permuting group 41: Sat Sep 14 00:30:54 2024

    ## Permuting group 42: Sat Sep 14 00:30:54 2024

    ## Permuting group 43: Sat Sep 14 00:30:55 2024

    ## Permuting group 44: Sat Sep 14 00:30:55 2024

    ## Permuting group 45: Sat Sep 14 00:30:55 2024

    ## Permuting group 46: Sat Sep 14 00:30:56 2024

    ## Permuting group 47: Sat Sep 14 00:30:56 2024

    ## Permuting group 48: Sat Sep 14 00:30:57 2024

    ## Permuting group 49: Sat Sep 14 00:30:57 2024

    ## Permuting group 50: Sat Sep 14 00:30:57 2024

    ## Permuting group 51: Sat Sep 14 00:30:58 2024

    ## Permuting group 52: Sat Sep 14 00:30:58 2024

    ## Permuting group 53: Sat Sep 14 00:30:58 2024

    ## Permuting group 54: Sat Sep 14 00:30:58 2024

    ## Permuting group 55: Sat Sep 14 00:30:59 2024

    ## Permuting group 56: Sat Sep 14 00:30:59 2024

    ## Permuting group 57: Sat Sep 14 00:31:00 2024

    ## Permuting group 58: Sat Sep 14 00:31:00 2024

    ## Permuting group 59: Sat Sep 14 00:31:00 2024

    ## Permuting group 60: Sat Sep 14 00:31:01 2024

    ## Permuting group 61: Sat Sep 14 00:31:01 2024

    ## Permuting group 62: Sat Sep 14 00:31:02 2024

    ## Permuting group 63: Sat Sep 14 00:31:02 2024

    ## Permuting group 64: Sat Sep 14 00:31:03 2024

    ## Permuting group 65: Sat Sep 14 00:31:04 2024

    ## Permuting group 66: Sat Sep 14 00:31:04 2024

    ## Permuting group 67: Sat Sep 14 00:31:05 2024

    ## Permuting group 68: Sat Sep 14 00:31:05 2024

    ## Permuting group 69: Sat Sep 14 00:31:06 2024

    ## Permuting group 70: Sat Sep 14 00:31:06 2024

    ## Permuting group 71: Sat Sep 14 00:31:07 2024

    ## Permuting group 72: Sat Sep 14 00:31:07 2024

    ## Permuting group 73: Sat Sep 14 00:31:07 2024

    ## Permuting group 74: Sat Sep 14 00:31:08 2024

    ## Permuting group 75: Sat Sep 14 00:31:08 2024

    ## Permuting group 76: Sat Sep 14 00:31:08 2024

    ## Permuting group 77: Sat Sep 14 00:31:09 2024

    ## Permuting group 78: Sat Sep 14 00:31:09 2024

    ## Permuting group 79: Sat Sep 14 00:31:09 2024

    ## Permuting group 80: Sat Sep 14 00:31:10 2024

    ## Permuting group 81: Sat Sep 14 00:31:10 2024

    ## Permuting group 82: Sat Sep 14 00:31:10 2024

    ## Permuting group 83: Sat Sep 14 00:31:11 2024

    ## Permuting group 84: Sat Sep 14 00:31:11 2024

    ## Permuting group 85: Sat Sep 14 00:31:11 2024

    ## Permuting group 86: Sat Sep 14 00:31:11 2024

    ## Permuting group 87: Sat Sep 14 00:31:12 2024

    ## Permuting group 88: Sat Sep 14 00:31:13 2024

    ## Permuting group 89: Sat Sep 14 00:31:13 2024

    ## Permuting group 90: Sat Sep 14 00:31:13 2024

    ## Permuting group 91: Sat Sep 14 00:31:13 2024

    ## Permuting group 92: Sat Sep 14 00:31:14 2024

    ## Permuting group 93: Sat Sep 14 00:31:14 2024

    ## Permuting group 94: Sat Sep 14 00:31:14 2024

    ## Permuting group 95: Sat Sep 14 00:31:15 2024

    ## Permuting group 96: Sat Sep 14 00:31:15 2024

    ## Permuting group 97: Sat Sep 14 00:31:15 2024

    ## Permuting group 98: Sat Sep 14 00:31:16 2024

    ## Permuting group 99: Sat Sep 14 00:31:16 2024

    ## Permuting group 100: Sat Sep 14 00:31:17 2024

    ## Permutation completed: Sat Sep 14 00:31:17 2024

    ## Computing pathway p-value: Sat Sep 14 00:31:17 2024

    ## sARTP files for warm.gene with  AMR is done!

    ## Loading definition of pathway: Sat Sep 14 00:31:17 2024

    ## Loading summary statistics: Sat Sep 14 00:31:17 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_EAS1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 00:31:18 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 00:31:18 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 00:31:18 2024

    ## Realigning allele information of reference: Sat Sep 14 00:31:18 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 00:31:18 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 00:31:18 2024

    ## Removing constant SNPs: Sat Sep 14 00:31:18 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 00:31:18 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 00:31:27 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 00:31:34 2024

    ## Calculating P and SE if not provided: Sat Sep 14 00:31:34 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 00:31:34 2024

    ## Recovering test statistics: Sat Sep 14 00:31:34 2024

    ## Setup completed: Sat Sep 14 00:31:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting group 1: Sat Sep 14 00:31:34 2024

    ## Permuting group 2: Sat Sep 14 00:31:35 2024

    ## Permuting group 3: Sat Sep 14 00:31:35 2024

    ## Permuting group 4: Sat Sep 14 00:31:35 2024

    ## Permuting group 5: Sat Sep 14 00:31:36 2024

    ## Permuting group 6: Sat Sep 14 00:31:36 2024

    ## Permuting group 7: Sat Sep 14 00:31:37 2024

    ## Permuting group 8: Sat Sep 14 00:31:37 2024

    ## Permuting group 9: Sat Sep 14 00:31:37 2024

    ## Permuting group 10: Sat Sep 14 00:31:38 2024

    ## Permuting group 11: Sat Sep 14 00:31:38 2024

    ## Permuting group 12: Sat Sep 14 00:31:38 2024

    ## Permuting group 13: Sat Sep 14 00:31:38 2024

    ## Permuting group 14: Sat Sep 14 00:31:39 2024

    ## Permuting group 15: Sat Sep 14 00:31:39 2024

    ## Permuting group 16: Sat Sep 14 00:31:39 2024

    ## Permuting group 17: Sat Sep 14 00:31:39 2024

    ## Permuting group 18: Sat Sep 14 00:31:40 2024

    ## Permuting group 19: Sat Sep 14 00:31:40 2024

    ## Permuting group 20: Sat Sep 14 00:31:40 2024

    ## Permuting group 21: Sat Sep 14 00:31:41 2024

    ## Permuting group 22: Sat Sep 14 00:31:41 2024

    ## Permuting group 23: Sat Sep 14 00:31:41 2024

    ## Permuting group 24: Sat Sep 14 00:31:41 2024

    ## Permuting group 25: Sat Sep 14 00:31:42 2024

    ## Permuting group 26: Sat Sep 14 00:31:42 2024

    ## Permuting group 27: Sat Sep 14 00:31:42 2024

    ## Permuting group 28: Sat Sep 14 00:31:43 2024

    ## Permuting group 29: Sat Sep 14 00:31:43 2024

    ## Permuting group 30: Sat Sep 14 00:31:43 2024

    ## Permuting group 31: Sat Sep 14 00:31:43 2024

    ## Permuting group 32: Sat Sep 14 00:31:43 2024

    ## Permuting group 33: Sat Sep 14 00:31:44 2024

    ## Permuting group 34: Sat Sep 14 00:31:44 2024

    ## Permuting group 35: Sat Sep 14 00:31:44 2024

    ## Permuting group 36: Sat Sep 14 00:31:44 2024

    ## Permuting group 37: Sat Sep 14 00:31:45 2024

    ## Permuting group 38: Sat Sep 14 00:31:45 2024

    ## Permuting group 39: Sat Sep 14 00:31:45 2024

    ## Permuting group 40: Sat Sep 14 00:31:46 2024

    ## Permuting group 41: Sat Sep 14 00:31:46 2024

    ## Permuting group 42: Sat Sep 14 00:31:46 2024

    ## Permuting group 43: Sat Sep 14 00:31:46 2024

    ## Permuting group 44: Sat Sep 14 00:31:47 2024

    ## Permuting group 45: Sat Sep 14 00:31:47 2024

    ## Permuting group 46: Sat Sep 14 00:31:47 2024

    ## Permuting group 47: Sat Sep 14 00:31:48 2024

    ## Permuting group 48: Sat Sep 14 00:31:48 2024

    ## Permuting group 49: Sat Sep 14 00:31:48 2024

    ## Permuting group 50: Sat Sep 14 00:31:48 2024

    ## Permuting group 51: Sat Sep 14 00:31:49 2024

    ## Permuting group 52: Sat Sep 14 00:31:49 2024

    ## Permuting group 53: Sat Sep 14 00:31:49 2024

    ## Permuting group 54: Sat Sep 14 00:31:50 2024

    ## Permuting group 55: Sat Sep 14 00:31:50 2024

    ## Permuting group 56: Sat Sep 14 00:31:50 2024

    ## Permuting group 57: Sat Sep 14 00:31:51 2024

    ## Permuting group 58: Sat Sep 14 00:31:52 2024

    ## Permuting group 59: Sat Sep 14 00:31:52 2024

    ## Permuting group 60: Sat Sep 14 00:31:53 2024

    ## Permuting group 61: Sat Sep 14 00:31:53 2024

    ## Permuting group 62: Sat Sep 14 00:31:53 2024

    ## Permuting group 63: Sat Sep 14 00:31:54 2024

    ## Permuting group 64: Sat Sep 14 00:31:54 2024

    ## Permuting group 65: Sat Sep 14 00:31:55 2024

    ## Permuting group 66: Sat Sep 14 00:31:55 2024

    ## Permuting group 67: Sat Sep 14 00:31:56 2024

    ## Permuting group 68: Sat Sep 14 00:31:56 2024

    ## Permuting group 69: Sat Sep 14 00:31:56 2024

    ## Permuting group 70: Sat Sep 14 00:31:56 2024

    ## Permuting group 71: Sat Sep 14 00:31:57 2024

    ## Permuting group 72: Sat Sep 14 00:31:57 2024

    ## Permuting group 73: Sat Sep 14 00:31:58 2024

    ## Permuting group 74: Sat Sep 14 00:31:58 2024

    ## Permuting group 75: Sat Sep 14 00:31:59 2024

    ## Permuting group 76: Sat Sep 14 00:31:59 2024

    ## Permuting group 77: Sat Sep 14 00:32:00 2024

    ## Permuting group 78: Sat Sep 14 00:32:00 2024

    ## Permuting group 79: Sat Sep 14 00:32:01 2024

    ## Permuting group 80: Sat Sep 14 00:32:01 2024

    ## Permuting group 81: Sat Sep 14 00:32:01 2024

    ## Permuting group 82: Sat Sep 14 00:32:02 2024

    ## Permuting group 83: Sat Sep 14 00:32:02 2024

    ## Permuting group 84: Sat Sep 14 00:32:03 2024

    ## Permuting group 85: Sat Sep 14 00:32:04 2024

    ## Permuting group 86: Sat Sep 14 00:32:04 2024

    ## Permuting group 87: Sat Sep 14 00:32:04 2024

    ## Permuting group 88: Sat Sep 14 00:32:05 2024

    ## Permuting group 89: Sat Sep 14 00:32:06 2024

    ## Permuting group 90: Sat Sep 14 00:32:06 2024

    ## Permuting group 91: Sat Sep 14 00:32:07 2024

    ## Permuting group 92: Sat Sep 14 00:32:07 2024

    ## Permuting group 93: Sat Sep 14 00:32:08 2024

    ## Permuting group 94: Sat Sep 14 00:32:08 2024

    ## Permuting group 95: Sat Sep 14 00:32:08 2024

    ## Permuting group 96: Sat Sep 14 00:32:08 2024

    ## Permuting group 97: Sat Sep 14 00:32:09 2024

    ## Permuting group 98: Sat Sep 14 00:32:09 2024

    ## Permuting group 99: Sat Sep 14 00:32:10 2024

    ## Permuting group 100: Sat Sep 14 00:32:10 2024

    ## Permutation completed: Sat Sep 14 00:32:11 2024

    ## Computing pathway p-value: Sat Sep 14 00:32:11 2024

    ## sARTP files for warm.gene with  EAS is done!

    ## Loading definition of pathway: Sat Sep 14 00:32:11 2024

    ## Loading summary statistics: Sat Sep 14 00:32:11 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_EUR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 00:32:11 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 00:32:11 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 00:32:11 2024

    ## Realigning allele information of reference: Sat Sep 14 00:32:12 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 00:32:12 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 00:32:12 2024

    ## Removing constant SNPs: Sat Sep 14 00:32:12 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 00:32:12 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 00:32:20 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 00:32:33 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 00:32:38 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 00:32:49 2024

    ## Calculating P and SE if not provided: Sat Sep 14 00:32:49 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 00:32:49 2024

    ## Recovering test statistics: Sat Sep 14 00:32:49 2024

    ## Setup completed: Sat Sep 14 00:32:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting group 1: Sat Sep 14 00:32:50 2024

    ## Permuting group 2: Sat Sep 14 00:32:50 2024

    ## Permuting group 3: Sat Sep 14 00:32:50 2024

    ## Permuting group 4: Sat Sep 14 00:32:51 2024

    ## Permuting group 5: Sat Sep 14 00:32:51 2024

    ## Permuting group 6: Sat Sep 14 00:32:52 2024

    ## Permuting group 7: Sat Sep 14 00:32:52 2024

    ## Permuting group 8: Sat Sep 14 00:32:52 2024

    ## Permuting group 9: Sat Sep 14 00:32:53 2024

    ## Permuting group 10: Sat Sep 14 00:32:53 2024

    ## Permuting group 11: Sat Sep 14 00:32:53 2024

    ## Permuting group 12: Sat Sep 14 00:32:54 2024

    ## Permuting group 13: Sat Sep 14 00:32:54 2024

    ## Permuting group 14: Sat Sep 14 00:32:54 2024

    ## Permuting group 15: Sat Sep 14 00:32:55 2024

    ## Permuting group 16: Sat Sep 14 00:32:55 2024

    ## Permuting group 17: Sat Sep 14 00:32:55 2024

    ## Permuting group 18: Sat Sep 14 00:32:55 2024

    ## Permuting group 19: Sat Sep 14 00:32:56 2024

    ## Permuting group 20: Sat Sep 14 00:32:56 2024

    ## Permuting group 21: Sat Sep 14 00:32:56 2024

    ## Permuting group 22: Sat Sep 14 00:32:57 2024

    ## Permuting group 23: Sat Sep 14 00:32:57 2024

    ## Permuting group 24: Sat Sep 14 00:32:57 2024

    ## Permuting group 25: Sat Sep 14 00:32:57 2024

    ## Permuting group 26: Sat Sep 14 00:32:58 2024

    ## Permuting group 27: Sat Sep 14 00:32:58 2024

    ## Permuting group 28: Sat Sep 14 00:32:58 2024

    ## Permuting group 29: Sat Sep 14 00:32:58 2024

    ## Permuting group 30: Sat Sep 14 00:32:59 2024

    ## Permuting group 31: Sat Sep 14 00:32:59 2024

    ## Permuting group 32: Sat Sep 14 00:33:00 2024

    ## Permuting group 33: Sat Sep 14 00:33:00 2024

    ## Permuting group 34: Sat Sep 14 00:33:00 2024

    ## Permuting group 35: Sat Sep 14 00:33:01 2024

    ## Permuting group 36: Sat Sep 14 00:33:01 2024

    ## Permuting group 37: Sat Sep 14 00:33:02 2024

    ## Permuting group 38: Sat Sep 14 00:33:02 2024

    ## Permuting group 39: Sat Sep 14 00:33:02 2024

    ## Permuting group 40: Sat Sep 14 00:33:02 2024

    ## Permuting group 41: Sat Sep 14 00:33:03 2024

    ## Permuting group 42: Sat Sep 14 00:33:03 2024

    ## Permuting group 43: Sat Sep 14 00:33:04 2024

    ## Permuting group 44: Sat Sep 14 00:33:04 2024

    ## Permuting group 45: Sat Sep 14 00:33:04 2024

    ## Permuting group 46: Sat Sep 14 00:33:05 2024

    ## Permuting group 47: Sat Sep 14 00:33:05 2024

    ## Permuting group 48: Sat Sep 14 00:33:05 2024

    ## Permuting group 49: Sat Sep 14 00:33:06 2024

    ## Permuting group 50: Sat Sep 14 00:33:06 2024

    ## Permuting group 51: Sat Sep 14 00:33:06 2024

    ## Permuting group 52: Sat Sep 14 00:33:07 2024

    ## Permuting group 53: Sat Sep 14 00:33:07 2024

    ## Permuting group 54: Sat Sep 14 00:33:07 2024

    ## Permuting group 55: Sat Sep 14 00:33:07 2024

    ## Permuting group 56: Sat Sep 14 00:33:07 2024

    ## Permuting group 57: Sat Sep 14 00:33:08 2024

    ## Permuting group 58: Sat Sep 14 00:33:09 2024

    ## Permuting group 59: Sat Sep 14 00:33:09 2024

    ## Permuting group 60: Sat Sep 14 00:33:09 2024

    ## Permuting group 61: Sat Sep 14 00:33:10 2024

    ## Permuting group 62: Sat Sep 14 00:33:10 2024

    ## Permuting group 63: Sat Sep 14 00:33:11 2024

    ## Permuting group 64: Sat Sep 14 00:33:11 2024

    ## Permuting group 65: Sat Sep 14 00:33:12 2024

    ## Permuting group 66: Sat Sep 14 00:33:12 2024

    ## Permuting group 67: Sat Sep 14 00:33:13 2024

    ## Permuting group 68: Sat Sep 14 00:33:13 2024

    ## Permuting group 69: Sat Sep 14 00:33:14 2024

    ## Permuting group 70: Sat Sep 14 00:33:14 2024

    ## Permuting group 71: Sat Sep 14 00:33:15 2024

    ## Permuting group 72: Sat Sep 14 00:33:15 2024

    ## Permuting group 73: Sat Sep 14 00:33:15 2024

    ## Permuting group 74: Sat Sep 14 00:33:16 2024

    ## Permuting group 75: Sat Sep 14 00:33:16 2024

    ## Permuting group 76: Sat Sep 14 00:33:16 2024

    ## Permuting group 77: Sat Sep 14 00:33:16 2024

    ## Permuting group 78: Sat Sep 14 00:33:17 2024

    ## Permuting group 79: Sat Sep 14 00:33:17 2024

    ## Permuting group 80: Sat Sep 14 00:33:18 2024

    ## Permuting group 81: Sat Sep 14 00:33:18 2024

    ## Permuting group 82: Sat Sep 14 00:33:18 2024

    ## Permuting group 83: Sat Sep 14 00:33:19 2024

    ## Permuting group 84: Sat Sep 14 00:33:19 2024

    ## Permuting group 85: Sat Sep 14 00:33:20 2024

    ## Permuting group 86: Sat Sep 14 00:33:20 2024

    ## Permuting group 87: Sat Sep 14 00:33:20 2024

    ## Permuting group 88: Sat Sep 14 00:33:21 2024

    ## Permuting group 89: Sat Sep 14 00:33:21 2024

    ## Permuting group 90: Sat Sep 14 00:33:21 2024

    ## Permuting group 91: Sat Sep 14 00:33:22 2024

    ## Permuting group 92: Sat Sep 14 00:33:22 2024

    ## Permuting group 93: Sat Sep 14 00:33:22 2024

    ## Permuting group 94: Sat Sep 14 00:33:23 2024

    ## Permuting group 95: Sat Sep 14 00:33:23 2024

    ## Permuting group 96: Sat Sep 14 00:33:23 2024

    ## Permuting group 97: Sat Sep 14 00:33:24 2024

    ## Permuting group 98: Sat Sep 14 00:33:24 2024

    ## Permuting group 99: Sat Sep 14 00:33:25 2024

    ## Permuting group 100: Sat Sep 14 00:33:25 2024

    ## Permutation completed: Sat Sep 14 00:33:25 2024

    ## Computing pathway p-value: Sat Sep 14 00:33:25 2024

    ## sARTP files for warm.gene with  EUR is done!

    ## Loading definition of pathway: Sat Sep 14 00:33:25 2024

    ## Loading summary statistics: Sat Sep 14 00:33:26 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_SAS1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 00:33:26 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 00:33:26 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 00:33:26 2024

    ## Realigning allele information of reference: Sat Sep 14 00:33:26 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 00:33:26 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 00:33:26 2024

    ## Removing constant SNPs: Sat Sep 14 00:33:26 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 00:33:27 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 00:33:36 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 00:33:48 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 00:33:53 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 00:34:04 2024

    ## Calculating P and SE if not provided: Sat Sep 14 00:34:04 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 00:34:04 2024

    ## Recovering test statistics: Sat Sep 14 00:34:04 2024

    ## Setup completed: Sat Sep 14 00:34:04 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting group 1: Sat Sep 14 00:34:04 2024

    ## Permuting group 2: Sat Sep 14 00:34:05 2024

    ## Permuting group 3: Sat Sep 14 00:34:05 2024

    ## Permuting group 4: Sat Sep 14 00:34:05 2024

    ## Permuting group 5: Sat Sep 14 00:34:06 2024

    ## Permuting group 6: Sat Sep 14 00:34:07 2024

    ## Permuting group 7: Sat Sep 14 00:34:08 2024

    ## Permuting group 8: Sat Sep 14 00:34:08 2024

    ## Permuting group 9: Sat Sep 14 00:34:08 2024

    ## Permuting group 10: Sat Sep 14 00:34:08 2024

    ## Permuting group 11: Sat Sep 14 00:34:09 2024

    ## Permuting group 12: Sat Sep 14 00:34:09 2024

    ## Permuting group 13: Sat Sep 14 00:34:10 2024

    ## Permuting group 14: Sat Sep 14 00:34:10 2024

    ## Permuting group 15: Sat Sep 14 00:34:10 2024

    ## Permuting group 16: Sat Sep 14 00:34:11 2024

    ## Permuting group 17: Sat Sep 14 00:34:11 2024

    ## Permuting group 18: Sat Sep 14 00:34:11 2024

    ## Permuting group 19: Sat Sep 14 00:34:12 2024

    ## Permuting group 20: Sat Sep 14 00:34:12 2024

    ## Permuting group 21: Sat Sep 14 00:34:12 2024

    ## Permuting group 22: Sat Sep 14 00:34:13 2024

    ## Permuting group 23: Sat Sep 14 00:34:13 2024

    ## Permuting group 24: Sat Sep 14 00:34:13 2024

    ## Permuting group 25: Sat Sep 14 00:34:13 2024

    ## Permuting group 26: Sat Sep 14 00:34:14 2024

    ## Permuting group 27: Sat Sep 14 00:34:14 2024

    ## Permuting group 28: Sat Sep 14 00:34:15 2024

    ## Permuting group 29: Sat Sep 14 00:34:15 2024

    ## Permuting group 30: Sat Sep 14 00:34:15 2024

    ## Permuting group 31: Sat Sep 14 00:34:15 2024

    ## Permuting group 32: Sat Sep 14 00:34:16 2024

    ## Permuting group 33: Sat Sep 14 00:34:16 2024

    ## Permuting group 34: Sat Sep 14 00:34:16 2024

    ## Permuting group 35: Sat Sep 14 00:34:16 2024

    ## Permuting group 36: Sat Sep 14 00:34:17 2024

    ## Permuting group 37: Sat Sep 14 00:34:17 2024

    ## Permuting group 38: Sat Sep 14 00:34:17 2024

    ## Permuting group 39: Sat Sep 14 00:34:18 2024

    ## Permuting group 40: Sat Sep 14 00:34:18 2024

    ## Permuting group 41: Sat Sep 14 00:34:18 2024

    ## Permuting group 42: Sat Sep 14 00:34:19 2024

    ## Permuting group 43: Sat Sep 14 00:34:19 2024

    ## Permuting group 44: Sat Sep 14 00:34:20 2024

    ## Permuting group 45: Sat Sep 14 00:34:20 2024

    ## Permuting group 46: Sat Sep 14 00:34:20 2024

    ## Permuting group 47: Sat Sep 14 00:34:20 2024

    ## Permuting group 48: Sat Sep 14 00:34:21 2024

    ## Permuting group 49: Sat Sep 14 00:34:21 2024

    ## Permuting group 50: Sat Sep 14 00:34:21 2024

    ## Permuting group 51: Sat Sep 14 00:34:22 2024

    ## Permuting group 52: Sat Sep 14 00:34:22 2024

    ## Permuting group 53: Sat Sep 14 00:34:22 2024

    ## Permuting group 54: Sat Sep 14 00:34:22 2024

    ## Permuting group 55: Sat Sep 14 00:34:22 2024

    ## Permuting group 56: Sat Sep 14 00:34:23 2024

    ## Permuting group 57: Sat Sep 14 00:34:24 2024

    ## Permuting group 58: Sat Sep 14 00:34:24 2024

    ## Permuting group 59: Sat Sep 14 00:34:24 2024

    ## Permuting group 60: Sat Sep 14 00:34:25 2024

    ## Permuting group 61: Sat Sep 14 00:34:25 2024

    ## Permuting group 62: Sat Sep 14 00:34:26 2024

    ## Permuting group 63: Sat Sep 14 00:34:26 2024

    ## Permuting group 64: Sat Sep 14 00:34:26 2024

    ## Permuting group 65: Sat Sep 14 00:34:27 2024

    ## Permuting group 66: Sat Sep 14 00:34:27 2024

    ## Permuting group 67: Sat Sep 14 00:34:28 2024

    ## Permuting group 68: Sat Sep 14 00:34:28 2024

    ## Permuting group 69: Sat Sep 14 00:34:29 2024

    ## Permuting group 70: Sat Sep 14 00:34:29 2024

    ## Permuting group 71: Sat Sep 14 00:34:30 2024

    ## Permuting group 72: Sat Sep 14 00:34:30 2024

    ## Permuting group 73: Sat Sep 14 00:34:30 2024

    ## Permuting group 74: Sat Sep 14 00:34:31 2024

    ## Permuting group 75: Sat Sep 14 00:34:31 2024

    ## Permuting group 76: Sat Sep 14 00:34:31 2024

    ## Permuting group 77: Sat Sep 14 00:34:32 2024

    ## Permuting group 78: Sat Sep 14 00:34:32 2024

    ## Permuting group 79: Sat Sep 14 00:34:32 2024

    ## Permuting group 80: Sat Sep 14 00:34:33 2024

    ## Permuting group 81: Sat Sep 14 00:34:33 2024

    ## Permuting group 82: Sat Sep 14 00:34:33 2024

    ## Permuting group 83: Sat Sep 14 00:34:33 2024

    ## Permuting group 84: Sat Sep 14 00:34:34 2024

    ## Permuting group 85: Sat Sep 14 00:34:34 2024

    ## Permuting group 86: Sat Sep 14 00:34:34 2024

    ## Permuting group 87: Sat Sep 14 00:34:35 2024

    ## Permuting group 88: Sat Sep 14 00:34:35 2024

    ## Permuting group 89: Sat Sep 14 00:34:36 2024

    ## Permuting group 90: Sat Sep 14 00:34:36 2024

    ## Permuting group 91: Sat Sep 14 00:34:36 2024

    ## Permuting group 92: Sat Sep 14 00:34:37 2024

    ## Permuting group 93: Sat Sep 14 00:34:37 2024

    ## Permuting group 94: Sat Sep 14 00:34:37 2024

    ## Permuting group 95: Sat Sep 14 00:34:38 2024

    ## Permuting group 96: Sat Sep 14 00:34:38 2024

    ## Permuting group 97: Sat Sep 14 00:34:38 2024

    ## Permuting group 98: Sat Sep 14 00:34:39 2024

    ## Permuting group 99: Sat Sep 14 00:34:39 2024

    ## Permuting group 100: Sat Sep 14 00:34:40 2024

    ## Permutation completed: Sat Sep 14 00:34:40 2024

    ## Computing pathway p-value: Sat Sep 14 00:34:40 2024

    ## sARTP files for warm.gene with  SAS is done!

``` r
res.warm.FISHER<-warm.start.multiPop.gene(sARTP2.list=res_warm, method = "FISHER", nthread=nthread,
                                          delete.gene.files=FALSE, expand.fold = 0, use.ranks.pathway=FALSE)
```

    ## Computing pathway p-value: Sat Sep 14 00:34:40 2024

``` r
res.warm.META<-warm.start.multiPop.gene(sARTP2.list=res_warm, method = "META", nthread=nthread,
                                        delete.gene.files=FALSE, expand.fold = 0, use.ranks.pathway=FALSE)
```

    ## Computing pathway p-value: Sat Sep 14 00:34:41 2024

``` r
res.warm.MAX<-warm.start.multiPop.gene(sARTP2.list=res_warm, method = "MAX", nthread=nthread,
                                       delete.gene.files=FALSE, expand.fold = 0, use.ranks.pathway=FALSE)
```

    ## Computing pathway p-value: Sat Sep 14 00:34:42 2024

``` r
res.warm.WGTFISHER<-warm.start.multiPop.gene(sARTP2.list=res_warm, method = "WGTFISHER", nthread=nthread,
                                             delete.gene.files=FALSE, expand.fold = 0, use.ranks.pathway=FALSE)
```

    ## Computing pathway p-value: Sat Sep 14 00:34:43 2024

``` r
res.warm.ACAT<-warm.start.multiPop.gene(sARTP2.list=res_warm, method = "ACAT", nthread=nthread,
                                        delete.gene.files=TRUE, expand.fold = 0, use.ranks.pathway=FALSE)
```

    ## Computing pathway p-value: Sat Sep 14 00:34:57 2024

``` r
########### please update the output
path_pvalue[npop+(1:5)] = c(res.warm.FISHER$pathway.pvalue,res.warm.META$pathway.pvalue,
                            res.warm.MAX$pathway.pvalue,res.warm.WGTFISHER$pathway.pvalue,res.warm.ACAT$pathway.pvalue)
gene_pvalue[npop+(1:5)] = list(res.warm.FISHER$gene.pvalue,res.warm.META$gene.pvalue,
                               res.warm.MAX$gene.pvalue,res.warm.WGTFISHER$gene.pvalue,res.warm.ACAT$gene.pvalue)
names(path_pvalue) = c(text,"Fisher","Meta","Max","WGTFISHER","Acat")
names(gene_pvalue) = c(text,"Fisher","Meta","Max","WGTFISHER","Acat")

res = list(gene_pvalue=gene_pvalue,path_pvalue=path_pvalue)
names(res$gene_pvalue)
```

    ##  [1] "AFR"       "AMR"       "EAS"       "EUR"       "SAS"       "Fisher"   
    ##  [7] "Meta"      "Max"       "WGTFISHER" "Acat"

``` r
res$path_pvalue
```

    ##          AFR          AMR          EAS          EUR          SAS       Fisher 
    ## 0.6721282787 0.0148148519 0.0243447566 0.0038799612 0.1650383496 0.0002099979 
    ##         Meta          Max    WGTFISHER         Acat 
    ## 0.0002599974 0.0019599804 0.0001299987 0.0007849922
