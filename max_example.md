Max
================
2024-09-14

## Max

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
.libPaths(c("/data/BB_Bioinformatics/Kevin/tools/Rpackages",.libPaths()))

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

lambda.list = as.list(rep(lambda,npop))

family <- 'binomial'

ncases <- split(ncase, seq_along(ncase))
nctrls <- split(nctrl, seq_along(nctrl))
ncases.list <- split(ncases, seq_along(ncases))
nctrls.list <- split(nctrls, seq_along(ncases))

study.list = as.list(studys)
reference.list <- split(reference, 1:nrow(reference))
running_time = NULL

## preparing sARTP files for warm.start.multiPop.gene
s0 = 10000
jobid <- "max_example"
outdir=paste0("./", jobid)
if (!dir.exists(outdir)) dir.create(outdir)

options_max <- list()
for(j in 1:npop){
  options_max[[j]] <- list(maf = .02, HWE.p = 0, min.marg.p=1e-20,
                           gene.R2=0.9, chr.R2=0.9, huge.gene.R2=0.8, huge.chr.R2=0.8,
                           id.str = paste0(text[j],"-max"),
                           out.dir=outdir, seed = s0+j,
                           nthread = nthread, save.setup = FALSE)
 
}
options.merged1 = list(nperm = nperm, merge.snp.method = "MAX", 
                       id.str = paste0("MAX-"),
                       nthread = nthread, expand.fold = 0, 
                       out.dir=outdir,
                       seed = s0, inspect.snp.n=2, group.gap=9E6)

# start <- Sys.time()
res_temp <- sARTP.multiPop.SNP(summary.files.list=study.list, pathway=pathway_file, family=family, 
                               reference.list=reference.list,
                              lambda.list=lambda.list, ncases.list=ncases.list, ncontrols.list = nctrls.list, 
                              options.list = options_max, 
                              options.merged=options.merged1)
```

    ## Loading definition of pathway: Sat Sep 14 01:35:40 2024

    ## Loading summary statistics: Sat Sep 14 01:35:40 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_AFR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:35:40 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:35:40 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:35:40 2024

    ## Realigning allele information of reference: Sat Sep 14 01:35:41 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:35:41 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:35:41 2024

    ## Removing constant SNPs: Sat Sep 14 01:35:41 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:35:41 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:35:51 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 01:36:16 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 01:36:24 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:36:44 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:36:44 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:36:44 2024

    ## Loading definition of pathway: Sat Sep 14 01:36:45 2024

    ## Loading summary statistics: Sat Sep 14 01:36:45 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_AMR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:36:45 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:36:45 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:36:45 2024

    ## Realigning allele information of reference: Sat Sep 14 01:36:46 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:36:46 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:36:46 2024

    ## Removing constant SNPs: Sat Sep 14 01:36:46 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:36:46 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:36:57 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 01:37:13 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 01:37:23 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:37:34 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:37:34 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:37:34 2024

    ## Loading definition of pathway: Sat Sep 14 01:37:35 2024

    ## Loading summary statistics: Sat Sep 14 01:37:35 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_EAS1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:37:35 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:37:35 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:37:35 2024

    ## Realigning allele information of reference: Sat Sep 14 01:37:35 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:37:36 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:37:36 2024

    ## Removing constant SNPs: Sat Sep 14 01:37:36 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:37:36 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:37:47 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:37:54 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:37:54 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:37:54 2024

    ## Loading definition of pathway: Sat Sep 14 01:37:54 2024

    ## Loading summary statistics: Sat Sep 14 01:37:55 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_EUR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:37:55 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:37:55 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:37:55 2024

    ## Realigning allele information of reference: Sat Sep 14 01:37:55 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:37:55 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:37:55 2024

    ## Removing constant SNPs: Sat Sep 14 01:37:56 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:37:56 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:38:07 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 01:38:19 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 01:38:26 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:38:36 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:38:36 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:38:36 2024

    ## Loading definition of pathway: Sat Sep 14 01:38:36 2024

    ## Loading summary statistics: Sat Sep 14 01:38:37 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_SAS1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:38:37 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:38:37 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:38:37 2024

    ## Realigning allele information of reference: Sat Sep 14 01:38:37 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:38:37 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:38:37 2024

    ## Removing constant SNPs: Sat Sep 14 01:38:37 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:38:38 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:38:47 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 01:38:59 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 01:39:06 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:39:15 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:39:15 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:39:15 2024

    ## Recovering test statistics: Sat Sep 14 01:39:15 2024

    ## Recovering test statistics: Sat Sep 14 01:39:16 2024

    ## Recovering test statistics: Sat Sep 14 01:39:17 2024

    ## Recovering test statistics: Sat Sep 14 01:39:18 2024
    ## Recovering test statistics: Sat Sep 14 01:39:18 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 1: Sat Sep 14 01:39:21 2024

    ## Permutation completed: Sat Sep 14 01:39:22 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 1: Sat Sep 14 01:39:22 2024

    ## Permutation completed: Sat Sep 14 01:39:22 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 1: Sat Sep 14 01:39:22 2024

    ## Permutation completed: Sat Sep 14 01:39:22 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 1: Sat Sep 14 01:39:22 2024

    ## Permutation completed: Sat Sep 14 01:39:22 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 1: Sat Sep 14 01:39:22 2024

    ## Permutation completed: Sat Sep 14 01:39:22 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 2: Sat Sep 14 01:39:23 2024

    ## Permutation completed: Sat Sep 14 01:39:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 2: Sat Sep 14 01:39:23 2024

    ## Permutation completed: Sat Sep 14 01:39:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 2: Sat Sep 14 01:39:23 2024

    ## Permutation completed: Sat Sep 14 01:39:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 2: Sat Sep 14 01:39:24 2024

    ## Permutation completed: Sat Sep 14 01:39:24 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 2: Sat Sep 14 01:39:24 2024

    ## Permutation completed: Sat Sep 14 01:39:24 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 3: Sat Sep 14 01:39:24 2024

    ## Permutation completed: Sat Sep 14 01:39:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 3: Sat Sep 14 01:39:25 2024

    ## Permutation completed: Sat Sep 14 01:39:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 3: Sat Sep 14 01:39:25 2024

    ## Permutation completed: Sat Sep 14 01:39:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 3: Sat Sep 14 01:39:25 2024

    ## Permutation completed: Sat Sep 14 01:39:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 3: Sat Sep 14 01:39:25 2024

    ## Permutation completed: Sat Sep 14 01:39:25 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 4: Sat Sep 14 01:39:25 2024

    ## Permutation completed: Sat Sep 14 01:39:26 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 4: Sat Sep 14 01:39:26 2024

    ## Permutation completed: Sat Sep 14 01:39:26 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 4: Sat Sep 14 01:39:27 2024

    ## Permutation completed: Sat Sep 14 01:39:27 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 4: Sat Sep 14 01:39:27 2024

    ## Permutation completed: Sat Sep 14 01:39:27 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 4: Sat Sep 14 01:39:27 2024

    ## Permutation completed: Sat Sep 14 01:39:28 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 5: Sat Sep 14 01:39:29 2024

    ## Permutation completed: Sat Sep 14 01:39:29 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 5: Sat Sep 14 01:39:29 2024

    ## Permutation completed: Sat Sep 14 01:39:29 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 5: Sat Sep 14 01:39:29 2024

    ## Permutation completed: Sat Sep 14 01:39:29 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 5: Sat Sep 14 01:39:29 2024

    ## Permutation completed: Sat Sep 14 01:39:30 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 5: Sat Sep 14 01:39:30 2024

    ## Permutation completed: Sat Sep 14 01:39:30 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 6: Sat Sep 14 01:39:30 2024

    ## Permutation completed: Sat Sep 14 01:39:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 6: Sat Sep 14 01:39:31 2024

    ## Permutation completed: Sat Sep 14 01:39:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 6: Sat Sep 14 01:39:31 2024

    ## Permutation completed: Sat Sep 14 01:39:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 6: Sat Sep 14 01:39:31 2024

    ## Permutation completed: Sat Sep 14 01:39:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 6: Sat Sep 14 01:39:31 2024

    ## Permutation completed: Sat Sep 14 01:39:32 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 7: Sat Sep 14 01:39:33 2024

    ## Permutation completed: Sat Sep 14 01:39:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 7: Sat Sep 14 01:39:33 2024

    ## Permutation completed: Sat Sep 14 01:39:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 7: Sat Sep 14 01:39:33 2024

    ## Permutation completed: Sat Sep 14 01:39:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 7: Sat Sep 14 01:39:33 2024

    ## Permutation completed: Sat Sep 14 01:39:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 7: Sat Sep 14 01:39:34 2024

    ## Permutation completed: Sat Sep 14 01:39:34 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 8: Sat Sep 14 01:39:34 2024

    ## Permutation completed: Sat Sep 14 01:39:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 8: Sat Sep 14 01:39:35 2024

    ## Permutation completed: Sat Sep 14 01:39:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 8: Sat Sep 14 01:39:35 2024

    ## Permutation completed: Sat Sep 14 01:39:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 8: Sat Sep 14 01:39:35 2024

    ## Permutation completed: Sat Sep 14 01:39:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 8: Sat Sep 14 01:39:35 2024

    ## Permutation completed: Sat Sep 14 01:39:35 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 9: Sat Sep 14 01:39:35 2024

    ## Permutation completed: Sat Sep 14 01:39:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 9: Sat Sep 14 01:39:35 2024

    ## Permutation completed: Sat Sep 14 01:39:36 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 9: Sat Sep 14 01:39:36 2024

    ## Permutation completed: Sat Sep 14 01:39:36 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 9: Sat Sep 14 01:39:36 2024

    ## Permutation completed: Sat Sep 14 01:39:36 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 9: Sat Sep 14 01:39:36 2024

    ## Permutation completed: Sat Sep 14 01:39:36 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 10: Sat Sep 14 01:39:36 2024

    ## Permutation completed: Sat Sep 14 01:39:37 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 10: Sat Sep 14 01:39:37 2024

    ## Permutation completed: Sat Sep 14 01:39:37 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 10: Sat Sep 14 01:39:37 2024

    ## Permutation completed: Sat Sep 14 01:39:37 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 10: Sat Sep 14 01:39:37 2024

    ## Permutation completed: Sat Sep 14 01:39:37 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 10: Sat Sep 14 01:39:37 2024

    ## Permutation completed: Sat Sep 14 01:39:37 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 11: Sat Sep 14 01:39:38 2024

    ## Permutation completed: Sat Sep 14 01:39:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 11: Sat Sep 14 01:39:38 2024

    ## Permutation completed: Sat Sep 14 01:39:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 11: Sat Sep 14 01:39:38 2024

    ## Permutation completed: Sat Sep 14 01:39:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 11: Sat Sep 14 01:39:38 2024

    ## Permutation completed: Sat Sep 14 01:39:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 11: Sat Sep 14 01:39:38 2024

    ## Permutation completed: Sat Sep 14 01:39:39 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 12: Sat Sep 14 01:39:39 2024

    ## Permutation completed: Sat Sep 14 01:39:39 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 12: Sat Sep 14 01:39:40 2024

    ## Permutation completed: Sat Sep 14 01:39:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 12: Sat Sep 14 01:39:40 2024

    ## Permutation completed: Sat Sep 14 01:39:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 12: Sat Sep 14 01:39:40 2024

    ## Permutation completed: Sat Sep 14 01:39:41 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 12: Sat Sep 14 01:39:41 2024

    ## Permutation completed: Sat Sep 14 01:39:41 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 13: Sat Sep 14 01:39:42 2024

    ## Permutation completed: Sat Sep 14 01:39:42 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 13: Sat Sep 14 01:39:42 2024

    ## Permutation completed: Sat Sep 14 01:39:42 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 13: Sat Sep 14 01:39:42 2024

    ## Permutation completed: Sat Sep 14 01:39:43 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 13: Sat Sep 14 01:39:43 2024

    ## Permutation completed: Sat Sep 14 01:39:43 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 13: Sat Sep 14 01:39:43 2024

    ## Permutation completed: Sat Sep 14 01:39:43 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 14: Sat Sep 14 01:39:43 2024

    ## Permutation completed: Sat Sep 14 01:39:43 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 14: Sat Sep 14 01:39:43 2024

    ## Permutation completed: Sat Sep 14 01:39:43 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 14: Sat Sep 14 01:39:44 2024

    ## Permutation completed: Sat Sep 14 01:39:44 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 14: Sat Sep 14 01:39:44 2024

    ## Permutation completed: Sat Sep 14 01:39:44 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 14: Sat Sep 14 01:39:44 2024

    ## Permutation completed: Sat Sep 14 01:39:44 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 15: Sat Sep 14 01:39:44 2024

    ## Permutation completed: Sat Sep 14 01:39:44 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 15: Sat Sep 14 01:39:44 2024

    ## Permutation completed: Sat Sep 14 01:39:45 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 15: Sat Sep 14 01:39:45 2024

    ## Permutation completed: Sat Sep 14 01:39:45 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 15: Sat Sep 14 01:39:45 2024

    ## Permutation completed: Sat Sep 14 01:39:45 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 15: Sat Sep 14 01:39:45 2024

    ## Permutation completed: Sat Sep 14 01:39:45 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 16: Sat Sep 14 01:39:46 2024

    ## Permutation completed: Sat Sep 14 01:39:46 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 16: Sat Sep 14 01:39:46 2024

    ## Permutation completed: Sat Sep 14 01:39:46 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 16: Sat Sep 14 01:39:46 2024

    ## Permutation completed: Sat Sep 14 01:39:46 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 16: Sat Sep 14 01:39:46 2024

    ## Permutation completed: Sat Sep 14 01:39:46 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 16: Sat Sep 14 01:39:46 2024

    ## Permutation completed: Sat Sep 14 01:39:46 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 17: Sat Sep 14 01:39:47 2024

    ## Permutation completed: Sat Sep 14 01:39:47 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 17: Sat Sep 14 01:39:47 2024

    ## Permutation completed: Sat Sep 14 01:39:47 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 17: Sat Sep 14 01:39:47 2024

    ## Permutation completed: Sat Sep 14 01:39:47 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 17: Sat Sep 14 01:39:47 2024

    ## Permutation completed: Sat Sep 14 01:39:47 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 17: Sat Sep 14 01:39:47 2024

    ## Permutation completed: Sat Sep 14 01:39:47 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 18: Sat Sep 14 01:39:48 2024

    ## Permutation completed: Sat Sep 14 01:39:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 18: Sat Sep 14 01:39:48 2024

    ## Permutation completed: Sat Sep 14 01:39:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 18: Sat Sep 14 01:39:48 2024

    ## Permutation completed: Sat Sep 14 01:39:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 18: Sat Sep 14 01:39:48 2024

    ## Permutation completed: Sat Sep 14 01:39:49 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 18: Sat Sep 14 01:39:49 2024

    ## Permutation completed: Sat Sep 14 01:39:49 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 19: Sat Sep 14 01:39:49 2024

    ## Permutation completed: Sat Sep 14 01:39:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 19: Sat Sep 14 01:39:50 2024

    ## Permutation completed: Sat Sep 14 01:39:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 19: Sat Sep 14 01:39:50 2024

    ## Permutation completed: Sat Sep 14 01:39:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 19: Sat Sep 14 01:39:50 2024

    ## Permutation completed: Sat Sep 14 01:39:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 19: Sat Sep 14 01:39:50 2024

    ## Permutation completed: Sat Sep 14 01:39:50 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 20: Sat Sep 14 01:39:51 2024

    ## Permutation completed: Sat Sep 14 01:39:51 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 20: Sat Sep 14 01:39:51 2024

    ## Permutation completed: Sat Sep 14 01:39:51 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 20: Sat Sep 14 01:39:51 2024

    ## Permutation completed: Sat Sep 14 01:39:51 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 20: Sat Sep 14 01:39:51 2024

    ## Permutation completed: Sat Sep 14 01:39:51 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 20: Sat Sep 14 01:39:51 2024

    ## Permutation completed: Sat Sep 14 01:39:51 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 21: Sat Sep 14 01:39:52 2024

    ## Permutation completed: Sat Sep 14 01:39:52 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 21: Sat Sep 14 01:39:52 2024

    ## Permutation completed: Sat Sep 14 01:39:52 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 21: Sat Sep 14 01:39:52 2024

    ## Permutation completed: Sat Sep 14 01:39:53 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 21: Sat Sep 14 01:39:53 2024

    ## Permutation completed: Sat Sep 14 01:39:53 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 21: Sat Sep 14 01:39:53 2024

    ## Permutation completed: Sat Sep 14 01:39:53 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 22: Sat Sep 14 01:39:54 2024

    ## Permutation completed: Sat Sep 14 01:39:54 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 22: Sat Sep 14 01:39:54 2024

    ## Permutation completed: Sat Sep 14 01:39:54 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 22: Sat Sep 14 01:39:54 2024

    ## Permutation completed: Sat Sep 14 01:39:55 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 22: Sat Sep 14 01:39:55 2024

    ## Permutation completed: Sat Sep 14 01:39:55 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 22: Sat Sep 14 01:39:55 2024

    ## Permutation completed: Sat Sep 14 01:39:55 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 23: Sat Sep 14 01:39:55 2024

    ## Permutation completed: Sat Sep 14 01:39:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 23: Sat Sep 14 01:39:56 2024

    ## Permutation completed: Sat Sep 14 01:39:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 23: Sat Sep 14 01:39:56 2024

    ## Permutation completed: Sat Sep 14 01:39:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 23: Sat Sep 14 01:39:56 2024

    ## Permutation completed: Sat Sep 14 01:39:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 23: Sat Sep 14 01:39:56 2024

    ## Permutation completed: Sat Sep 14 01:39:56 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 24: Sat Sep 14 01:39:57 2024

    ## Permutation completed: Sat Sep 14 01:39:57 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 24: Sat Sep 14 01:39:57 2024

    ## Permutation completed: Sat Sep 14 01:39:57 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 24: Sat Sep 14 01:39:57 2024

    ## Permutation completed: Sat Sep 14 01:39:57 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 24: Sat Sep 14 01:39:57 2024

    ## Permutation completed: Sat Sep 14 01:39:57 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 24: Sat Sep 14 01:39:57 2024

    ## Permutation completed: Sat Sep 14 01:39:57 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 25: Sat Sep 14 01:39:58 2024

    ## Permutation completed: Sat Sep 14 01:39:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 25: Sat Sep 14 01:39:58 2024

    ## Permutation completed: Sat Sep 14 01:39:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 25: Sat Sep 14 01:39:58 2024

    ## Permutation completed: Sat Sep 14 01:39:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 25: Sat Sep 14 01:39:58 2024

    ## Permutation completed: Sat Sep 14 01:39:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 25: Sat Sep 14 01:39:58 2024

    ## Permutation completed: Sat Sep 14 01:39:58 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 26: Sat Sep 14 01:39:59 2024

    ## Permutation completed: Sat Sep 14 01:39:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 26: Sat Sep 14 01:39:59 2024

    ## Permutation completed: Sat Sep 14 01:39:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 26: Sat Sep 14 01:40:00 2024

    ## Permutation completed: Sat Sep 14 01:40:00 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 26: Sat Sep 14 01:40:00 2024

    ## Permutation completed: Sat Sep 14 01:40:00 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 26: Sat Sep 14 01:40:00 2024

    ## Permutation completed: Sat Sep 14 01:40:00 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 27: Sat Sep 14 01:40:01 2024

    ## Permutation completed: Sat Sep 14 01:40:01 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 27: Sat Sep 14 01:40:01 2024

    ## Permutation completed: Sat Sep 14 01:40:01 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 27: Sat Sep 14 01:40:01 2024

    ## Permutation completed: Sat Sep 14 01:40:02 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 27: Sat Sep 14 01:40:02 2024

    ## Permutation completed: Sat Sep 14 01:40:02 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 27: Sat Sep 14 01:40:02 2024

    ## Permutation completed: Sat Sep 14 01:40:02 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 28: Sat Sep 14 01:40:02 2024

    ## Permutation completed: Sat Sep 14 01:40:02 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 28: Sat Sep 14 01:40:02 2024

    ## Permutation completed: Sat Sep 14 01:40:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 28: Sat Sep 14 01:40:03 2024

    ## Permutation completed: Sat Sep 14 01:40:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 28: Sat Sep 14 01:40:03 2024

    ## Permutation completed: Sat Sep 14 01:40:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 28: Sat Sep 14 01:40:03 2024

    ## Permutation completed: Sat Sep 14 01:40:03 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 29: Sat Sep 14 01:40:03 2024

    ## Permutation completed: Sat Sep 14 01:40:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 29: Sat Sep 14 01:40:03 2024

    ## Permutation completed: Sat Sep 14 01:40:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 29: Sat Sep 14 01:40:03 2024

    ## Permutation completed: Sat Sep 14 01:40:04 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 29: Sat Sep 14 01:40:04 2024

    ## Permutation completed: Sat Sep 14 01:40:04 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 29: Sat Sep 14 01:40:04 2024

    ## Permutation completed: Sat Sep 14 01:40:04 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 30: Sat Sep 14 01:40:04 2024

    ## Permutation completed: Sat Sep 14 01:40:04 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 30: Sat Sep 14 01:40:04 2024

    ## Permutation completed: Sat Sep 14 01:40:05 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 30: Sat Sep 14 01:40:05 2024

    ## Permutation completed: Sat Sep 14 01:40:05 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 30: Sat Sep 14 01:40:05 2024

    ## Permutation completed: Sat Sep 14 01:40:05 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 30: Sat Sep 14 01:40:05 2024

    ## Permutation completed: Sat Sep 14 01:40:05 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 31: Sat Sep 14 01:40:05 2024

    ## Permutation completed: Sat Sep 14 01:40:06 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 31: Sat Sep 14 01:40:06 2024

    ## Permutation completed: Sat Sep 14 01:40:06 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 31: Sat Sep 14 01:40:06 2024

    ## Permutation completed: Sat Sep 14 01:40:06 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 31: Sat Sep 14 01:40:06 2024

    ## Permutation completed: Sat Sep 14 01:40:06 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 31: Sat Sep 14 01:40:06 2024

    ## Permutation completed: Sat Sep 14 01:40:06 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 32: Sat Sep 14 01:40:06 2024

    ## Permutation completed: Sat Sep 14 01:40:07 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 32: Sat Sep 14 01:40:07 2024

    ## Permutation completed: Sat Sep 14 01:40:07 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 32: Sat Sep 14 01:40:07 2024

    ## Permutation completed: Sat Sep 14 01:40:07 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 32: Sat Sep 14 01:40:07 2024

    ## Permutation completed: Sat Sep 14 01:40:07 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 32: Sat Sep 14 01:40:07 2024

    ## Permutation completed: Sat Sep 14 01:40:07 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 33: Sat Sep 14 01:40:07 2024

    ## Permutation completed: Sat Sep 14 01:40:08 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 33: Sat Sep 14 01:40:08 2024

    ## Permutation completed: Sat Sep 14 01:40:08 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 33: Sat Sep 14 01:40:08 2024

    ## Permutation completed: Sat Sep 14 01:40:08 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 33: Sat Sep 14 01:40:08 2024

    ## Permutation completed: Sat Sep 14 01:40:08 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 33: Sat Sep 14 01:40:08 2024

    ## Permutation completed: Sat Sep 14 01:40:08 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 34: Sat Sep 14 01:40:09 2024

    ## Permutation completed: Sat Sep 14 01:40:09 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 34: Sat Sep 14 01:40:09 2024

    ## Permutation completed: Sat Sep 14 01:40:09 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 34: Sat Sep 14 01:40:09 2024

    ## Permutation completed: Sat Sep 14 01:40:09 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 34: Sat Sep 14 01:40:09 2024

    ## Permutation completed: Sat Sep 14 01:40:09 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 34: Sat Sep 14 01:40:09 2024

    ## Permutation completed: Sat Sep 14 01:40:09 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 35: Sat Sep 14 01:40:10 2024

    ## Permutation completed: Sat Sep 14 01:40:10 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 35: Sat Sep 14 01:40:10 2024

    ## Permutation completed: Sat Sep 14 01:40:10 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 35: Sat Sep 14 01:40:10 2024

    ## Permutation completed: Sat Sep 14 01:40:11 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 35: Sat Sep 14 01:40:11 2024

    ## Permutation completed: Sat Sep 14 01:40:11 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 35: Sat Sep 14 01:40:11 2024

    ## Permutation completed: Sat Sep 14 01:40:11 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 36: Sat Sep 14 01:40:11 2024

    ## Permutation completed: Sat Sep 14 01:40:12 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 36: Sat Sep 14 01:40:12 2024

    ## Permutation completed: Sat Sep 14 01:40:12 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 36: Sat Sep 14 01:40:12 2024

    ## Permutation completed: Sat Sep 14 01:40:12 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 36: Sat Sep 14 01:40:12 2024

    ## Permutation completed: Sat Sep 14 01:40:13 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 36: Sat Sep 14 01:40:13 2024

    ## Permutation completed: Sat Sep 14 01:40:13 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 37: Sat Sep 14 01:40:14 2024

    ## Permutation completed: Sat Sep 14 01:40:15 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 37: Sat Sep 14 01:40:15 2024

    ## Permutation completed: Sat Sep 14 01:40:15 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 37: Sat Sep 14 01:40:15 2024

    ## Permutation completed: Sat Sep 14 01:40:15 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 37: Sat Sep 14 01:40:15 2024

    ## Permutation completed: Sat Sep 14 01:40:15 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 37: Sat Sep 14 01:40:15 2024

    ## Permutation completed: Sat Sep 14 01:40:15 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 38: Sat Sep 14 01:40:16 2024

    ## Permutation completed: Sat Sep 14 01:40:16 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 38: Sat Sep 14 01:40:16 2024

    ## Permutation completed: Sat Sep 14 01:40:16 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 38: Sat Sep 14 01:40:16 2024

    ## Permutation completed: Sat Sep 14 01:40:16 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 38: Sat Sep 14 01:40:16 2024

    ## Permutation completed: Sat Sep 14 01:40:17 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 38: Sat Sep 14 01:40:17 2024

    ## Permutation completed: Sat Sep 14 01:40:17 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 39: Sat Sep 14 01:40:17 2024

    ## Permutation completed: Sat Sep 14 01:40:17 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 39: Sat Sep 14 01:40:17 2024

    ## Permutation completed: Sat Sep 14 01:40:17 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 39: Sat Sep 14 01:40:17 2024

    ## Permutation completed: Sat Sep 14 01:40:18 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 39: Sat Sep 14 01:40:18 2024

    ## Permutation completed: Sat Sep 14 01:40:18 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 39: Sat Sep 14 01:40:18 2024

    ## Permutation completed: Sat Sep 14 01:40:18 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 40: Sat Sep 14 01:40:18 2024

    ## Permutation completed: Sat Sep 14 01:40:18 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 40: Sat Sep 14 01:40:18 2024

    ## Permutation completed: Sat Sep 14 01:40:19 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 40: Sat Sep 14 01:40:19 2024

    ## Permutation completed: Sat Sep 14 01:40:19 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 40: Sat Sep 14 01:40:19 2024

    ## Permutation completed: Sat Sep 14 01:40:19 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 40: Sat Sep 14 01:40:19 2024

    ## Permutation completed: Sat Sep 14 01:40:19 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 41: Sat Sep 14 01:40:19 2024

    ## Permutation completed: Sat Sep 14 01:40:20 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 41: Sat Sep 14 01:40:20 2024

    ## Permutation completed: Sat Sep 14 01:40:20 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 41: Sat Sep 14 01:40:20 2024

    ## Permutation completed: Sat Sep 14 01:40:20 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 41: Sat Sep 14 01:40:21 2024

    ## Permutation completed: Sat Sep 14 01:40:21 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 41: Sat Sep 14 01:40:21 2024

    ## Permutation completed: Sat Sep 14 01:40:21 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 42: Sat Sep 14 01:40:23 2024

    ## Permutation completed: Sat Sep 14 01:40:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 42: Sat Sep 14 01:40:23 2024

    ## Permutation completed: Sat Sep 14 01:40:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 42: Sat Sep 14 01:40:23 2024

    ## Permutation completed: Sat Sep 14 01:40:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 42: Sat Sep 14 01:40:23 2024

    ## Permutation completed: Sat Sep 14 01:40:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 42: Sat Sep 14 01:40:23 2024

    ## Permutation completed: Sat Sep 14 01:40:23 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 43: Sat Sep 14 01:40:24 2024

    ## Permutation completed: Sat Sep 14 01:40:24 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 43: Sat Sep 14 01:40:24 2024

    ## Permutation completed: Sat Sep 14 01:40:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 43: Sat Sep 14 01:40:25 2024

    ## Permutation completed: Sat Sep 14 01:40:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 43: Sat Sep 14 01:40:25 2024

    ## Permutation completed: Sat Sep 14 01:40:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 43: Sat Sep 14 01:40:25 2024

    ## Permutation completed: Sat Sep 14 01:40:25 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 44: Sat Sep 14 01:40:26 2024

    ## Permutation completed: Sat Sep 14 01:40:26 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 44: Sat Sep 14 01:40:26 2024

    ## Permutation completed: Sat Sep 14 01:40:27 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 44: Sat Sep 14 01:40:27 2024

    ## Permutation completed: Sat Sep 14 01:40:27 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 44: Sat Sep 14 01:40:27 2024

    ## Permutation completed: Sat Sep 14 01:40:27 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 44: Sat Sep 14 01:40:27 2024

    ## Permutation completed: Sat Sep 14 01:40:27 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 45: Sat Sep 14 01:40:27 2024

    ## Permutation completed: Sat Sep 14 01:40:27 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 45: Sat Sep 14 01:40:28 2024

    ## Permutation completed: Sat Sep 14 01:40:28 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 45: Sat Sep 14 01:40:28 2024

    ## Permutation completed: Sat Sep 14 01:40:28 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 45: Sat Sep 14 01:40:28 2024

    ## Permutation completed: Sat Sep 14 01:40:28 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 45: Sat Sep 14 01:40:28 2024

    ## Permutation completed: Sat Sep 14 01:40:28 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 46: Sat Sep 14 01:40:29 2024

    ## Permutation completed: Sat Sep 14 01:40:29 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 46: Sat Sep 14 01:40:29 2024

    ## Permutation completed: Sat Sep 14 01:40:29 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 46: Sat Sep 14 01:40:29 2024

    ## Permutation completed: Sat Sep 14 01:40:30 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 46: Sat Sep 14 01:40:30 2024

    ## Permutation completed: Sat Sep 14 01:40:30 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 46: Sat Sep 14 01:40:30 2024

    ## Permutation completed: Sat Sep 14 01:40:30 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 47: Sat Sep 14 01:40:31 2024

    ## Permutation completed: Sat Sep 14 01:40:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 47: Sat Sep 14 01:40:31 2024

    ## Permutation completed: Sat Sep 14 01:40:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 47: Sat Sep 14 01:40:31 2024

    ## Permutation completed: Sat Sep 14 01:40:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 47: Sat Sep 14 01:40:31 2024

    ## Permutation completed: Sat Sep 14 01:40:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 47: Sat Sep 14 01:40:31 2024

    ## Permutation completed: Sat Sep 14 01:40:32 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 48: Sat Sep 14 01:40:33 2024

    ## Permutation completed: Sat Sep 14 01:40:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 48: Sat Sep 14 01:40:33 2024

    ## Permutation completed: Sat Sep 14 01:40:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 48: Sat Sep 14 01:40:33 2024

    ## Permutation completed: Sat Sep 14 01:40:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 48: Sat Sep 14 01:40:33 2024

    ## Permutation completed: Sat Sep 14 01:40:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 48: Sat Sep 14 01:40:34 2024

    ## Permutation completed: Sat Sep 14 01:40:34 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 49: Sat Sep 14 01:40:34 2024

    ## Permutation completed: Sat Sep 14 01:40:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 49: Sat Sep 14 01:40:34 2024

    ## Permutation completed: Sat Sep 14 01:40:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 49: Sat Sep 14 01:40:34 2024

    ## Permutation completed: Sat Sep 14 01:40:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 49: Sat Sep 14 01:40:34 2024

    ## Permutation completed: Sat Sep 14 01:40:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 49: Sat Sep 14 01:40:35 2024

    ## Permutation completed: Sat Sep 14 01:40:35 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 50: Sat Sep 14 01:40:35 2024

    ## Permutation completed: Sat Sep 14 01:40:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 50: Sat Sep 14 01:40:35 2024

    ## Permutation completed: Sat Sep 14 01:40:36 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 50: Sat Sep 14 01:40:36 2024

    ## Permutation completed: Sat Sep 14 01:40:36 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 50: Sat Sep 14 01:40:36 2024

    ## Permutation completed: Sat Sep 14 01:40:36 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 50: Sat Sep 14 01:40:36 2024

    ## Permutation completed: Sat Sep 14 01:40:36 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 51: Sat Sep 14 01:40:37 2024

    ## Permutation completed: Sat Sep 14 01:40:37 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 51: Sat Sep 14 01:40:37 2024

    ## Permutation completed: Sat Sep 14 01:40:37 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 51: Sat Sep 14 01:40:37 2024

    ## Permutation completed: Sat Sep 14 01:40:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 51: Sat Sep 14 01:40:38 2024

    ## Permutation completed: Sat Sep 14 01:40:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 51: Sat Sep 14 01:40:38 2024

    ## Permutation completed: Sat Sep 14 01:40:38 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 52: Sat Sep 14 01:40:38 2024

    ## Permutation completed: Sat Sep 14 01:40:39 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 52: Sat Sep 14 01:40:39 2024

    ## Permutation completed: Sat Sep 14 01:40:39 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 52: Sat Sep 14 01:40:39 2024

    ## Permutation completed: Sat Sep 14 01:40:39 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 52: Sat Sep 14 01:40:39 2024

    ## Permutation completed: Sat Sep 14 01:40:39 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 52: Sat Sep 14 01:40:39 2024

    ## Permutation completed: Sat Sep 14 01:40:39 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 53: Sat Sep 14 01:40:39 2024

    ## Permutation completed: Sat Sep 14 01:40:39 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 53: Sat Sep 14 01:40:40 2024

    ## Permutation completed: Sat Sep 14 01:40:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 53: Sat Sep 14 01:40:40 2024

    ## Permutation completed: Sat Sep 14 01:40:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 53: Sat Sep 14 01:40:40 2024

    ## Permutation completed: Sat Sep 14 01:40:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 53: Sat Sep 14 01:40:40 2024

    ## Permutation completed: Sat Sep 14 01:40:40 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 54: Sat Sep 14 01:40:41 2024

    ## Permutation completed: Sat Sep 14 01:40:41 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 54: Sat Sep 14 01:40:41 2024

    ## Permutation completed: Sat Sep 14 01:40:41 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 54: Sat Sep 14 01:40:41 2024

    ## Permutation completed: Sat Sep 14 01:40:41 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 54: Sat Sep 14 01:40:41 2024

    ## Permutation completed: Sat Sep 14 01:40:41 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 54: Sat Sep 14 01:40:41 2024

    ## Permutation completed: Sat Sep 14 01:40:41 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 55: Sat Sep 14 01:40:42 2024

    ## Permutation completed: Sat Sep 14 01:40:42 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 55: Sat Sep 14 01:40:42 2024

    ## Permutation completed: Sat Sep 14 01:40:42 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 55: Sat Sep 14 01:40:42 2024

    ## Permutation completed: Sat Sep 14 01:40:42 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 55: Sat Sep 14 01:40:42 2024

    ## Permutation completed: Sat Sep 14 01:40:42 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 55: Sat Sep 14 01:40:42 2024

    ## Permutation completed: Sat Sep 14 01:40:43 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 56: Sat Sep 14 01:40:43 2024

    ## Permutation completed: Sat Sep 14 01:40:44 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 56: Sat Sep 14 01:40:44 2024

    ## Permutation completed: Sat Sep 14 01:40:45 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 56: Sat Sep 14 01:40:45 2024

    ## Permutation completed: Sat Sep 14 01:40:45 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 56: Sat Sep 14 01:40:45 2024

    ## Permutation completed: Sat Sep 14 01:40:46 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 56: Sat Sep 14 01:40:46 2024

    ## Permutation completed: Sat Sep 14 01:40:46 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 57: Sat Sep 14 01:40:48 2024

    ## Permutation completed: Sat Sep 14 01:40:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 57: Sat Sep 14 01:40:48 2024

    ## Permutation completed: Sat Sep 14 01:40:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 57: Sat Sep 14 01:40:48 2024

    ## Permutation completed: Sat Sep 14 01:40:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 57: Sat Sep 14 01:40:48 2024

    ## Permutation completed: Sat Sep 14 01:40:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 57: Sat Sep 14 01:40:48 2024

    ## Permutation completed: Sat Sep 14 01:40:48 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 58: Sat Sep 14 01:40:49 2024

    ## Permutation completed: Sat Sep 14 01:40:49 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 58: Sat Sep 14 01:40:49 2024

    ## Permutation completed: Sat Sep 14 01:40:49 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 58: Sat Sep 14 01:40:49 2024

    ## Permutation completed: Sat Sep 14 01:40:49 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 58: Sat Sep 14 01:40:49 2024

    ## Permutation completed: Sat Sep 14 01:40:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 58: Sat Sep 14 01:40:50 2024

    ## Permutation completed: Sat Sep 14 01:40:50 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 59: Sat Sep 14 01:40:50 2024

    ## Permutation completed: Sat Sep 14 01:40:51 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 59: Sat Sep 14 01:40:51 2024

    ## Permutation completed: Sat Sep 14 01:40:51 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 59: Sat Sep 14 01:40:51 2024

    ## Permutation completed: Sat Sep 14 01:40:51 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 59: Sat Sep 14 01:40:51 2024

    ## Permutation completed: Sat Sep 14 01:40:52 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 59: Sat Sep 14 01:40:52 2024

    ## Permutation completed: Sat Sep 14 01:40:52 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 60: Sat Sep 14 01:40:53 2024

    ## Permutation completed: Sat Sep 14 01:40:53 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 60: Sat Sep 14 01:40:53 2024

    ## Permutation completed: Sat Sep 14 01:40:53 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 60: Sat Sep 14 01:40:53 2024

    ## Permutation completed: Sat Sep 14 01:40:54 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 60: Sat Sep 14 01:40:54 2024

    ## Permutation completed: Sat Sep 14 01:40:54 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 60: Sat Sep 14 01:40:54 2024

    ## Permutation completed: Sat Sep 14 01:40:54 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 61: Sat Sep 14 01:40:55 2024

    ## Permutation completed: Sat Sep 14 01:40:55 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 61: Sat Sep 14 01:40:55 2024

    ## Permutation completed: Sat Sep 14 01:40:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 61: Sat Sep 14 01:40:56 2024

    ## Permutation completed: Sat Sep 14 01:40:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 61: Sat Sep 14 01:40:56 2024

    ## Permutation completed: Sat Sep 14 01:40:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 61: Sat Sep 14 01:40:56 2024

    ## Permutation completed: Sat Sep 14 01:40:56 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 62: Sat Sep 14 01:40:57 2024

    ## Permutation completed: Sat Sep 14 01:40:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 62: Sat Sep 14 01:40:58 2024

    ## Permutation completed: Sat Sep 14 01:40:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 62: Sat Sep 14 01:40:58 2024

    ## Permutation completed: Sat Sep 14 01:40:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 62: Sat Sep 14 01:40:58 2024

    ## Permutation completed: Sat Sep 14 01:40:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 62: Sat Sep 14 01:40:58 2024

    ## Permutation completed: Sat Sep 14 01:40:58 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 63: Sat Sep 14 01:40:59 2024

    ## Permutation completed: Sat Sep 14 01:40:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 63: Sat Sep 14 01:40:59 2024

    ## Permutation completed: Sat Sep 14 01:40:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 63: Sat Sep 14 01:40:59 2024

    ## Permutation completed: Sat Sep 14 01:40:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 63: Sat Sep 14 01:40:59 2024

    ## Permutation completed: Sat Sep 14 01:40:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 63: Sat Sep 14 01:40:59 2024

    ## Permutation completed: Sat Sep 14 01:41:00 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 64: Sat Sep 14 01:41:00 2024

    ## Permutation completed: Sat Sep 14 01:41:01 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 64: Sat Sep 14 01:41:01 2024

    ## Permutation completed: Sat Sep 14 01:41:02 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 64: Sat Sep 14 01:41:02 2024

    ## Permutation completed: Sat Sep 14 01:41:02 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 64: Sat Sep 14 01:41:02 2024

    ## Permutation completed: Sat Sep 14 01:41:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 64: Sat Sep 14 01:41:03 2024

    ## Permutation completed: Sat Sep 14 01:41:04 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 65: Sat Sep 14 01:41:05 2024

    ## Permutation completed: Sat Sep 14 01:41:05 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 65: Sat Sep 14 01:41:05 2024

    ## Permutation completed: Sat Sep 14 01:41:05 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 65: Sat Sep 14 01:41:05 2024

    ## Permutation completed: Sat Sep 14 01:41:05 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 65: Sat Sep 14 01:41:05 2024

    ## Permutation completed: Sat Sep 14 01:41:06 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 65: Sat Sep 14 01:41:06 2024

    ## Permutation completed: Sat Sep 14 01:41:06 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 66: Sat Sep 14 01:41:06 2024

    ## Permutation completed: Sat Sep 14 01:41:07 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 66: Sat Sep 14 01:41:07 2024

    ## Permutation completed: Sat Sep 14 01:41:07 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 66: Sat Sep 14 01:41:07 2024

    ## Permutation completed: Sat Sep 14 01:41:08 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 66: Sat Sep 14 01:41:08 2024

    ## Permutation completed: Sat Sep 14 01:41:08 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 66: Sat Sep 14 01:41:08 2024

    ## Permutation completed: Sat Sep 14 01:41:09 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 67: Sat Sep 14 01:41:10 2024

    ## Permutation completed: Sat Sep 14 01:41:10 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 67: Sat Sep 14 01:41:10 2024

    ## Permutation completed: Sat Sep 14 01:41:11 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 67: Sat Sep 14 01:41:11 2024

    ## Permutation completed: Sat Sep 14 01:41:11 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 67: Sat Sep 14 01:41:11 2024

    ## Permutation completed: Sat Sep 14 01:41:11 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 67: Sat Sep 14 01:41:11 2024

    ## Permutation completed: Sat Sep 14 01:41:11 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 68: Sat Sep 14 01:41:12 2024

    ## Permutation completed: Sat Sep 14 01:41:12 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 68: Sat Sep 14 01:41:12 2024

    ## Permutation completed: Sat Sep 14 01:41:13 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 68: Sat Sep 14 01:41:13 2024

    ## Permutation completed: Sat Sep 14 01:41:13 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 68: Sat Sep 14 01:41:13 2024

    ## Permutation completed: Sat Sep 14 01:41:13 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 68: Sat Sep 14 01:41:13 2024

    ## Permutation completed: Sat Sep 14 01:41:13 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 69: Sat Sep 14 01:41:14 2024

    ## Permutation completed: Sat Sep 14 01:41:14 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 69: Sat Sep 14 01:41:14 2024

    ## Permutation completed: Sat Sep 14 01:41:14 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 69: Sat Sep 14 01:41:15 2024

    ## Permutation completed: Sat Sep 14 01:41:15 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 69: Sat Sep 14 01:41:15 2024

    ## Permutation completed: Sat Sep 14 01:41:15 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 69: Sat Sep 14 01:41:15 2024

    ## Permutation completed: Sat Sep 14 01:41:15 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 70: Sat Sep 14 01:41:16 2024

    ## Permutation completed: Sat Sep 14 01:41:16 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 70: Sat Sep 14 01:41:16 2024

    ## Permutation completed: Sat Sep 14 01:41:16 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 70: Sat Sep 14 01:41:16 2024

    ## Permutation completed: Sat Sep 14 01:41:17 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 70: Sat Sep 14 01:41:17 2024

    ## Permutation completed: Sat Sep 14 01:41:17 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 70: Sat Sep 14 01:41:17 2024

    ## Permutation completed: Sat Sep 14 01:41:17 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 71: Sat Sep 14 01:41:19 2024

    ## Permutation completed: Sat Sep 14 01:41:19 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 71: Sat Sep 14 01:41:19 2024

    ## Permutation completed: Sat Sep 14 01:41:19 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 71: Sat Sep 14 01:41:19 2024

    ## Permutation completed: Sat Sep 14 01:41:19 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 71: Sat Sep 14 01:41:19 2024

    ## Permutation completed: Sat Sep 14 01:41:19 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 71: Sat Sep 14 01:41:19 2024

    ## Permutation completed: Sat Sep 14 01:41:19 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 72: Sat Sep 14 01:41:20 2024

    ## Permutation completed: Sat Sep 14 01:41:20 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 72: Sat Sep 14 01:41:20 2024

    ## Permutation completed: Sat Sep 14 01:41:20 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 72: Sat Sep 14 01:41:20 2024

    ## Permutation completed: Sat Sep 14 01:41:20 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 72: Sat Sep 14 01:41:20 2024

    ## Permutation completed: Sat Sep 14 01:41:20 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 72: Sat Sep 14 01:41:21 2024

    ## Permutation completed: Sat Sep 14 01:41:21 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 73: Sat Sep 14 01:41:21 2024

    ## Permutation completed: Sat Sep 14 01:41:21 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 73: Sat Sep 14 01:41:21 2024

    ## Permutation completed: Sat Sep 14 01:41:21 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 73: Sat Sep 14 01:41:21 2024

    ## Permutation completed: Sat Sep 14 01:41:22 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 73: Sat Sep 14 01:41:22 2024

    ## Permutation completed: Sat Sep 14 01:41:22 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 73: Sat Sep 14 01:41:22 2024

    ## Permutation completed: Sat Sep 14 01:41:22 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 74: Sat Sep 14 01:41:23 2024

    ## Permutation completed: Sat Sep 14 01:41:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 74: Sat Sep 14 01:41:23 2024

    ## Permutation completed: Sat Sep 14 01:41:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 74: Sat Sep 14 01:41:23 2024

    ## Permutation completed: Sat Sep 14 01:41:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 74: Sat Sep 14 01:41:23 2024

    ## Permutation completed: Sat Sep 14 01:41:23 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 74: Sat Sep 14 01:41:23 2024

    ## Permutation completed: Sat Sep 14 01:41:24 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 75: Sat Sep 14 01:41:24 2024

    ## Permutation completed: Sat Sep 14 01:41:24 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 75: Sat Sep 14 01:41:24 2024

    ## Permutation completed: Sat Sep 14 01:41:24 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 75: Sat Sep 14 01:41:25 2024

    ## Permutation completed: Sat Sep 14 01:41:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 75: Sat Sep 14 01:41:25 2024

    ## Permutation completed: Sat Sep 14 01:41:25 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 75: Sat Sep 14 01:41:25 2024

    ## Permutation completed: Sat Sep 14 01:41:25 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 76: Sat Sep 14 01:41:26 2024

    ## Permutation completed: Sat Sep 14 01:41:26 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 76: Sat Sep 14 01:41:26 2024

    ## Permutation completed: Sat Sep 14 01:41:26 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 76: Sat Sep 14 01:41:26 2024

    ## Permutation completed: Sat Sep 14 01:41:26 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 76: Sat Sep 14 01:41:26 2024

    ## Permutation completed: Sat Sep 14 01:41:26 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 76: Sat Sep 14 01:41:26 2024

    ## Permutation completed: Sat Sep 14 01:41:27 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 77: Sat Sep 14 01:41:27 2024

    ## Permutation completed: Sat Sep 14 01:41:27 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 77: Sat Sep 14 01:41:27 2024

    ## Permutation completed: Sat Sep 14 01:41:27 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 77: Sat Sep 14 01:41:27 2024

    ## Permutation completed: Sat Sep 14 01:41:28 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 77: Sat Sep 14 01:41:28 2024

    ## Permutation completed: Sat Sep 14 01:41:28 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 77: Sat Sep 14 01:41:28 2024

    ## Permutation completed: Sat Sep 14 01:41:28 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 78: Sat Sep 14 01:41:28 2024

    ## Permutation completed: Sat Sep 14 01:41:29 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 78: Sat Sep 14 01:41:29 2024

    ## Permutation completed: Sat Sep 14 01:41:29 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 78: Sat Sep 14 01:41:29 2024

    ## Permutation completed: Sat Sep 14 01:41:29 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 78: Sat Sep 14 01:41:29 2024

    ## Permutation completed: Sat Sep 14 01:41:30 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 78: Sat Sep 14 01:41:30 2024

    ## Permutation completed: Sat Sep 14 01:41:30 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 79: Sat Sep 14 01:41:31 2024

    ## Permutation completed: Sat Sep 14 01:41:31 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 79: Sat Sep 14 01:41:31 2024

    ## Permutation completed: Sat Sep 14 01:41:32 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 79: Sat Sep 14 01:41:32 2024

    ## Permutation completed: Sat Sep 14 01:41:32 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 79: Sat Sep 14 01:41:32 2024

    ## Permutation completed: Sat Sep 14 01:41:32 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 79: Sat Sep 14 01:41:32 2024

    ## Permutation completed: Sat Sep 14 01:41:32 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 80: Sat Sep 14 01:41:32 2024

    ## Permutation completed: Sat Sep 14 01:41:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 80: Sat Sep 14 01:41:33 2024

    ## Permutation completed: Sat Sep 14 01:41:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 80: Sat Sep 14 01:41:33 2024

    ## Permutation completed: Sat Sep 14 01:41:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 80: Sat Sep 14 01:41:33 2024

    ## Permutation completed: Sat Sep 14 01:41:33 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 80: Sat Sep 14 01:41:33 2024

    ## Permutation completed: Sat Sep 14 01:41:33 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 81: Sat Sep 14 01:41:34 2024

    ## Permutation completed: Sat Sep 14 01:41:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 81: Sat Sep 14 01:41:34 2024

    ## Permutation completed: Sat Sep 14 01:41:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 81: Sat Sep 14 01:41:34 2024

    ## Permutation completed: Sat Sep 14 01:41:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 81: Sat Sep 14 01:41:34 2024

    ## Permutation completed: Sat Sep 14 01:41:34 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 81: Sat Sep 14 01:41:34 2024

    ## Permutation completed: Sat Sep 14 01:41:34 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 82: Sat Sep 14 01:41:35 2024

    ## Permutation completed: Sat Sep 14 01:41:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 82: Sat Sep 14 01:41:35 2024

    ## Permutation completed: Sat Sep 14 01:41:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 82: Sat Sep 14 01:41:35 2024

    ## Permutation completed: Sat Sep 14 01:41:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 82: Sat Sep 14 01:41:35 2024

    ## Permutation completed: Sat Sep 14 01:41:35 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 82: Sat Sep 14 01:41:35 2024

    ## Permutation completed: Sat Sep 14 01:41:36 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 83: Sat Sep 14 01:41:36 2024

    ## Permutation completed: Sat Sep 14 01:41:36 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 83: Sat Sep 14 01:41:36 2024

    ## Permutation completed: Sat Sep 14 01:41:36 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 83: Sat Sep 14 01:41:36 2024

    ## Permutation completed: Sat Sep 14 01:41:37 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 83: Sat Sep 14 01:41:37 2024

    ## Permutation completed: Sat Sep 14 01:41:37 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 83: Sat Sep 14 01:41:37 2024

    ## Permutation completed: Sat Sep 14 01:41:37 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 84: Sat Sep 14 01:41:38 2024

    ## Permutation completed: Sat Sep 14 01:41:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 84: Sat Sep 14 01:41:38 2024

    ## Permutation completed: Sat Sep 14 01:41:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 84: Sat Sep 14 01:41:38 2024

    ## Permutation completed: Sat Sep 14 01:41:38 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 84: Sat Sep 14 01:41:38 2024

    ## Permutation completed: Sat Sep 14 01:41:39 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 84: Sat Sep 14 01:41:39 2024

    ## Permutation completed: Sat Sep 14 01:41:39 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 85: Sat Sep 14 01:41:40 2024

    ## Permutation completed: Sat Sep 14 01:41:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 85: Sat Sep 14 01:41:40 2024

    ## Permutation completed: Sat Sep 14 01:41:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 85: Sat Sep 14 01:41:40 2024

    ## Permutation completed: Sat Sep 14 01:41:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 85: Sat Sep 14 01:41:40 2024

    ## Permutation completed: Sat Sep 14 01:41:40 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 85: Sat Sep 14 01:41:40 2024

    ## Permutation completed: Sat Sep 14 01:41:40 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 86: Sat Sep 14 01:41:41 2024

    ## Permutation completed: Sat Sep 14 01:41:41 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 86: Sat Sep 14 01:41:41 2024

    ## Permutation completed: Sat Sep 14 01:41:41 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 86: Sat Sep 14 01:41:41 2024

    ## Permutation completed: Sat Sep 14 01:41:41 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 86: Sat Sep 14 01:41:41 2024

    ## Permutation completed: Sat Sep 14 01:41:42 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 86: Sat Sep 14 01:41:42 2024

    ## Permutation completed: Sat Sep 14 01:41:42 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 87: Sat Sep 14 01:41:42 2024

    ## Permutation completed: Sat Sep 14 01:41:43 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 87: Sat Sep 14 01:41:43 2024

    ## Permutation completed: Sat Sep 14 01:41:43 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 87: Sat Sep 14 01:41:43 2024

    ## Permutation completed: Sat Sep 14 01:41:44 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 87: Sat Sep 14 01:41:44 2024

    ## Permutation completed: Sat Sep 14 01:41:44 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 87: Sat Sep 14 01:41:44 2024

    ## Permutation completed: Sat Sep 14 01:41:45 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 88: Sat Sep 14 01:41:46 2024

    ## Permutation completed: Sat Sep 14 01:41:46 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 88: Sat Sep 14 01:41:46 2024

    ## Permutation completed: Sat Sep 14 01:41:46 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 88: Sat Sep 14 01:41:46 2024

    ## Permutation completed: Sat Sep 14 01:41:46 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 88: Sat Sep 14 01:41:46 2024

    ## Permutation completed: Sat Sep 14 01:41:47 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 88: Sat Sep 14 01:41:47 2024

    ## Permutation completed: Sat Sep 14 01:41:47 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 89: Sat Sep 14 01:41:47 2024

    ## Permutation completed: Sat Sep 14 01:41:47 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 89: Sat Sep 14 01:41:48 2024

    ## Permutation completed: Sat Sep 14 01:41:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 89: Sat Sep 14 01:41:48 2024

    ## Permutation completed: Sat Sep 14 01:41:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 89: Sat Sep 14 01:41:48 2024

    ## Permutation completed: Sat Sep 14 01:41:48 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 89: Sat Sep 14 01:41:48 2024

    ## Permutation completed: Sat Sep 14 01:41:48 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 90: Sat Sep 14 01:41:48 2024

    ## Permutation completed: Sat Sep 14 01:41:49 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 90: Sat Sep 14 01:41:49 2024

    ## Permutation completed: Sat Sep 14 01:41:49 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 90: Sat Sep 14 01:41:49 2024

    ## Permutation completed: Sat Sep 14 01:41:49 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 90: Sat Sep 14 01:41:49 2024

    ## Permutation completed: Sat Sep 14 01:41:49 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 90: Sat Sep 14 01:41:49 2024

    ## Permutation completed: Sat Sep 14 01:41:49 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 91: Sat Sep 14 01:41:49 2024

    ## Permutation completed: Sat Sep 14 01:41:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 91: Sat Sep 14 01:41:50 2024

    ## Permutation completed: Sat Sep 14 01:41:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 91: Sat Sep 14 01:41:50 2024

    ## Permutation completed: Sat Sep 14 01:41:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 91: Sat Sep 14 01:41:50 2024

    ## Permutation completed: Sat Sep 14 01:41:50 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 91: Sat Sep 14 01:41:51 2024

    ## Permutation completed: Sat Sep 14 01:41:51 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 92: Sat Sep 14 01:41:52 2024

    ## Permutation completed: Sat Sep 14 01:41:52 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 92: Sat Sep 14 01:41:52 2024

    ## Permutation completed: Sat Sep 14 01:41:52 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 92: Sat Sep 14 01:41:52 2024

    ## Permutation completed: Sat Sep 14 01:41:52 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 92: Sat Sep 14 01:41:52 2024

    ## Permutation completed: Sat Sep 14 01:41:53 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 92: Sat Sep 14 01:41:53 2024

    ## Permutation completed: Sat Sep 14 01:41:53 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 93: Sat Sep 14 01:41:53 2024

    ## Permutation completed: Sat Sep 14 01:41:54 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 93: Sat Sep 14 01:41:54 2024

    ## Permutation completed: Sat Sep 14 01:41:54 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 93: Sat Sep 14 01:41:54 2024

    ## Permutation completed: Sat Sep 14 01:41:54 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 93: Sat Sep 14 01:41:54 2024

    ## Permutation completed: Sat Sep 14 01:41:54 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 93: Sat Sep 14 01:41:54 2024

    ## Permutation completed: Sat Sep 14 01:41:55 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 94: Sat Sep 14 01:41:55 2024

    ## Permutation completed: Sat Sep 14 01:41:55 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 94: Sat Sep 14 01:41:55 2024

    ## Permutation completed: Sat Sep 14 01:41:55 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 94: Sat Sep 14 01:41:55 2024

    ## Permutation completed: Sat Sep 14 01:41:55 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 94: Sat Sep 14 01:41:56 2024

    ## Permutation completed: Sat Sep 14 01:41:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 94: Sat Sep 14 01:41:56 2024

    ## Permutation completed: Sat Sep 14 01:41:56 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 95: Sat Sep 14 01:41:56 2024

    ## Permutation completed: Sat Sep 14 01:41:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 95: Sat Sep 14 01:41:56 2024

    ## Permutation completed: Sat Sep 14 01:41:56 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 95: Sat Sep 14 01:41:56 2024

    ## Permutation completed: Sat Sep 14 01:41:57 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 95: Sat Sep 14 01:41:57 2024

    ## Permutation completed: Sat Sep 14 01:41:57 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 95: Sat Sep 14 01:41:57 2024

    ## Permutation completed: Sat Sep 14 01:41:57 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 96: Sat Sep 14 01:41:57 2024

    ## Permutation completed: Sat Sep 14 01:41:57 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 96: Sat Sep 14 01:41:57 2024

    ## Permutation completed: Sat Sep 14 01:41:57 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 96: Sat Sep 14 01:41:58 2024

    ## Permutation completed: Sat Sep 14 01:41:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 96: Sat Sep 14 01:41:58 2024

    ## Permutation completed: Sat Sep 14 01:41:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 96: Sat Sep 14 01:41:58 2024

    ## Permutation completed: Sat Sep 14 01:41:58 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 97: Sat Sep 14 01:41:58 2024

    ## Permutation completed: Sat Sep 14 01:41:58 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 97: Sat Sep 14 01:41:58 2024

    ## Permutation completed: Sat Sep 14 01:41:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 97: Sat Sep 14 01:41:59 2024

    ## Permutation completed: Sat Sep 14 01:41:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 97: Sat Sep 14 01:41:59 2024

    ## Permutation completed: Sat Sep 14 01:41:59 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 97: Sat Sep 14 01:41:59 2024

    ## Permutation completed: Sat Sep 14 01:41:59 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 98: Sat Sep 14 01:42:00 2024

    ## Permutation completed: Sat Sep 14 01:42:01 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 98: Sat Sep 14 01:42:01 2024

    ## Permutation completed: Sat Sep 14 01:42:01 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 98: Sat Sep 14 01:42:01 2024

    ## Permutation completed: Sat Sep 14 01:42:01 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 98: Sat Sep 14 01:42:01 2024

    ## Permutation completed: Sat Sep 14 01:42:02 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 98: Sat Sep 14 01:42:02 2024

    ## Permutation completed: Sat Sep 14 01:42:02 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 99: Sat Sep 14 01:42:03 2024

    ## Permutation completed: Sat Sep 14 01:42:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 99: Sat Sep 14 01:42:03 2024

    ## Permutation completed: Sat Sep 14 01:42:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 99: Sat Sep 14 01:42:03 2024

    ## Permutation completed: Sat Sep 14 01:42:03 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 99: Sat Sep 14 01:42:04 2024

    ## Permutation completed: Sat Sep 14 01:42:04 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 99: Sat Sep 14 01:42:04 2024

    ## Permutation completed: Sat Sep 14 01:42:04 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 100: Sat Sep 14 01:42:04 2024

    ## Permutation completed: Sat Sep 14 01:42:04 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 100: Sat Sep 14 01:42:04 2024

    ## Permutation completed: Sat Sep 14 01:42:04 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 100: Sat Sep 14 01:42:04 2024

    ## Permutation completed: Sat Sep 14 01:42:04 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 100: Sat Sep 14 01:42:04 2024

    ## Permutation completed: Sat Sep 14 01:42:05 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting chromosome 100: Sat Sep 14 01:42:05 2024

    ## Permutation completed: Sat Sep 14 01:42:05 2024

    ## ARTP2 supports multi-threading on this OS
    ## ARTP2 supports multi-threading on this OS

    ## Computing pathway p-value: Sat Sep 14 01:42:05 2024

``` r
res=list(gene_pvalue=res_temp$gene.pvalue,path_pvalue=res_temp$pathway.pvalue)

res$path_pvalue
```

    ## [1] 0.003779962
