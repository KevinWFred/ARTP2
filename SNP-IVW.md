Union
================
2024-09-14

## Union

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

# ncases <- list()
# nctrls <- list()
# ncases[[1]] <- ncase
# nctrls[[1]] <- nctrl
family <- 'binomial'

## preparing sARTP files for warm.start.multiPop.gene

jobid ="union_example"
outdir=paste0("./", jobid)
if (!dir.exists(outdir)) dir.create(outdir)

# running sARTP.multiPop.SNP.intersect_union
study.list = as.list(studys)
reference.list <- split(reference, 1:nrow(reference))
lambda.list = as.list(rep(lambda,npop))

ncases <- split(ncase, seq_along(ncase))
nctrls <- split(nctrl, seq_along(nctrl))
ncases.list <- split(ncases, seq_along(ncases))
nctrls.list <- split(nctrls, seq_along(ncases))


options_union <- list()

s0 = 10000
for(j in 1:npop){
  options_union[[j]] <- list(maf = .02, HWE.p = 0, min.marg.p=1e-20,
                             gene.R2=0.9, chr.R2=0.9, huge.gene.R2=0.8, huge.chr.R2=0.8,
                             id.str = paste0(text[j],"-union"),
                             out.dir=outdir,
                             nthread = nthread, save.setup = FALSE)
}
options.merged=list(nperm = nperm, id.str = paste0("merge","-union"),
                    out.dir=outdir,
                    nthread = nthread, expand.fold = 0, 
                    seed = s0, inspect.snp.n=2, group.gap=9*10^6)

# Record the starting time
start_time <- Sys.time()
#### running sARTP.multiPop.SNP.union_intersect using union method
res_union <- sARTP.multiPop.SNP.union_intersect(method="union", summary.files.list=study.list,
                                                pathway=pathway_file, family=family, 
                                                reference.list=reference.list, 
                                                lambda.list=lambda.list, 
                                                ncases.list = ncases.list,
                                                ncontrols.list = nctrls.list, 
                                                options.list = options_union,
                                                options.merged=options.merged)
```

    ## Loading definition of pathway: Sat Sep 14 01:25:59 2024

    ## Loading summary statistics: Sat Sep 14 01:25:59 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_AFR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:25:59 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:25:59 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:25:59 2024

    ## Realigning allele information of reference: Sat Sep 14 01:26:00 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:26:00 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:26:00 2024

    ## Removing constant SNPs: Sat Sep 14 01:26:00 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:26:00 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:26:09 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 01:26:33 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 01:26:39 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:27:00 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:27:00 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:27:00 2024

    ## Loading definition of pathway: Sat Sep 14 01:27:00 2024

    ## Loading summary statistics: Sat Sep 14 01:27:00 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_AMR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:27:01 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:27:01 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:27:01 2024

    ## Realigning allele information of reference: Sat Sep 14 01:27:01 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:27:01 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:27:01 2024

    ## Removing constant SNPs: Sat Sep 14 01:27:01 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:27:01 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:27:10 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 01:27:26 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 01:27:33 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:27:44 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:27:44 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:27:44 2024

    ## Loading definition of pathway: Sat Sep 14 01:27:45 2024

    ## Loading summary statistics: Sat Sep 14 01:27:45 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_EAS1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:27:45 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:27:45 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:27:45 2024

    ## Realigning allele information of reference: Sat Sep 14 01:27:45 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:27:45 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:27:46 2024

    ## Removing constant SNPs: Sat Sep 14 01:27:46 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:27:46 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:27:54 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:28:01 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:28:01 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:28:01 2024

    ## Loading definition of pathway: Sat Sep 14 01:28:02 2024

    ## Loading summary statistics: Sat Sep 14 01:28:02 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_EUR1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:28:02 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:28:02 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:28:02 2024

    ## Realigning allele information of reference: Sat Sep 14 01:28:02 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:28:02 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:28:02 2024

    ## Removing constant SNPs: Sat Sep 14 01:28:03 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:28:03 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:28:11 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 01:28:24 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 01:28:29 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:28:39 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:28:39 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:28:39 2024

    ## Loading definition of pathway: Sat Sep 14 01:28:40 2024

    ## Loading summary statistics: Sat Sep 14 01:28:40 2024

    ## Warning in load.summary.statistics(summary.files, pathway, options): Direction
    ## is not found in
    ## /gpfs/gsfs12/users/BB_Bioinformatics/Kevin/tools/Rpackages/ARTP2/extdata/study_SAS1.txt.gz.
    ## ARTP2 assumes equal sample size of SNPs in the study.

    ## Loading allele information from PLINK files: Sat Sep 14 01:28:40 2024

    ## Removing SNPs with conflictive allele information: Sat Sep 14 01:28:40 2024

    ## Loading genotypes from PLINK files: Sat Sep 14 01:28:40 2024

    ## Realigning allele information of reference: Sat Sep 14 01:28:40 2024

    ## Removing SNPs with high missing rate: Sat Sep 14 01:28:41 2024

    ## Removing SNPs with low MAFs: Sat Sep 14 01:28:41 2024

    ## Removing constant SNPs: Sat Sep 14 01:28:41 2024

    ## Removing high LD SNPs within genes: Sat Sep 14 01:28:41 2024

    ## Removing SNPs in high LD within chromosomes: Sat Sep 14 01:28:49 2024

    ## Removing high LD SNPs within genes of huge chromosome: Sat Sep 14 01:29:01 2024

    ## Removing SNPs in high LD within huge chromosomes: Sat Sep 14 01:29:08 2024

    ## Removing genes which are subsets of other genes: Sat Sep 14 01:29:17 2024

    ## Calculating P and SE if not provided: Sat Sep 14 01:29:17 2024

    ## Removing SNPs close to marginal signals: Sat Sep 14 01:29:17 2024

    ## Recovering test statistics: Sat Sep 14 01:29:18 2024

    ## Recovering test statistics: Sat Sep 14 01:29:19 2024
    ## Recovering test statistics: Sat Sep 14 01:29:19 2024

    ## Recovering test statistics: Sat Sep 14 01:29:20 2024
    ## Recovering test statistics: Sat Sep 14 01:29:20 2024

    ## ARTP2 supports multi-threading on this OS

    ## Permuting group 1: Sat Sep 14 01:29:22 2024

    ## Permuting group 2: Sat Sep 14 01:29:23 2024

    ## Permuting group 3: Sat Sep 14 01:29:23 2024

    ## Permuting group 4: Sat Sep 14 01:29:23 2024

    ## Permuting group 5: Sat Sep 14 01:29:25 2024

    ## Permuting group 6: Sat Sep 14 01:29:25 2024

    ## Permuting group 7: Sat Sep 14 01:29:26 2024

    ## Permuting group 8: Sat Sep 14 01:29:27 2024

    ## Permuting group 9: Sat Sep 14 01:29:27 2024

    ## Permuting group 10: Sat Sep 14 01:29:28 2024

    ## Permuting group 11: Sat Sep 14 01:29:28 2024

    ## Permuting group 12: Sat Sep 14 01:29:29 2024

    ## Permuting group 13: Sat Sep 14 01:29:31 2024

    ## Permuting group 14: Sat Sep 14 01:29:31 2024

    ## Permuting group 15: Sat Sep 14 01:29:31 2024

    ## Permuting group 16: Sat Sep 14 01:29:32 2024

    ## Permuting group 17: Sat Sep 14 01:29:32 2024

    ## Permuting group 18: Sat Sep 14 01:29:32 2024

    ## Permuting group 19: Sat Sep 14 01:29:33 2024

    ## Permuting group 20: Sat Sep 14 01:29:34 2024

    ## Permuting group 21: Sat Sep 14 01:29:34 2024

    ## Permuting group 22: Sat Sep 14 01:29:35 2024

    ## Permuting group 23: Sat Sep 14 01:29:35 2024

    ## Permuting group 24: Sat Sep 14 01:29:36 2024

    ## Permuting group 25: Sat Sep 14 01:29:36 2024

    ## Permuting group 26: Sat Sep 14 01:29:37 2024

    ## Permuting group 27: Sat Sep 14 01:29:37 2024

    ## Permuting group 28: Sat Sep 14 01:29:38 2024

    ## Permuting group 29: Sat Sep 14 01:29:38 2024

    ## Permuting group 30: Sat Sep 14 01:29:38 2024

    ## Permuting group 31: Sat Sep 14 01:29:39 2024

    ## Permuting group 32: Sat Sep 14 01:29:39 2024

    ## Permuting group 33: Sat Sep 14 01:29:39 2024

    ## Permuting group 34: Sat Sep 14 01:29:40 2024

    ## Permuting group 35: Sat Sep 14 01:29:41 2024

    ## Permuting group 36: Sat Sep 14 01:29:41 2024

    ## Permuting group 37: Sat Sep 14 01:29:43 2024

    ## Permuting group 38: Sat Sep 14 01:29:44 2024

    ## Permuting group 39: Sat Sep 14 01:29:44 2024

    ## Permuting group 40: Sat Sep 14 01:29:44 2024

    ## Permuting group 41: Sat Sep 14 01:29:45 2024

    ## Permuting group 42: Sat Sep 14 01:29:46 2024

    ## Permuting group 43: Sat Sep 14 01:29:47 2024

    ## Permuting group 44: Sat Sep 14 01:29:48 2024

    ## Permuting group 45: Sat Sep 14 01:29:48 2024

    ## Permuting group 46: Sat Sep 14 01:29:49 2024

    ## Permuting group 47: Sat Sep 14 01:29:49 2024

    ## Permuting group 48: Sat Sep 14 01:29:50 2024

    ## Permuting group 49: Sat Sep 14 01:29:50 2024

    ## Permuting group 50: Sat Sep 14 01:29:51 2024

    ## Permuting group 51: Sat Sep 14 01:29:51 2024

    ## Permuting group 52: Sat Sep 14 01:29:52 2024

    ## Permuting group 53: Sat Sep 14 01:29:52 2024

    ## Permuting group 54: Sat Sep 14 01:29:53 2024

    ## Permuting group 55: Sat Sep 14 01:29:54 2024

    ## Permuting group 56: Sat Sep 14 01:29:54 2024

    ## Permuting group 57: Sat Sep 14 01:29:56 2024

    ## Permuting group 58: Sat Sep 14 01:29:57 2024

    ## Permuting group 59: Sat Sep 14 01:29:57 2024

    ## Permuting group 60: Sat Sep 14 01:29:58 2024

    ## Permuting group 61: Sat Sep 14 01:29:59 2024

    ## Permuting group 62: Sat Sep 14 01:30:00 2024

    ## Permuting group 63: Sat Sep 14 01:30:00 2024

    ## Permuting group 64: Sat Sep 14 01:30:00 2024

    ## Permuting group 65: Sat Sep 14 01:30:02 2024

    ## Permuting group 66: Sat Sep 14 01:30:02 2024

    ## Permuting group 67: Sat Sep 14 01:30:04 2024

    ## Permuting group 68: Sat Sep 14 01:30:05 2024

    ## Permuting group 69: Sat Sep 14 01:30:05 2024

    ## Permuting group 70: Sat Sep 14 01:30:06 2024

    ## Permuting group 71: Sat Sep 14 01:30:07 2024

    ## Permuting group 72: Sat Sep 14 01:30:07 2024

    ## Permuting group 73: Sat Sep 14 01:30:08 2024

    ## Permuting group 74: Sat Sep 14 01:30:08 2024

    ## Permuting group 75: Sat Sep 14 01:30:09 2024

    ## Permuting group 76: Sat Sep 14 01:30:10 2024

    ## Permuting group 77: Sat Sep 14 01:30:10 2024

    ## Permuting group 78: Sat Sep 14 01:30:11 2024

    ## Permuting group 79: Sat Sep 14 01:30:12 2024

    ## Permuting group 80: Sat Sep 14 01:30:12 2024

    ## Permuting group 81: Sat Sep 14 01:30:12 2024

    ## Permuting group 82: Sat Sep 14 01:30:13 2024

    ## Permuting group 83: Sat Sep 14 01:30:13 2024

    ## Permuting group 84: Sat Sep 14 01:30:14 2024

    ## Permuting group 85: Sat Sep 14 01:30:14 2024

    ## Permuting group 86: Sat Sep 14 01:30:15 2024

    ## Permuting group 87: Sat Sep 14 01:30:15 2024

    ## Permuting group 88: Sat Sep 14 01:30:16 2024

    ## Permuting group 89: Sat Sep 14 01:30:17 2024

    ## Permuting group 90: Sat Sep 14 01:30:17 2024

    ## Permuting group 91: Sat Sep 14 01:30:17 2024

    ## Permuting group 92: Sat Sep 14 01:30:18 2024

    ## Permuting group 93: Sat Sep 14 01:30:19 2024

    ## Permuting group 94: Sat Sep 14 01:30:19 2024

    ## Permuting group 95: Sat Sep 14 01:30:19 2024

    ## Permuting group 96: Sat Sep 14 01:30:20 2024

    ## Permuting group 97: Sat Sep 14 01:30:20 2024

    ## Permuting group 98: Sat Sep 14 01:30:21 2024

    ## Permuting group 99: Sat Sep 14 01:30:22 2024

    ## Permuting group 100: Sat Sep 14 01:30:22 2024

    ## Permutation completed: Sat Sep 14 01:30:23 2024

    ## Computing pathway p-value: Sat Sep 14 01:30:23 2024

``` r
# Record the ending time
end_time <- Sys.time()
# Calculate the time taken
t.union <- end_time - start_time

res=list(gene_pvalue=res_union$gene.pvalue,path_pvalue=res_union$pathway.pvalue,
          running_time=as.numeric(t.union))

res$path_pvalue
```

    ## [1] 0.101349
