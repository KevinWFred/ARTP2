SNP-IVW
================
2024-09-14

## SNP centric SNP-IVW method (see Fu et al. for more details on the method)

An example to run the SNP-IVW procedure on data from five populations:

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

library("ARTP3")

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

``` r
# Record the ending time
end_time <- Sys.time()
# Calculate the time taken
t.union <- end_time - start_time

res=list(gene_pvalue=res_union$gene.pvalue,path_pvalue=res_union$pathway.pvalue,
          running_time=as.numeric(t.union))

res$path_pvalue
```

    ## [1] 0.1756824
