# The R package ARTP2
It is increasingly recognized that pathway analyses, a joint test of association between the outcome
and a group of single nucleotide polymorphisms (SNPs) within a biological pathway, could potentially
complement single-SNP analysis and provide additional insights for the genetic architecture
of complex diseases. Building upon existing P-value combining methods, we propose a class of
highly flexible pathway analysis approaches based on an adaptive rank truncated product statistic
that can effectively combine evidence of associations over different SNPs and genes within a
pathway. The statistical significance of the pathway-level test statistics is evaluated using a highly
efficient permutation algorithm that remains computationally feasible irrespective of the size of
the pathway and complexity of the underlying test statistics for summarizing SNP- and gene-level
associations.
The main functions in this package for a single population are sARTP when only summary level data
are available, rARTP when genotype data are available, and warm.start for computing gene and
pathway p-values when previously saved information is available.
This package can also merge results from multiple independent populations and then perform the
ARTP test using the merged data. The merging of results can be done at the gene level and
also at the SNP level. At the gene level, the function warm.start.multiPop.gene is called;
at the SNP level the functions sARTP.multiPop.SNP.union, warm.start.multiPop.SNP.union,
sARTP.multiPop.SNP.intersect and warm.start.multiPop.SNP.intersect are used.
