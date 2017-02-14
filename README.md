# gsasnp2
GSA-SNP2 is successor of GSA-SNP (Nam et al. 2010, NAR web server issue). GSA-SNP2 accepts human GWAS summary data (rs numbers, p-values) or gene-wise p-values (possibly obtained from VEGAS or GATES) and outputs pathway gene sets ‘enriched’ with genes associated with the given phenotype. 

Project website: https://sites.google.com/view/gsasnp2

# Five features of GSA-SNP2:

1/  Reasonable type I error control by the following two processes:
  A)  Gene scores are ‘adjusted’ to the number of SNPs assigned to each gene using monotone cubic spline trend curve.
  B)  Adjacent genes with high inter-gene correlations within each gene set were removed: Inter-gene genotype correlations were provided for five races from 1000 genome data.

2/ High power based on the random set model (Newton et al. 2007)

3/ Without any critical free parameter

4/ Protein interaction networks among the member genes were visualized for the significant pathways. This function enables the user to prioritize core sub-networks within each pathway. STRING and HIPPIE network data are currently provided.  

5/ Easy to use: Only requires GWAS summary data (or gene p-values) and takes only seconds to minutes to get results. Other powerful self-contained pathway tools also require SNP correlation input and take from several days to two weeks. User can also upload their own gene-sets and protein interaction networks.

# Why competitive pathway analysis for GWAS data? 

Competitive methods directly target pathway-level aberrations by testing the ‘enrichment’ of association signals within each pathway set. On the other hand, self-contained methods test the ‘existence’ of association signal in each pathway set. Self-contained methods are in general highly sensitive, so are useful in finding novel pathways. However, genes typically have multiple functions and mere existence of association signal does not necessarily imply the pathway-level aberration. So, both the approaches are useful and complementary to each other. Unfortunately, many competitive methods for GWAS data so far exhibited low powers and were susceptible to some free parameters. 

# Comparison with other competitive methods 

Performance of GSA-SNP2 was compared with those of five existing competitive methods, GSA-SNP, iGSEA4GWAS MAGENTA, INRICH and GOWINDA. GSA-SNP2 was a little liberal in false positive control compared to others, but exhibited high power and best discriminatory ability. 

1/ Comparison of type I error control: twenty null genotype data sets were generated using 1000 genome (European) and GWAsimulator tool and corresponding p-values were input to each program. GSA-SNP2 exhibited greatly improved type I error control compared to GSA-SNP.
  
2/ Power comparison: DIAGRAM consortium GWAS p-values (European) were used to compare the statistical power. 16 curated T2D related pathways (Morris et al. Nat. Genetics 2012) as well as the terms including ‘diabetes’ were regarded as true positives (TPs). GSA-SNP2 exhibited high power and best ranks of TPs. See the results here.
