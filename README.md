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
