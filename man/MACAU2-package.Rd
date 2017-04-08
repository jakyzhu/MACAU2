\name{MACAU2-package}
\alias{MACAU2-package}
\alias{MACAU2}
\docType{package}
\title{
Efficient Mixed Model Analysis of Count Data in Large-Scale Genomic Sequencing Studies
}
\description{
MACAU2 is an R package for efficient differerntial analysis of large-scale RNA sequencing (RNAseq) data and Bisulfite sequencing (BSseq) data in the presence of individual relatedness and population structure. MACAU2 first fits a GLMM with adjusted covariates, predictor of interest and random effects to account for population structure and individual relatedness, and then performs wald tests for each gene in RNAseq or site for BSseq. 
}
\details{
\tabular{ll}{
Package: \tab MACAU2\cr
Type: \tab Package\cr
Version: \tab 1.10\cr
Date: \tab 2017-03-31\cr
License: \tab GPL-3\cr
}
}
\author{
Shiquan Sun, Jiaqiang Zhu, Xiang Zhou

Maintainer: Shiquan Sun <shiquans@umich.edu>
}
\references{

Breslow, N.E. and Clayton, D.G. (1993) Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88, 9-25.

Chen, H., Wang, C., Conomos, M.P., Stilp, A.M., Li, Z., Sofer, T., Szpiro, A.A., Chen, W., Brehm, J.M., Celedón, J.C., Redline, S., Papanicolaou, G.J., Thornton, T.A., Laurie, C.C., Rice, K. and Lin, X (2016). Control for population structure and relatedness for binary traits in genetic association studies using logistic mixed models. The American Journal of Human Genetics, 98, 653-666.

Gilmour, A.R., Thompson, R. and Cullis, B.R. (1995) Average Information REML: An Efficient Algorithm for Variance Parameter Estimation in Linear Mixed Models. Biometrics 51, 1440-1450.

Lea, A.J., Tung, J. and Zhou,X. (2015) A flexible, effcient binomial mixed model for identifying differential DNA methylation in bisulfite sequencing data. PLoS Genetics. 11: e1005650.

Rue H., Martino S. , and Chopin N.(2009) Approximate Bayesian inference for latent Gaussian models using integrated nested Laplace approximations (with discussion). Journal of the Royal Statistical Society, Series B,71(2):319-392. 

Yang, J., Lee, S.H., Goddard, M.E. and Visscher, P.M. (2011) GCTA: A Tool for Genome-wide Complex Trait Analysis. The American Journal of Human Genetics 88, 76-82.

Zhou, X., Carbonetto, P. and Stephens,M. (2013) Polygenic modeling with Bayesian sparse linear mixed models. PLoS Genetics. 9(2): e1003264.

Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nature Genetics 44, 821-824.

Zhou, X. and Stephens, M.(2014) Efficient multivariate linear mixed model algorithms for genome-wide association studies. Nature Methods.

}
\keyword{ package }
\keyword{ GLMMs, MACAU2}


