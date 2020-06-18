---
title: "sparseWeightBasedPCA Package"
author: "Niek C. de Schipper <n.c.deschipper@uvt.nl>"
geometry: margin=1cm
date: "2020-06-17"
output:
    html_document:
        keep_md: true
    pdf_document:
        toc: true
    github_document:
        toc: true
    rmarkdown::html_vignette:
        toc: true
vignette: >
  %\VignetteIndexEntry{bayespca Package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# sparseWeightBasedPCA: A package for Regularized weight based Simultaneous Component Analysis (SCA) and Principal Component Analysis (PCA)

## Theoretical background
Principal component analysis (PCA) is a widely used analysis technique for data reduction. It can give crucial insights in the underlying structure of the data when used as a latent variable model.

Given a data matrix $\mathbf{X}$ that contains the scores for $i = 1...I$ observations on $j = 1...J$ variables; we follow the convention to present the $J$ variable scores of observation $i$ in row $i$ and thus $\mathbf{X}$ has size $I \times J$. PCA decomposes the data into $Q$ components as follows,
\begin{equation}
    \label{PCA}
    \begin{aligned}
        &\mathbf{X} = \mathbf{XWP}^T + \mathbf{E} \\
        &\text{subject to } \mathbf{P}^T\mathbf{P} = \mathbf{I},
      \end{aligned}
\end{equation}
where $\mathbf{W}$ is a $J \times Q$ component weight matrix, $\mathbf{P}$ is a $J \times Q$ loading matrix and $\mathbf{E}$ is a $I \times J$ residual matrix. The component weight matrix $\mathbf{W}$ will be the focus of this package, note that $\mathbf{XW}$ represent the component scores. 

The decomposition in \ref{PCA}  can be extended to the case of multi-block data by taking $\mathbf{X} = [\mathbf{X}_{1} \hdots \mathbf{X}_{K}]$; this is concatenating the $K$ data blocks composed of different sets of variables for the same observation units. The decomposition of $\mathbf{X}$ has the same block structured decomposition with $\mathbf{W} = [\mathbf{W}^{T}_{1} \hdots \mathbf{W}^{T}_{K}]^{T}$ and  $\mathbf{P} = [\mathbf{P}^{T}_{1} \hdots \mathbf{P}^{T}_{K}]^{T}$. This multi-block formulation of PCA is known as simultaneous component analysis (SCA).

Also in the multi-bock case $\mathbf{W}$ can be penalized to obtain sparse weights, we will call this variant sparse SCA When analyzing multi-block data with sparse SCA, we can search for blockwise structures in the component weights that tell us whether a component is uniquely determined by variables from one single data block (distinctive component), or whether it is a component that is determined by variables from multiple data blocks (common component). In other words, a distinctive component is a linear combination of variables of a particular data block only, whereas a common component is a linear combination of variables of multiple data blocks. An example of common and distinctive components in the situation with two data blocks is given below. The first two components are distinctive components, the third component is a common component, 
\[ \mathbf{T}
    =
    \begin{bmatrix}
        \mathbf{X}_1 & \mathbf{X}_2
    \end{bmatrix}
    \begin{bmatrix}
        \mathbf{W}_1 \\
        \mathbf{W}_2
    \end{bmatrix}
    =
    \begin{bmatrix}
        \mathbf{X}_1 & \mathbf{X}_2
    \end{bmatrix}
\begin{bmatrix}
    0 & w_{1_{1}2} & w_{1_{1}3}   \\
    0 & w_{2_{1}2} & w_{2_{1}3}   \\
    0 & w_{3_{1}2} & w_{3_{1}3}  
    \\[5pt]
    w_{1_{2}1} & 0 & w_{1_{2}3} \\
    w_{2_{2}1} & 0 & w_{2_{2}3} \\
    w_{3_{2}1} & 0 & w_{2_{2}3} 
\end{bmatrix}. \]

The `sparseWeightBasedPCA` package performs the following tasks:

1. Task 1
2. Task 2
3. Task 3

## The `sparseWeightBasedPCA` package


```r
set.seed(1)
X <- matrix(rnorm(100), 20, 5)
print(X)
```

```
##              [,1]        [,2]       [,3]         [,4]       [,5]
##  [1,] -0.62645381  0.91897737 -0.1645236  2.401617761 -0.5686687
##  [2,]  0.18364332  0.78213630 -0.2533617 -0.039240003 -0.1351786
##  [3,] -0.83562861  0.07456498  0.6969634  0.689739362  1.1780870
##  [4,]  1.59528080 -1.98935170  0.5566632  0.028002159 -1.5235668
##  [5,]  0.32950777  0.61982575 -0.6887557 -0.743273209  0.5939462
##  [6,] -0.82046838 -0.05612874 -0.7074952  0.188792300  0.3329504
##  [7,]  0.48742905 -0.15579551  0.3645820 -1.804958629  1.0630998
##  [8,]  0.73832471 -1.47075238  0.7685329  1.465554862 -0.3041839
##  [9,]  0.57578135 -0.47815006 -0.1123462  0.153253338  0.3700188
## [10,] -0.30538839  0.41794156  0.8811077  2.172611670  0.2670988
## [11,]  1.51178117  1.35867955  0.3981059  0.475509529 -0.5425200
## [12,]  0.38984324 -0.10278773 -0.6120264 -0.709946431  1.2078678
## [13,] -0.62124058  0.38767161  0.3411197  0.610726353  1.1604026
## [14,] -2.21469989 -0.05380504 -1.1293631 -0.934097632  0.7002136
## [15,]  1.12493092 -1.37705956  1.4330237 -1.253633400  1.5868335
## [16,] -0.04493361 -0.41499456  1.9803999  0.291446236  0.5584864
## [17,] -0.01619026 -0.39428995 -0.3672215 -0.443291873 -1.2765922
## [18,]  0.94383621 -0.05931340 -1.0441346  0.001105352 -0.5732654
## [19,]  0.82122120  1.10002537  0.5697196  0.074341324 -1.2246126
## [20,]  0.59390132  0.76317575 -0.1350546 -0.589520946 -0.4734006
```











