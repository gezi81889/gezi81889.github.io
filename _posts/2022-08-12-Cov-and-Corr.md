---
layout: post
title: Covariance and correlation
date: 2022-08-12 11:12:00-0400
comments: true
description: 
tags: math covariance correlation
categories: math-notes
---
  

### Variance vs. Covariance  

Variance refers to the spread of a variable around its mean value, while a covariance refers to the measure of the directional relationship between two random variables. 

Suppose we have a dataset with variable $$X$$ and $$Y$$. For a variable $$X$$, it has n values $$[x_1,x_2...x_n]$$ , 

$$
Var(X)= \frac{1}{n-1}\sum(x_i-\overline x)^2
$$

For another variable $$Y$$ in this dataset, we can also calculate its variance. But $$Var(X)$$ and $$Var(Y)$$ don't give us information about the distribution of the data in 2D space, meaning the trend of data in higher space. The covariance quantifies the **direction of the trend**, but it does not say how strong the trend is and covariance is very **sensitive to scaling**. 

$$
Cov(X,Y) = \frac{1}{n-1}\sum(x_i-\overline x)(y_i-\overline y)=E[(X-E(X))(Y-E(Y))]=E(XY)-E(X)E(Y)
$$

For a positive covariance, $$X$$ and $$Y$$ are incrasing together; For a negative value, $$X$$ is increasing but $$ Y $$ is decreasing. 

  

### Correlation 

Correlation between two variables measures the **strength of the trend**. After normalization, correlation coefficient does not get affected by scaling. Both covariance and correlation are used for quantifying the "Independence" between ramdon variables. 

$$
Cor(X,Y) = \frac{Cor(X,Y)}{\sqrt{Var(X)Var(Y)}}
$$
  

### Covariance matrix 

Covariance matrix contain variances and covariance between variables in our data, which means we have variables in both column and row. Suppose we have $$p$$ cells $$\times$$ $$n$$ genes in dataset $$X$$ and $$q$$ cells $$\times $$ $$m$$ genes in dataset $$Y$$ (Note that $$X$$ and $$Y$$ are now datasets not variables ): 

$$
X=\begin {pmatrix} x_{11} & x_{12} & ... & x_{1n} \\
                   x_{21} & x_{22} & ... & x_{2n} \\
                   . & . & ... & . \\
                   . & . & ... & . \\
                   x_{p1} & x_{p2} & ... & x_{pn}  \end{pmatrix} = [\vec x_1,\vec x_2,...\vec x_p]’
$$

$$
Y=\begin {pmatrix} y_{11} & y_{12} & ... & y_{1n} \\
                   y_{21} & y_{22} & ... & y_{2n} \\
                   . & . & ... & . \\
                   . & . & ... & . \\
                   y_{q1} & y_{q2} & ... & y_{qn}  \end{pmatrix} = [\vec y_1,\vec y_2,...\vec y_q]’ 
$$

Mean-centered (each column has zero expectation):

$$
\begin{align*}\label{1}
\overline x &=\frac{1}{p}\sum_{i=1}^p \vec x_{i} = [E(gene_1),E(gene_2),...E(gene_n)]
\\
\widetilde X &=[\vec x_1-\overline x,\vec x_2-\overline x,...\vec x_p-\overline x]'=[gene_1-E(gene_1),gene_2-E(gene_2),....gene_n-E(gene_n)]
\\
\overline y &=\frac{1}{q}\sum_{i=1}^q \vec y_{i}= [E(gene_1),E(gene_2),...E(gene_m)]
\\
\widetilde Y &=[\vec y_1-\overline y,\vec y_2-\overline y,...\vec y_q-\overline y]'=[gene_1-E(gene_1),gene_2-E(gene_2),....gene_m-E(gene_m)]
\end{align*}
$$

Covariance matrix is the matrix whose $(i,j)$ entry is the covariance between the variables $gene_i, gene_j$:

$$
\begin{align*}\label{2}
Var(\vec X)_{i,j} &= Cov(gene_i, gene_j)=>Var(\vec X)= Cov(\vec X,\vec X)=\frac{1}{p-1}\widetilde X^T\widetilde X
\\
Var(\vec Y)_{i,j} &=Cov(gene_i, gene_j)=>Var(\vec Y)=Cov(\vec Y,\vec Y) = \frac{1}{q-1}\widetilde Y^T\widetilde Y
\end{align*}
$$

If we denote the variable vector as $$\vec X=[gene_1,gene_2,...gene_n]'$$ (sometime also called random vectors), the covariance matrix of dataset $$X$$ can also be considered as covariance matrix of the random vector $$\vec X$$ , which is typically denoted by $$K_{XX}$$ or $$\Sigma_{XX}$$ . 

  
***Characteristics of Covariance matrix $$A=Var(\vec X)$$ ***

1. $$n \times n$$ Square symmetric -> $$A = A^T$$

2. $$Var(\vec X)_{ii}=Cov(gene_i,gene_i)=Var(gene_i)$$. 

3. Positive semi-definite -> For every non-zero vector $$\vec z$$,  $$z'Az≥0$$  

4. If $$gene_1,gene_2,....gene_n$$ are all independent,$$Var(\vec X)=diag(Var(gene_1),Var(gene_2),...Var(gene_n))$$

    

Similarly we can also have $$K_{XY}$$ or $$\Sigma_{XY}$$ , which is the **cross-covariance matrix** of random vector $$\vec X$$ and $$\vec Y$$ of dataset $$X$$ and $$Y$$. $$(i,j)$$ entry is the covariance between the variables $$gene_i$$ from dataset $$X$$ and  $$gene_j$$ from dataset $$Y$$ (note that for calculation, we have to subset two datasets to make them have same observations/cells $$r$$): 

$$
Cov(\vec X,\vec Y)_{i,j} =Cov(gene_i, gene_j)=>Cov(\vec X,\vec Y)= \frac{1}{r-1}\widetilde X^T\widetilde Y
$$

Cross-covariance can measure how dependent two variable set of two datasets are, but it also depends on the unit of $$\vec X$$ and $$\vec Y$$, that is to say, they are also sensitive to scaling. 

  

***Characteristics of Cross-coariance matrix $$A=Cov(\vec X, \vec Y), B=Cov(\vec Y,\vec X)^T$$***

1. $$n \times m$$ matrix, in general not symmetric
2. $$A= B^T$$.
3. $$Cov(\vec X_1+\vec X_2,\vec Y)=Cov(\vec X_1,\vec Y)+Cov(\vec X_2,\vec Y)$$. 
4. If $$a$$ is a constant, $$Cov(a\vec X,\vec Y)=Cov(\vec X,a\vec Y)=aA$$ ; $$Cov(a+\vec X,\vec Y)=Cov(\vec X,a +\vec Y)=A$$ (Linear recombination is discussed below)
5. If $$n=m$$, $$Var(\vec X+\vec Y)=Var(\vec X)+Var(\vec Y)+A+B$$

  

### Covariance of linear recombinations

If constants $$a_i,b_i \in R$$ , $$i=1,2,...,n$$ , then the following is true:

$$
\begin{align*}\label{3}
& Cov(a_1gene_1+a_2gene_2+....+a_ngene_n,b_1gene1+b_2gene2+...+b_ngene_n) 
\\
& =Cov(\vec a \cdot \vec X, \vec b \cdot \vec X) = 
\sum_{i=1}^n\sum_{j=1}^n a_ib_jCov(gene_i,gene_j)=\vec a^T\Sigma_{XX}\vec b
\end{align*}
$$

$$\vec a=[a_1,a_2,...a_n]', \vec b=[b_1,b_2,...b_n]'$$. Similarly,

$$
\begin{align*}\label{4}
& Cov(a_1gene_1+a_2gene_2+....+a_ngene_n,a_1gene1+a_2gene2+...+a_ngene_n) 
\\
& =Cov(\vec a \cdot \vec X, \vec a \cdot \vec X) = Var(\vec a \cdot \vec X)=
\sum_{i=1}^n\sum_{j=1}^n a_ia_jCov(gene_i,gene_j)=\vec a^TVar(\vec X)\vec a 
\end{align*}
$$

Recalled $$Var(aX)=a^2Var(X)$$, $$X$$ is a random variable. 

  

For a square $$n \times n$$ matrix $$A=[\vec a_1',\vec a_2',...\vec a_n']'$$ and a square $$n \times n$$ matrix $$B=[\vec  b_1',\vec b_2',...\vec b_n']'$$, $$A\vec X$$ is the vector of linear recombination of $$X_i$$ variables, we have: 

$$
Var(A\vec X)= AVar(\vec X)A^T
\\
Cor(A\vec X,B\vec X) = AVar(\vec X)B^T
$$


In summary, for ramdon vectors $$\vec X$$ and $$\vec Y$$, $$\vec X=[X_1,X_2,...X_n], \vec Y=[Y_1,Y_2,...Y_m]$$ , $$X_i$$ and $$Y_i$$ are random variables (In our cases above, gene variable). 

1. $$Cov(X_i,X_i)=Var(X_i)$$ is a value

2. $$Cov(\vec X,\vec X)=Var(\vec X)=K_{XX}=\Sigma_{XX}$$ is covariance matrix  

3. $$Cov(X_i,Y_i)$$ is a value

4. $$Cov(\vec X,\vec Y)=K_{XY}=\Sigma_{XY}$$ is covariance matrix  

5. $$Cov(\vec a^T \vec X,\vec b^T \vec X)= Cov(X^*,X^*)$$ is a value. $$X^*$$ is the linear recombination of variable $$X_i$$ ; Similarly, $$Cov(\vec a^T \vec X,\vec b^T \vec Y)$$ is a value

   
  

### Correlation matrix

Correlation matrix contain correlations between variables in our data, which means we have variables in both column and row.$$(i,j)$$ entry is the correlation between the variables $$gene_i, gene_j$$:

Standardized matrix: 

$$
\begin{align*}\label{5}
\sigma(gene_i) &=\sqrt{Var(gene_i)}
\\
X_s &=[\frac{gene_1-E(gene_1)}{\sigma(gene_1)},\frac{gene_2-E(gene_2)}{\sigma(gene_2)},....\frac{gene_n-E(gene_n)}{\sigma(gene_n)}]
\\
Y_s &=[\frac{gene_1-E(gene_1)}{\sigma(gene_1)},\frac{gene_2-E(gene_2)}{\sigma(gene_2)},....\frac{gene_n-E(gene_m)}{\sigma(gene_m)}]
\end{align*}
$$

Correlation matrix: 

$$
\begin{align*}\label{6}
&Corr(\vec X)_{ij}=Cor(gene_i,gene_j)->Corr(\vec X)=\frac{1}{p-1}X_s^TX_s
\\
&Corr(\vec Y)_{ij}=Cor(gene_i,gene_j)->Corr(\vec Y)=\frac{1}{q-1}Y_s^TY_s
\end{align*}
$$


Similarly, we also have **Cross-correlation**. Similarly we need to subset the dataset to make them have same observations/cells $$r$$. 

$$
\begin{align*}\label{7}
Corr(\vec X,\vec Y)_{ij} &=Cor(gene_i,gene_j)-> Corr(\vec X,\vec Y)=\frac{1}{r-1}X^T_sY_s
\\
Corr(\vec Y,\vec X)_{ij} &=Cor(gene_i,gene_j)-> Corr(\vec Y,\vec X)=\frac{1}{r-1}Y^T_sX_s
\end{align*}
$$

  

### Reference 

This blog heavily referenced the Course of Youtuber molypath. 

https://www.youtube.com/watch?v=QptKkD__k-c

