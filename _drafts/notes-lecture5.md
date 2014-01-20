---
layout: post
title: "Lecture 5: More on inference"
category: 'Lecture'
---

Instructor: Alexandre Bouchard-C&ocirc;t&eacute;

Editor: TBA

***Schedule office hours to make sure nobody is lost***

### Theoretical properties of DPs (continued)

#### Summary of material from last time



#### Moments of a DP

In this section, we derive the  first and second moments of $G(A)$, for $G \sim \dirp(\alpha\_0, G\_0)$ and  $A \subset \Omega$.  To do that, we use the Kolmogorov definition and consider the partition $(A, A^c)$ of $\Omega$.  We get:
\\begin{eqnarray}
(G(A), G(A^c)) \sim \mbox{Dir}(\alpha\_0G\_0(A), \alpha\_0G\_0(A^c)).
\\end{eqnarray}
This implies that:
\\begin{eqnarray}
G(A) \sim \mbox{Beta}(F, G),
\\end{eqnarray}
where $x$ denotes $\alpha\_0G\_0(A)$, and $y$ denotes $\alpha\_0G\_0(A^c)$.

The first moment of  $G(A)$ is therefore
\\begin{eqnarray}
\E[G(A)]=\frac{x}{x+y}=\frac{\alpha\_0G\_0(A)}{\alpha\_0G\_0(A)+\alpha\_0G\_0(A^c)}=G\_0(A), 
\\end{eqnarray}
and the second moment of  $G(A)$ is 
\\begin{eqnarray}
\var[G(A)] &=& \frac{xy}{(x+y)^2(1+x+y)} \\\\
&=&\frac{\alpha\_0^2G\_0(A)(1-G\_0(A))}{\alpha\_0^2(\alpha\_0+1)} \\\\
&=&\frac{G\_0(A)(1-G\_0(A))}{\alpha\_0+1}.
\\end{eqnarray}

This gives an interpretation of $\alpha\_0$ as a precision parameter for the Dirichlet process. 

### Inference using DPMs: an end-to-end example

**Goal:** show the detail of using a DP to answer a statistical question. We will introduce along the way two MCMC techniques that can be used to approximate the posterior.  We will also look at the question of how to use the posterior to answer the task at hand.

**Tasks:** 

- The task in density estimation is to give an estimate, based on observed data, of an unobservable underlying probability density function. The unobservable density function is thought of as the density according to which a large population is distributed; the data are usually thought of as a random sample from that population.  
- In cases where the population is thought as being the union of sub-populations, the task of cluster analysis is to find the sub-population structure (usually without labeled data).  Let us assume for simplicity that we wish to separate the data into two clusters. <sub>Note that the Dirichlet process is still a useful tool, even when the number of desired cluster is fixed.  This is because each cluster that is output may need internally more than one mixture to be explained adequately under the likelihood model at hand.  </sub>

***Cautionary note regarding using DP for clustering***



**Bayesian approach:** Let us take a Bayesian approach to these problems.  This means that we (the modeler) need to pick:

- A joint probability distribution over knowns and unknowns (prior+likelihood): we will  pick the DPM as defined last week.
- A loss function: let us look at a detailed example. 

#### Examples of loss functions

In the case of clustering, a popular choice is the rand loss between a true and putative labeled partitions  $\rhot, \rhop$, denoted by $\randindex(\rhot, \rhop)$.<sub>Note that we turn the standard notion of rand index into a loss by taking 1 - the rand index.</sub>

---

**Definition:**
The rand loss is defined as the number of (unordered) pairs of data points indices $\{i,j\}$ such that $(i \sim\_{\rhot} j) \neq (i \sim\_{\rhop} j)$, i.e.:
\begin{eqnarray}
\sum\_{1\le i < j \le n} \1[(i \sim\_\rhot j) \neq (i \sim\_\rhop j)],
\end{eqnarray}
where:
\begin{eqnarray} 
(i \sim\_\rho j) = \left\{ \begin{array}{ll}1\ \ &\textrm{if there is a }B\in\rho\ \textrm{s.t.}\{i,j\}\subseteq B \\ 0&\textrm{o.w.}  \end{array} \right.
\end{eqnarray}
In other words, a loss of one is incurred each time either: (1) two points are assumed to be in the same cluster when they should not, or (2) two points are assumed to be in different clusters when they should be in the same cluster.

---

The rand loss has several problems, motivating other clustering losses such as the adjusted rand index, but we will look at the rand loss here since the derivation of the Bayes estimator is easy for that particular loss.

***Give references found by Andy***

In the case of density estimation, if the task is to reconstruct the density itself, examples of loss functions include the Hellinger and KL losses.  However, density estimation is usually an intermediate step for another task, and the loss should then be defined on this final task rather than on the intermediate density estimation task.  For example, in the first part of the second assignment, the task under consideration is to construct a visualization of the predictive density using a large but finite number of points $\tilde y\_j\in S$ on a grid $S$.  In this case, an example of loss function is:
\begin{eqnarray} 
L(p,p') = || p\_j - p'\_j ||\_2^2,
\end{eqnarray}
where $p,p'$ are $|S|-$dimensional vectors of nonnegative real numbers (where $p\_j$ corresponds to the density sampled at $\tilde y\_j$).

#### Combining probability models and loss functions to do inference

**Warning:** in this lecture, we use a slightly different notation for the cluster indicators and observations of a DP. We make this change to avoid using the letter $Z$ for random variables, as it is often use in the sampling literature as a normalization constant.

- We now use the random vector $X$ to represent the cluster indicators.
- We now use the random vector $Y$ to represent the observations.  

As reviewed earlier, the Bayesian framework is reductionist: given a loss function $L$ and a probability model $(X, Y) \sim \P$, it prescribes the following estimator:
\begin{eqnarray} 
\argmin\_{x} \E[L(x, X) | Y].
\end{eqnarray}

We will revisit the examples of clustering and density estimation with the loss function defined in Section~\ref{sec:losses} to see how this abstract quantity can be computed or approximated in practice.

First, for the rand loss, we can write:
\\begin{eqnarray}
\argmin\_{\textrm{partition }\rho} \E\left[\randindex(x, \rho)|y\right] & = &
\argmin\_{\textrm{partition }\rho} \sum\_{i<j} \E \left[\1 \left[\1[x\_i = x\_j]\right] \neq \rho\_{ij}|y\right] \\\\
&=&\argmin\_{\textrm{partition }\rho} \sum\_{i<j} \left\{(1-\rho\_{ij})\P(x\_i= x\_j|y) + \rho\_{ij} \left(1- \P(x\_i = x\_j |y)\right)\right\}
\\end{eqnarray}
where $\rho\_{i,j} = (i \sim\_{\rho} j)$.

This means that computing an optimal bipartition of the data into two clusters can be done in two steps:

1. Simulating a Markov chain, and use the samples to estimate $\partstrength\_{i,j} = \P(x\_i = x\_j | y)$ via  Monte Carlo averages.
2. Minimize the linear objective function $\sum\_{i<j} \left\{(1-\rho\_{ij})\partstrength\_{i,j} + \rho\_{ij} \left(1- \partstrength\_{i,j}\right)\right\}$ over bipartitions $\rho$.

Note that the second step can be efficiently computed using min-flow/max-cut algorithms (understanding how this algorithm works is outside of the scope of this lecture, but if you are curious, see \cite{clrs}).  Our focus will be on computing the first step, i.e. the posterior over the random cluster membership variables $x\_i$.  Note that the $\partstrength\_{i,j}$ are easy to compute from samples since the Monte carlo average of a function $f$ applied to MCMC samples converges to the expectation of the function under the stationary distribution (as long as $f$ is integrable, which is the case here since the indicator function is bounded).  Sampling will be the topic of the next section.

For density estimation, going over the same process for the special loss defined Section~\ref{sec:losses}, we get:
\begin{eqnarray} 
\min\_{p\in [0, \infty)^{|S|}} \E \left[ \sum\_{j=1}^{|S|} (p\_j - F(\tilde y\_j))^2 \Big| y\right] =  
 \sum\_{j=1}^{|S|} \min\_{p\_j\in [0, \infty)} \E \left[(p\_j - F(\tilde y\_j))^2 \Big| y\right],
\end{eqnarray}
so that the estimator is $\hat p\_j = \E[ F(\tilde y\_j) | y]$, which we will approximate using again a MC average:
\begin{eqnarray}\label{eq:clust-loss-min}
\hat p\_j &=& \E[ F(\tilde y\_j) | y]\notag \\\\
&\approx& \frac{1}{T} \sum\_{t = 1}^{T} \E[ F(\tilde y\_j) | x^{(t)}, y].
\end{eqnarray}
We will come back on how to compute $\E[ F(\tilde y\_j) | x^{(t)}, y]$ for each sample $x^{(t)}$ later on, but the important point here is that again, we need to simulate cluster membership variables $x$ from the posterior distribution.

#### Sampling

At each iteration, the collapsed sampler maintains values only for the cluster membership variables $x$, or more precisely, a labeled partition $\rho$ over the datapoints, which, as will is see, is sufficient thanks to the results on the Chinese Restaurant Process representation of the part~2 of this set of notes.   We will write $(\rho(x) = \rho)$ for the labeled partition induced by the cluster membership variables (overloading $\rho(\cdot)$ to denote also the function that extracts the labeled partition induced by the cluster membership variables).

***Recall/review of CRP and general sampling algo, connecting the two***

***Simplification of the formula***

***Note on Peskun paper***


### Another DPM posterior sampling method


<!-- Hierarchical DP? Keep for later? -->


### Supplementary references and notes

**Under construction**
