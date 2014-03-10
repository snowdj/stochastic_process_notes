---
layout: post
title: "Lecture 5: More on inference"
category: 'Lecture'
---

Instructor: Alexandre Bouchard-C&ocirc;t&eacute;

Editor: Neil Spencer

***TODO: Office hours scheduling***

### Theoretical properties of DPs (continued)

#### Summary of material from last time

**FDDs of a process:**

  - Finite Dimensional Distribution
  - In the case of a random measure: the measure evaluated at each of the blocks in a partition.
  
**Kolmogorov definition/characterization of DPs:**  a random measure is a Dirichlet process iff the FDDs are Dirichlet distributions ($\star$)

**Conjugacy:** of DP $G$ and (generalized) multinomial sampling $\utheta | G$ (i.e., $G|\utheta$ is a DP with updated hyperparameters). 

Why this is true:

- FDDs are Dirichlet distributions (by $(\Rightarrow)$ of ($\star$))
- Therefore, FDDs of $G|\utheta$ are Dirichlet distributions
- Therefore, $G|\utheta$ is a Dirichlet process (by $(\Leftarrow)$ of ($\star$))

**CRP follows from the same argument:** via the [predictive distribution](http://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution#Conditional_distribution) of a Dirichlet-multinomial model.

#### Moments of a DP

In this section, we derive the  first and second moments of $G(A)$, for:

- $G \sim \dirp(\alpha\_0, G\_0)$ and  
- $A \subset \Omega$ (fixed, deterministic).  

To do this derivation, we use the Kolmogorov definition and consider the partition $(A, A^c)$ of $\Omega$.  

We get:
\\begin{eqnarray}
(G(A), G(A^c)) \sim \mbox{Dir}(\alpha\_0G\_0(A), \alpha\_0G\_0(A^c)).
\\end{eqnarray}
This implies that:
\\begin{eqnarray}
G(A) \sim \mbox{Beta}(x, y),
\\end{eqnarray}
where $x$ denotes $\alpha\_0G\_0(A)$, and $y$ denotes $\alpha\_0G\_0(A^c)$.

**The first moment:** of  $G(A)$ is [therefore:](http://en.wikipedia.org/wiki/Dirichlet_distribution#Moments)
\\begin{eqnarray}
\E[G(A)]=\frac{x}{x+y}=\frac{\alpha\_0G\_0(A)}{\alpha\_0G\_0(A)+\alpha\_0G\_0(A^c)}=G\_0(A), 
\\end{eqnarray}
and **the second moment** of  $G(A)$ is 
\\begin{eqnarray}
\var[G(A)] &=& \frac{xy}{(x+y)^2(1+x+y)} \\\\
&=&\frac{\alpha\_0^2G\_0(A)(1-G\_0(A))}{\alpha\_0^2(\alpha\_0+1)} \\\\
&=&\frac{G\_0(A)(1-G\_0(A))}{\alpha\_0+1}.
\\end{eqnarray}

This gives an interpretation of $\alpha\_0$ as a precision parameter for the Dirichlet process. 

### Inference using DPMs (Dirichlet Process Mixture model): an end-to-end example

**Goals of this section:** 

- show the detail of using a DP to answer a statistical question. 
- We will introduce along the way two MCMC techniques that can be used to approximate the posterior.  
- We will also look at the question of how to use the posterior to answer the task at hand.

**Statistical tasks considered:** 

- **[Density estimation](http://www.foxgo.net/uploads/2/1/3/8/2138775/foxgodensityestimation_for_statistics_and_data.pdf):** The task in density estimation is to give an estimate, based on observed data, of an unobservable underlying probability density function. The unobservable density function is thought of as the density according to which a large population is distributed; the data are usually thought of as a random sample from that population.  
- **[Clustering](http://www.norusis.com/pdf/SPC_v13.pdf):** In cases where the population is thought as being the union of sub-populations, the task of cluster analysis is to find the sub-population structure (usually without labeled data).  Let us assume for simplicity that we wish to separate the data into two clusters. <sub>Note that the Dirichlet process is still a useful tool, even when the number of desired cluster is fixed.  This is because each cluster that is output may need internally more than one mixture to be explained adequately under the likelihood model at hand. Note also that the DPM is not necessarily the right tool for finding a point estimate on the number of cluster components, in contrast to popular belief: see [Miller and Harrison, 2013](http://media.nips.cc/nipsbooks/nipspapers/paper_files/nips26/173.pdf) for a simple cautionary example of inconsistency of the DPM for identifying the number of components.</sub>

**Bayesian approach:** Let us take a Bayesian approach to these problems.  This means that we need to pick:

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
where $(i \sim\_\rho j) = 1$ if there is a $B\in\rho$ s.t. $\\{i,j\\}\subseteq B$, and $(i \sim\_\rho j) = 0$ otherwise.

In other words, a loss of one is incurred each time either: (1) two points are assumed to be in the same cluster when they should not, or (2) two points are assumed to be in different clusters when they should be in the same cluster.

<img src="{{ site.url }}/images/Randloss.jpg" alt="Drawing" style="width: 300px;"/> 

---

The rand loss has several problems, motivating other clustering losses such as the adjusted rand index, but we will look at the rand loss here since the derivation of the Bayes estimator is easy for that particular loss. (See [Fritsch and Ickstadt, 2009](http://ba.stat.cmu.edu/journal/2009/vol04/issue02/fritsch.pdf) for discussion on other losses and how to approach the Bayes estimator optimization problem for these other losses.)

In the case of density estimation, if the task is to reconstruct the density itself, examples of loss functions include the Hellinger and KL losses.  However, density estimation is usually an intermediate step for another task, and the loss should then be defined on this final task rather than on the intermediate density estimation task.  

#### Combining probability models and loss functions to do inference

**Warning:** in this lecture, we use a slightly different notation for the cluster indicators and observations of a DP. We make this change to avoid using the letter $Z$ for random variables, as it is often use in the sampling literature as a normalization constant.

- We now use the random vector $X$ to represent the cluster indicators.
- We now use the random vector $Y$ to represent the observations.  

As reviewed earlier, the Bayesian framework is reductionist: given a loss function $L$ and a probability model $(X, Y) \sim \P$, it prescribes the following estimator:
\begin{eqnarray} 
\argmin\_{x} \E[L(x, X) | Y].
\end{eqnarray}

We will now see with the current example how this abstract quantity can be computed or approximated in practice.

First, for the rand loss, we can write:
\\begin{eqnarray}
\argmin\_{\textrm{partition }\rho} \E\left[\randindex(X, \rho)|Y\right] & = &
\argmin\_{\textrm{partition }\rho} \sum\_{i<j} \E \left[\1 \left[\1[X\_i = X\_j] \neq \rho\_{ij}\right]|Y\right] \\\\
&=&\argmin\_{\textrm{partition }\rho} \sum\_{i<j} \left\\{(1-\rho\_{ij})\P(X\_i= X\_j|Y) + \rho\_{ij} \left(1- \P(X\_i = X\_j |y)\right)\right\\} \label{eq:loss-id}
\\end{eqnarray}
where $\rho\_{i,j} = (i \sim\_{\rho} j)$. 

The above identity comes from the fact that $\rho\_{i,j}$ is either one or zero, so:

- the first term in the the brackets of Equation~(\ref{eq:loss-id}) corresponds to the edges not in the partition $\rho$ (for which we are penalized if the posterior probability of the edge is large), and 
- the second term in the same brackets corresponds to the edges in the partition $\rho$ (for which we are penalized if the posterior probability of the edge is small).

This means that computing an optimal bipartition of the data into two clusters can be done in two steps:

1. Simulating a Markov chain, and use the samples to estimate $\partstrength\_{i,j} = \P(X\_i = X\_j | Y)$ via  Monte Carlo averages.
2. Minimize the linear objective function $\sum\_{i<j} \left\\{(1-\rho\_{ij})\partstrength\_{i,j} + \rho\_{ij} \left(1- \partstrength\_{i,j}\right)\right\\}$ over bipartitions $\rho$.

Note that the second step can be efficiently computed using min-flow/max-cut algorithms (understanding how this algorithm works is outside of the scope of this lecture, but if you are curious, see [CLRS](http://mitpress.mit.edu/books/introduction-algorithms), chapter 26).  

Our focus will be on computing the first step, i.e. the posterior over the random cluster membership variables $x\_i$.  Note that the $\partstrength\_{i,j}$ are easy to compute from samples since the Monte carlo average of a function $f$ applied to MCMC samples converges to the expectation of the function under the stationary distribution (as long as $f$ is integrable, which is the case here since the indicator function is bounded).  Sampling will be the topic of the next section.

#### Sampling via the collapsed representation

At each iteration, the collapsed sampler maintains values only for the cluster membership variables $x$, or more precisely, a labeled partition $\rho$ over the datapoints, which, as we saw last time, is sufficient thanks the Chinese Restaurant Process representation.   

**Recall:**  we write $(\rho(X) = \rho)$ for the labeled partition induced by the cluster membership variables (overloading $\rho(\cdot)$ to denote also the function that extracts the labeled partition induced by the cluster membership variables).

Note that ideally, we would like to know, for each $\rho$, the exact value of the posterior distribution:
\\begin{eqnarray}\label{eq:exact}
\pi(\rho) = \P(\rho(X) = \rho | Y = y) = \frac{\P(\rho(X) = \rho, Y = y)}{\P(Y = y)}.
\\end{eqnarray}
We have shown last week that we can compute the numerator for each partition $\rho$, via the following formula:
\\begin{eqnarray}
\P(\part(X) = \part, Y = y) & = & \CRP(\part; \alpha\_0) \left( \prod\_{B \in \part} m(y\_B) \right).
\\end{eqnarray}

The problem comes for the denominator, which involves summing over all partitions of the data:
\\begin{eqnarray}
Z = \P(Y = y) = \sum\_{\textrm{part} \rho'} \P(\rho(X) = \rho', Y = y).
\\end{eqnarray}

This is why we resort to an approximation of Equation~(\ref{eq:exact}). We show in detail how this problem can be approached via MCMC, removing the need to compute $Z$. MCMC only requires evaluation of the target density $\pi$ up to normalization, a quantity we call $\tilde \pi$ (the numerator of Equation~(\ref{eq:exact})).

---

**Quick review of MCMC:** (if you feel rusty on this, you should still be able to follow today, but come to office hour this week for a refresher).

- General idea: construct a Markov chain with stationary distribution $\pi$.
- The states visited (with multiplicities) provide a consistent approximation of posterior expectations  via the law of large number for Markov chains.
- Start with a proposal $q$, transform it into a transition $T$ that satisfies detailed balance $\pi(x) T(x \to x') = \pi(x') T(x' \to x)$ (which implies $\pi$-invariance) by increasing the self-transition probability (rejection).
- Concretely, accept a state $x'$ proposed from $x$ with probability given by the Metropolis-Hastings ratio:
\begin{eqnarray}\label{eq:mh}
\min\left\\{ 1, \frac{\pi(x')}{\pi(x)} \frac{q(x' \to x)}{q(x \to x')} \right\\} = \min\left\\{ 1, \frac{\tilde \pi(x')}{\tilde \pi(x)} \frac{q(x' \to x)}{q(x \to x')} \right\\}
\end{eqnarray}Importantly, note that the right-hand side sidesteps the difficulty of computing the normalization $Z$.
- To obtain irreducibility, mix or alternate transitions $T^{(i)}$ obtained from a collection of proposals $q^{(i)}$.

---

In the DPM collapsed Gibbs sampler, there will be one proposal $q^{(i)}$, each corresponding to a customer (datapoint). Each is a **Gibbs** move: this means that $q$ is proportional to $\pi$, but with a severely restricted support. 

In our setup, the support is that we only allow one customer $i$ to change table when we are using proposal $q^{(i)}$.

Let us denote the partitions that are neighbors to the current state $\rho$ by $\rho'\_1, \dots, \rho'\_K$. By neighbor, we mean the states that can be reached by changing only the seating of customer $i$.

The Gibbs move simply consists in picking one of the neighbor or outcome $\rho'\_k$ proportionally to the density of the joint:
\begin{eqnarray}\label{eq:naive-way}
p\_k = \crp(\rho'\_k; \alpha\_0) \prod\_{B\in \rho'\_k} m(y\_B).
\end{eqnarray}

<img src="{{ site.url }}/images/Amove.jpg" alt="Drawing" style="width: 300px;"/> 

Some observations:

- Since $K-1$ is the number of tables after removing customer $i$, it is computationally feasible to sample from a multinomial over these $K$ outcomes. In particular, we can easily normalize $p\_1, \dots, p\_K$.
- An important property of Gibbs proposals is that they are never rejected (the acceptance ratio is one). This comes from the fact that $q$ is proportional to $\pi$, creating cancellations in Equation~(\ref{eq:mh}).
- While naively each move would take a running time $O(K^2)$ to perform (because each $p\_k$ involves a product over $O(K)$ tables), this running time can be reduced to $O(K)$. This will be covered in the next exercise. Hint: the $p\_k$s are also proportional (and in fact, equal) to the CRP predictive distributions.


### Another DPM posterior sampling method

While the collapsed sampler is simple to implement and marginalizes several random variables, this sampler is harder to parallelize and distribute (why?), a major concern when scaling DPMs to large datasets. 

The stick breaking representation gives an alternative way to sample that is more amenable to parallelization. It also gives an easy framework to approach non-conjugate models. 

To avoid truncation, we will extend the idea you used for stick-breaking forward simulation (exercise 1, question 2). This will be done using auxiliary variables. 

We first introduce some auxiliary variables $U$, more precisely one for each datapoint:

<img src="{{ site.url }}/images/dp-auxv.jpg" alt="Drawing" style="width: 200px;"/> 

***Warning:*** there is a missing edge between $\theta$ and $y\_n$. Similarly for subsequent figures.

where $U\_i|(X\_i, \pi) \sim \textrm{Uni}(0, \pi\_{X\_i})$. <sub>Note: $\pi$ is used here as the distribution induced by stick, while in the previous section it $\pi$ was used as the target distribution. Since the two are of different types, this will not create ambiguities.</sub>  This yields the following joint probability (imagine an integral on both sides to interpret the $\ud \theta$, etc):
\begin{eqnarray}\label{eq:joint-auxv} 
\P(\ud \theta, \ud \pi, X = x, \ud y, \ud u) &=& G\_0^\infty(\ud \theta) \gem(\ud \pi; \alpha\_0) \P(X = x|\pi) \prod\_{i=1}^n \ell(\ud y\_i | \theta\_{x\_i}) \unif(\ud u\_i; 0, \pi\_{x\_i}).  \\\\
&=& G\_0^\infty(\ud \theta) \gem(\ud\pi; \alpha\_0)  \prod\_{i=1}^n \1[0 \leq u\_i \leq \pi\_{x\_i}] \ud u\_i \ell(dy\_i|\theta\_{x\_i})
\end{eqnarray}
where we used $\P(X\_i = c|\pi) \unif(\ud u\_i; 0, \pi\_{x\_i}) = \1[0 \leq u\_i \leq \pi\_{x\_i}] \ud u\_i$ by the definition of uniform distributions (a uniform on an interval of length $b$ has a normalization of $b$) and $\P(X\_i = c |\pi) = \pi\_c$, also by definition.

The slice sampler proceeds by resampling the variables in three blocks: (1) all the dishes $\theta$, (2) all the cluster membership variables, and (3), both all the stick lengths and all the auxiliary variables.

**Cluster membership variables:** This move resamples the following variables:

<img src="{{ site.url }}/images/auxv-move1.jpg" alt="Drawing" style="width: 200px; float: right"/> 

This can be done by resampling each one  independently (conditionally on the Markov blanket; this follows from [Bayes Ball](http://www.cs.ubc.ca/~murphyk/Bayes/bnintro.html)).

To find the conditional distribution of $X\_i = c$,  for each $c=1, 2, \dots$, we look at the factors in Equation~(\ref{eq:joint-auxv}) that depend on $x\_i$, and obtain:
\begin{eqnarray} 
\P(X\_i = c| \textrm{rest}) \propto \1[0 \le u\_i \le \pi\_c] \ell(\ud y\_i|\theta\_c)
\end{eqnarray}

Thanks to lazy computation, we do not have to instantiate an infinite list of $\pi\_c$ in order to compute this.  Instead, we find the smallest $N \in \ZZ^{+}$, such that $\sum\_{c=1}^{N} \pi\_c > 1- u\_i$. For $c>N$, $\sum\_{c=1}^{N} \pi\_c > 1- u\_i$, so $\1[0 \leq u\_i \leq \pi\_c]=0$.

**Dishes:** This move resamples the following variables:

<img src="{{ site.url }}/images/auxv-move2.jpg" alt="Drawing" style="width: 200px; float: right"/> 

Again, this can be done by resampling each $\theta\_c$  independently (conditionally on the Markov blanket).

The posterior distribution of $\theta\_c|\textrm{rest}$ is
\begin{eqnarray}
\P(\ud \theta\_c| \textrm{rest}) \propto G\_0(\ud\theta\_c) \prod\_{i: x\_i=c} \ell(\ud y\_i|\theta\_c),
\end{eqnarray}
which can be resampled by Gibbs sampling in some conjugate models, and can be resampled using a Metropolis-Hastings step more generally.

**Auxiliary variables and stick lengths:** This resamples both all the auxiliary variables and all the stick lengths in one step.  Using the chain rule, this step is subdivided in two substeps:

<img src="{{ site.url }}/images/auxv-move3.jpg" alt="Drawing" style="width: 500px;"/> 

In other words:
\begin{eqnarray}
\P\left(\ud \pi, \ud u \big| \textrm{rest}\right) =& \P\left(\ud \pi\big| \textrm{ rest except for } u \right) \P\left(\ud u | \textrm{rest} \right)\\
\end{eqnarray}

The second factor can be sampled from readily:
\begin{eqnarray}
\P(\ud u\_i | \textrm{rest}) =& \unif(\ud u\_i; 0, \pi\_{x\_i}).
\end{eqnarray}

To sample from the first factor, we look at an equivalent, **sequential decision** scheme for sampling sticks that is an alternative to dart throwing.  This  sequential scheme makes it easier to see how to resample $\pi$ while integrating out the auxiliary variables.

This alternative consists in visiting the sticks in order, and flipping a coin each time to determine if we are going to pick the current stick, or a stick of larger index:

<img src="{{ site.url }}/images/auxv-move4.jpg" alt="Drawing" style="width: 400px;"/> 

Here the persons represent datapoints, and the right-hand-side represents a decision tree.  Since the $\beta\_c \sim \betarv(1, \alpha\_0)$, and that each decision in the decision tree is multinomial, we get by multinomial-dirichlet conjugacy:
    \begin{eqnarray}
    \P(\ud \pi\_c| \text{ rest except for } u) = \betarv(\ud \pi\_i; a\_c, b\_c)
    \end{eqnarray}
where:
\begin{eqnarray} 
a\_c &=& 1+ \sum\_{i=1}^{n}\1[x\_i=c] \\\\
b\_c &=& \alpha\_0 + \sum\_{i=1}^m \1[x\_i >c]
\end{eqnarray}

#### Comparison

We conclude this section by summarizing the advantages and disadvantages of each methods:

<img src="{{ site.url }}/images/pros-cons.jpg" alt="Drawing" style="width: 300px;"/> 

The restriction on the loss is that the expected loss needs to be computable from only samples from the cluster indicators.  The restriction on the likelihood is the conjugacy assumption discussed in the section on the collapsed sampler.  Note that Rao-Blackwellization does not necessarily mean that the collapsed sampler will be more efficient, since each sampler resamples different blocks of variables with a different computation cost per sample.  

The memory needs of the slice sampler can get large in the case where the value of the auxiliary variables is low.  Note that using a non-uniform distribution on the auxiliary variables could potentially alleviate this problem.  

Note also that for some other prior distributions (for example general stick-breaking distributions, which are covered in the next set of notes), only the slice sampler may be applicable.  In other extensions of the DP, both slice and collapsed samplers are available.



### Supplementary references and notes

[Porteous, I., Ihler, A. T., Smyth, P., & Welling, M. (2012). Gibbs sampling for (coupled) infinite mixture models in the stick breaking representation. arXiv preprint arXiv:1206.6845.](http://arxiv.org/pdf/1206.6845.pdf)
[Kalli, M., Griffin, J. E., & Walker, S. G. (2011). Slice sampling mixture models. Statistics and computing, 21(1), 93-105.](http://kar.kent.ac.uk/24721/1/slice_mix.pdf)

