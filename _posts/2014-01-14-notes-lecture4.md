---
layout: post
title: "Lecture 4: Dirichlet process: inference, properties and extensions"
category: 'Lecture'
---

Instructor: Alexandre Bouchard-C&ocirc;t&eacute;

Editor: Sean Jewell

### Composing many parametric models into a larger, non-parametric model

Based on these [notes](http://www.stat.ubc.ca/~bouchard/courses/stat547-sp2011/notes-part3.pdf) from the previous time I taught the course.

#### Dirichlet process mixture model

We can now integrate these new concepts to make our picture of Bayesian non-parametric mixture models more precise. Let us start with a model based on the stick breaking representation. Later, we will connect it to the CRP representation.

We pick:

- A likelihood model with density $\ell(x|\theta)$ over each individual observation (a weight). For example, a normal distribution (a bit broken since weights are positive, but should suffice for the purpose of exposition).
- A conjugate base measure, $G\_0$, with density $p(\theta)$. As before, $\theta$ is an ordered pair containing a real number (modelling a sub-population mean) and a positive real number (modelling a sub-population variance (or equivalently, precision, the inverse of the variance)). A normal-inverse-gamma distribution is an example of such a prior.
- Some hyper-parameters for this parametric prior, as well as a hyper-parameter $\alpha\_0$ for the Dirichlet prior.

To simulate a dataset, use the following steps:

1. Break a stick $\pi$ according to the algorithm covered [last time]({{ site.url }}/lecture/2014/01/12/notes-lecture3.html).
2. Simulate an infinite sequence $(\theta\_1, \theta\_2, \dots)$ of iid normal-inverse-gamma random variables. The first one corresponds to the first stick segment, the second, to the second stick segment, etc.
3. For each datapoint $i$:
     1. Throw a dart on the stick, use this random stick index $Z\_i$. Grab the corresponding parameter $\theta\_{Z\_i}$.
     2. Simulate a new datapoint $X\_i$ according to $\ell(\cdot | \theta\_{Z\_i})$.
     
This algorithmic description has the following graphical model representation (for details on the plate notation see the [wiki page](https://en.wikipedia.org/wiki/Plate_notation)).  
<img src="{{ site.url }}/images/dpm-gm.jpg" alt="DPMGM" style="width: 300px;"/>
     
#### Truncated DPM posterior simulation
     
While forward simulation is easy, exact posterior simulation is computationally hard, as we will see later in this course. We will therefore need approximate posterior simulation. The good news is that we can show only one additional operation is required: computing the density of a realization of $\pi, Z, X$. 

Depending on the situation, computing the density of a realization can be easy or hard. It is hard when there is a complicated set of outcomes that lead to the same realization. Fortunately, this is not the case here. 

Another complication is that the space is infinite-dimensional, and taking an infinite product of densities is ill-behaved. We will see two elegant ways of resolving this shortly. For now, let us just **truncate** the DP: for some fixed $K$, set $\beta\_K = 1$. This gives probability zero to the subsequent components. We can then express the joint distribution of $\pi, \theta, Z, X$ via the density:

\\begin{eqnarray}\label{eq:joint}
p(\pi, \theta, z, x) = \left(\prod\_{k=1}^K \Beta(\beta\_k(\pi); 1, \alpha\_0) p(\theta\_k) \right) \left( \prod\_{i=1}^N \Mult(z\_i; \pi) \ell(x\_i | \theta\_{z\_i}) \right).
\\end{eqnarray}

Here $\beta\_k(\pi)$ is a deterministic function that returns the $k-th$ beta random variable used in the stick breaking process given a broken stick $\pi$ (convince yourself that these beta variables can be recovered deterministically from $\pi$. *There exists a one-to-one mapping between $\pi$ and $\beta$.*) 

Evaluation of Equation~(\ref{eq:joint}) is the main calculation needed to run a basic MCMC sampler for posterior simulation of a truncated DP.

#### Dish marginalization when conjugacy is available

**Rao-Blackwellization:** as we will see, it is often useful to marginalize some variables and to then define MCMC samplers on this smaller space (obtained through marginalization).

We will show it is possible to marginalize both $\pi$ and $\theta$. Today, we first consider marginalizing only $\theta$.

**Example:** Let us look at an example how this works: 

- Suppose that in a small dataset of size three, the current state of the sampler is set as follows: the current table assignment $z$ consists of 
  - two datapoints at one table, 
  - and one datapoint alone at its table. 
- Set the truncation level $K=4$.  

Let us group the terms in the second parenthesis of Equation~(\ref{eq:joint}) by $\theta\_k$:

\\begin{eqnarray}
p(\pi, \theta, z, x) & = & \left(\prod\_{k=1}^K \Beta(\beta\_k(\pi); 1, \alpha\_0) \right) \left( \prod\_{i=1}^N \Mult(z\_i; \pi) \right) \\\\
& & \times \left( p(\theta\_1) \ell(x\_1 | \theta\_1) \right) \\\\
& & \times \left( p(\theta\_2) \ell(x\_2 | \theta\_2) \ell(x\_3 | \theta\_2) \right) \\\\
& & \times \left( p(\theta\_3)  \right) \\\\
& & \times \left( p(\theta\_4)  \right).
\\end{eqnarray}

To get the marginal over $(\pi, Z, X)$, we need to integrate over $\theta\_1, \dots, \theta\_4$. For tables with nobody sitting at them, we get factors of the form:

\\begin{eqnarray}
\int   p(\theta\_3)   \ud \theta\_3 = 1.
\\end{eqnarray}

For the occupied tables, we get a product of  factors of the form:

\\begin{eqnarray}
\int   p(\theta\_2) \ell(x\_2 | \theta\_2) \ell(x\_3 | \theta\_2) \ud \theta\_2.
\\end{eqnarray}

But this is just $m(x\_2, x\_3)$, which we can handle by conjugacy! Recall from last [lecture]({{ site.url }}/lecture/2014/01/12/notes-lecture3.html): 

\\begin{eqnarray}
m(x) & = & \frac{p\_{h}(z) \ell(x | z)}{p(z | x)} \\\\
& = & \frac{p\_{h}(z) \ell(x | z)}{p\_{u(x, h)}(z)}.
\\end{eqnarray}

Here $x$ is a subset of the dataset, given by the subset of the points that share the same  table. 

**Notation:** Let us introduce some notation for the subset of people that share the same table. 

- Note first that the assignments $z$ induce a partition $\part(z)$ of the customer indices (datapoint indices $\\{1, 2, \dots, N\\}$). 
- Each block $B \in \part(z)$ then corresponds to the subset of *indices* of the customers sitting at a table. 
- Finally, we will write $x\_B$ for the subset of the observations themselves; for example in the above example $\part(z) = \\{\\{1\\}, \\{2,3\\}\\}$ and if $B = \\{2, 3\\}, x\_B = (x\_2, x\_3)$.

This gives us:

\\begin{eqnarray}
\int  p(\pi, \theta, z, x) \ud \theta & = & \left(\prod\_{k=1}^K \Beta(\beta\_k(\pi); 1, \alpha\_0) \right) \left( \prod\_{i=1}^N \Mult(z\_i; \pi) \right) \times \prod\_{B \in \part(z)} \left( \int p(\theta) \prod\_{i\in B} \ell(x\_i | \theta) \ud \theta  \right) \\\\
& = &\left(\prod\_{k=1}^K \Beta(\beta\_k(\pi); 1, \alpha\_0) \right) \left( \prod\_{i=1}^N \Mult(z\_i; \pi) \right)  \times \prod\_{B \in \part(z)} m(x\_B)
\\end{eqnarray}

#### Stick marginalization using CRP

We will now push this marginalization further, marginalizing over both $\theta$ and $\pi$ this time. In fact, the only latent variable we will keep around is $\part(Z)$, the partition induced by the repeated DP sampling. This allows us to get rid of the truncation at the same time.

In order to formalize this idea, we will avoid the density notation and work with Lebesgue integrals with respect to measures. Let us  denote the measure on sticks induced by the stick breaking process by $\\gem$. Assume for simplicity that the observations are discrete.

We want to compute, for any partition $\part = \\{B\_1, \dots, B\_K\\}$ of $\\{1, 2, \dots, N\\}$:

\\begin{eqnarray}
\P(\part(Z) = \part, X = x) & = & \int \int \P(\part(Z) = \part|\pi) \left( \prod\_{k=1}^K \prod\_{i \in B\_k} \ell(x\_i | \theta\_k) \right) G\_0^\infty(\ud \theta) \gem(\ud \pi; \alpha\_0) \\\\
& = & \label{eq:marg} \left( \int  \P(Z = z|\pi) \gem(\ud \pi; \alpha\_0) \right)  \left( \int \left( \prod\_{i=1}^N \ell( x\_i | \theta\_{z\_i}) \right) G_0^\infty(\ud \theta) \right) 
\\end{eqnarray}

We showed in the previous section how to compute the second parenthesis of Equation~(\ref{eq:marg}). 

For the first parenthesis of the same equation,  we need to compute the probability of observing a seating arrangement integrating over the stick lengths. We have discussed [last lecture]({{ site.url }}/lecture/2014/01/12/notes-lecture3.html) that this is given by the CRP (and we will explain why this is the case in the following section). Let us denote the probability of this seating arrangement (marginalizing over sticks) by $\CRP(\part(z))$. Recall that this probability can be computed by creating the partition one point at the time (in any order), multiplying conditional CRP probabilities. 

Combining these two facts yields:

\\begin{eqnarray}
\P(\part(Z) = \part, X = x) & = & \label{eq:crp-joint} \CRP(\part(z); \alpha\_0) \left( \prod\_{B \in \part(z)} m(x\_B) \right) 
\\end{eqnarray}

Surprisingly, we managed to get rid of all continuous latent variables! So, in principle, we could compute the exact posterior probability $\P(Z = z | X = x)$ by summing over all partitions $\part$. In practice, the number of partitions is super-exponential, so this is not practical for large problems. As a result we use an approximation algorithm such as MCMC. But the explicit exact formula can be useful on tiny examples, in particular, for debugging approximation algorithms.

### Theoretical properties of DPs

Let us look now at the theoretical justifications of some of the steps we took in the previous argument. First, we will need yet another characterization/definition of the DP.

#### FDDs of a DP

Consider the following thought experiment:

1. Fix a partition of $\Omega$: $A\_1, A\_2, A\_3$. 
<img src="{{ site.url }}/images/partition.jpg" alt="Drawing" style="width: 100px; float: right"/> 
2. Imagine that we simulate a realization of $G$. 
3. Let us look at $G(A\_1)$. What is that? Just the sum of the heights of the sticks that fall in $A\_1$. 
4. Do the same for $A\_2, A\_3$. We get a triplet of positive real numbers $(A\_1, A\_2, A\_3)$.
5. Now let's repeat the steps 2-4. We get a list of triplets. Let's ask the following question: what is the multivariate distribution of these triplets? 

**Equivalent definition of a Dirichlet Process:** We say that $G : \Omega' \to (\sa\_{\Omega} \to [0, 1])$ is distributed according to the Dirichlet process distribution, denoted $G \sim \DP(\alpha\_0, G\_0)$, if for all measurable partitions of $\Omega$, $(A\_1, \dots, A\_m)$, the FDDs are Dirichlet distributed:
\\begin{eqnarray}\label{eq:dp-marg}
(G(A\_1), \dots, G(A\_m)) \sim \Dir(\alpha\_0 G\_0(A\_1), \dots, \alpha\_0 G\_0(A\_m)).
\\end{eqnarray}

**Question:** First, why is this declarative definition well-defined? In other words, why are we guaranteed that there exists a unique process that satisfies Equation~(\ref{eq:dp-marg})?

**Answer:** This follows from the [Kolmogorov extension theorem](http://en.wikipedia.org/wiki/Kolmogorov_extension_theorem). We will revisit this theorem later in more detail, but informally it says that if we can show that the marginal specified by Equation~(\ref{eq:dp-marg}) are *consistent*, then existence and uniqueness of a process is guaranteed by the Kolmogorov extension theorem.

The main step can be illustrated by this example:

---

**Proposition:** 
Let $(A\_1, A\_2)$ and $(B\_1, B\_2, B\_3)$ be the following two measurable partitions of $\Omega$:

<img src="{{ site.url }}/images/lecture-6-9.jpg" alt="Drawing" style="width: 200px;"/>


I.e. $A\_1 = B\_1$, and $(B\_2, B\_3)$ is a partition of $A\_2$.  Let $U\_1, U\_2$ and $V\_1, V\_2, V\_3$ be the random variables as defined in the figure above ($U\_1 = G(A)$, etc.).  Then: $(U\_1, U\_2) \deq (V\_1, V\_2 + V\_3)$, where the special equality symbol denotes equality in distribution.

---

In order to prove this, we use an important tool in the study of Dirichlet processes: gamma representation of Dirichlet distributions:

---

**Lemma:**
If $Y\_1, \dots, Y\_K \sim \gammarv(\alpha\_i, \theta)$ are independent, where $\alpha\_i$ is a shape parameter and $\theta$ is the scale parameter, then:
\begin{eqnarray} 
\left(\frac{Y\_1}{\sum\_k Y\_k}, \dots, \frac{Y\_K}{\sum\_k Y\_k}\right) \sim \dir(\alpha\_1, \dots, \alpha\_K).
\end{eqnarray}

---

**Proof:** A standard change of variable problem.  See the wikipedia page on the [Dirichlet distribution](http://en.wikipedia.org/wiki/Dirichlet_distribution#Related_distributions) and citations there-in.

---

We now turn to the proof of the proposition:

---

**Proof:** Let $\theta > 0$ be an arbitrary scale parameter, and $Y\_i \sim \gammarv(\alpha\_0 G\_0(B\_i), \theta)$ be independent.
From the lemma: 
\\begin{eqnarray} 
(V\_1, V\_2 + V\_3) &\deq \left(\frac{Y\_1}{\sum\_k Y\_k}, \frac{Y\_2 + Y\_3}{\sum\_k Y\_k}\right) \\
&\deq \left(\frac{Y\_1}{\sum\_k Y\_k}, \frac{Y'}{\sum\_k Y\_k}\right) \\
&\deq (U\_1, U\_2),
\\end{eqnarray}
where $Y' \sim \gammarv(G\_0(B\_2) + G\_0(B\_3), \theta) = \gammarv(G\_0(A\_2), \theta)$ by [standard properties](http://en.wikipedia.org/wiki/Gamma_distribution#Summation) of Gamma random variables.

---

The full proof would consider any finite number of blocks, but follows the same argument.    Therefore, we indeed have a stochastic process.

#### Formal connection between the CRP and process definitions

The key step to show this connection is to prove that the Dirichlet process is conjugate to multinomial sampling. This is done by lifting the argument that Dirichlet-multinomial is conjugage via Kolmogorov consistency. 

Let $G\sim \dirp(\alpha\_0, G\_0)$.  

**Recall**:

- since $G$ is a measure-valued random variable, we can sample random variables $\utheta$ from realizations of $G$, i.e. we define $\utheta$ by $\utheta | G \sim G$. <sub>Recall that the notation $\utheta | G \sim G$ means: for all bounded $h$, $\E[h(\utheta) | G] = \int h(x) G(\ud x) = \sum\_{c=1}^\infty \pi\_c h(\theta\_c)$ for $\pi\sim\gem(\alpha\_0)$ and $\theta\_c \sim G\_0$ independent.</sub> 
- we are using the underscore to differentiate the random variables $\theta\_i$ used in the stick breaking construction (types) from the samples $\utheta\_i$ from the Dirichlet process (tokens).   
- $\utheta = \theta\_x$ for a random $x$ sampled from a multinomial distribution with parameters given by the sticks $x\sim\mult(\pi)$.<sub>Note that finite multinomial distributions over $\{1, 2, \dots, K\}$ can be extended to distributions over the infinite list $\{1, 2, \dots\}$, in which case they take as parameters an infinite list of non-negative real numbers that sum to one (e.g.: sticks from the stick breaking construction). We saw in the first exercise how to sample from such generalized multinomials in practice.</sub>

In this section we show that the Dirichlet process is conjugate in the following sense: $G | \utheta \sim \dirp(\alpha'\_0, G'\_0)$ for $\alpha'\_0, G'\_0$ defined below.

To prove this result, we first look at the posterior of the finite dimensional distributions:

---

**Lemma:** 
If $(B\_1, \dots, B\_K)$ is a measurable partition of $\Omega$, then: $$(G(B\_1), \dots, G(B\_K))|\utheta \sim \dir\big(\alpha\_0G\_0(B\_1) + \delta\_{\{\utheta\}}(B\_1), \dots, \alpha\_0G\_0(B\_K) + \delta\_{\{\utheta\}}(B\_K)\big).$$

---

**Proof:**
Define the random variables $Z = (G(B\_1), \dots, G(B\_K))$ and $X = (\delta\_{\utheta(1)}(B\_1), \dots, \delta\_{\utheta(1)}(B\_K)$.  We have:
\\begin{eqnarray}
Z&\sim\dir(\alpha\_0G\_0(B\_1) , \dots, \alpha\_0G\_0(B\_K)) \\ \quad
X | Z &\sim \mult(Z).
\\end{eqnarray}
The result follows by standard Dirichlet-multinomial conjugacy.

---

Since this result is true for all partitions, this means that the posterior is a Dirichlet process as well by the Kolmogorov consistency definition!

We can now obtain the parameters $\alpha'\_0, G'\_0$ of the updated Dirichlet process.  To get $\alpha\_0$, we take the sum of the parameters of any finite dimensional Dirichlet distribution, obtaining $\alpha'\_0 = \alpha\_0 + 1$.  To get $G'\_0$, we normalize the expression in the conclusion of the last lemma to get:
\\begin{eqnarray}
G'\_0 & = & \frac{\alpha\_0}{\alpha\_0 + 1}G\_0 + \frac{1}{\alpha\_0+1} \delta\_{\{\utheta\}} \\\\
\alpha'\_0 & = & \alpha\_0 + 1.
\\end{eqnarray}

This formula can be generalized to the case of multiple observations, by applying it $n$ times:

---

**Proposition:**
Suppose $G \sim \dirp(\alpha\_0, G\_0)$ and $\utheta\_i | G \sim G$ for $i\in\{1, \dots, n\}$, iid given $G$.  Then the posterior has the following distribution:
\\begin{eqnarray}
G | \utheta\_1, \dots, \utheta\_n \sim \dirp\left(\alpha\_0 + n, \frac{\alpha\_0}{\alpha\_0 + n} G\_0 + \frac{1}{\alpha\_0+n} \sum\_{i=1}^n \delta\_{\{\utheta\_i\}}\right).
\\end{eqnarray}

---


#### Equivalence of stick breaking and process definitions

Refer to section 2.5 of these [notes](http://www.stat.ubc.ca/~bouchard/courses/stat547-sp2011/notes-part2.pdf) from the first time I taught the course.



### Supplementary references and notes

- [Markov Chain Sampling Methods for Dirichlet Process Mixture Models](http://www.jstor.org/stable/1390653)
Radford M. Neal
Journal of Computational and Graphical Statistics , Vol. 9, No. 2 (Jun., 2000) , pp. 249-265 

- [Gibbs Sampling Methods for Stick-Breaking Priors](http://www.jstor.org/stable/2670356)
Hemant Ishwaran and Lancelot F. James
Journal of the American Statistical Association , Vol. 96, No. 453 (Mar., 2001) , pp. 161-173