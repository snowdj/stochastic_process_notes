---
layout: post
title: "Lecture 7: Related topics in computational statistics"
category: 'Lecture'
---

Instructor: Alexandre Bouchard-C&ocirc;t&eacute;

Editor: TBA

### Related topics in computational statistics

By now, you should all be familiar with hard inference problems of the form
\begin{eqnarray}\label{eq:problem}
\bar \pi(x) = \frac{\pi(x)}{Z},
\end{eqnarray}
where the goal is to compute an expectation with respect to the density $\bar \pi$ and/or estimate $Z = \normof{\pi}$. For simplicity, for the first part of today we will assume the space is discrete, so you can think of $\bar \pi$ as a p.m.f.

Our running example will be the posterior of seatings under a CRP.

We will call $\phi$ the function we want to take expectation of:

\begin{eqnarray}
\bar \pi \phi = \sum\_x \phi(x) \bar \pi(x).
\end{eqnarray}

Example of $\phi$?

MCMC is one solution but suffers from some problems:

- Fundamentally sequential, making it harder to adapt to parallel architectures (GPUs, multi-processors, clusters)
- Reversible moves can be difficult to construct in some situations (e.g. requiring reversible jumps)
- The normalization constant $Z$ is non-trivial to estimate from the MCMC output
- Non-trivial to analyze

Today we will discuss other approaches to Equation~(\ref{eq:problem}):

- Sequential Monte Carlo (SMC)
  - Generalized view
  - PMCMC
- Exact propagation of messages
- HMC (if time permits)

#### Sequential Monte Carlo (SMC)

**Background:** The foundation of SMC is importance sampling (IS). IS is a way of building a random measure $\tilde \pi^{(K)}$ such that:

- $\pi^{(K)}$ approximates $\pi$,
- it is easy to compute expectations and normalization of $\pi^{(K)}$,
- these two approximations converge to the true value as $k \to \infty$.

IS works by sampling some samples $X\_i$ according to a proposal with density $q$, and building a random measure with atoms at the $X\_i$s

\begin{eqnarray}
\pi^{(K)} = \sum\_k w^{(k)} \delta\_{X^{(k)}}
\end{eqnarray}

where the random weights are given by $w^{(k)} = w(X^{(k)})$, with: 

\begin{eqnarray}
w(x) = \frac{\pi(x)}{q(x)}
\end{eqnarray}

Note:

1. the normalization $\normof{\pi^{(K)}}$ of this approximation is easy to compute (what is it?)
2. expectations with respect to $\bar \pi^{(K)}$ are easy to compute (what do they look like?)

As long as the support of $q$ is no smaller than the support of $\pi$, we have convergence almost sure (a.s.) of the estimates 1 and 2 above to their targets $Z$ and $\bar \pi \phi$. The key to proving this argument (exercise) is a change of measure:

\begin{eqnarray}
Z &=& \sum\_x \pi(x) \\\\
&=& \sum\_x \left(\frac{\pi(x)}{q(x)}\right) q(x),
\end{eqnarray}

which can be interpreted as $\E[w(X)]$ if $X \sim q$. By the LLN, the approximation

\begin{eqnarray}
\frac{1}{K} \sum_k w\left(X^{(k)}\right)
\end{eqnarray}

will therefore converge to $Z$ if the $X^{(k)}$ are iid samples from $q$.

**Sequential Importance Sampling (SIS):** In the CRP case, it is more natural to write $q$ as a product of conditionals:

\begin{eqnarray}
q(x) = \prod\_n q(x\_{1:n} | x\_{1:n-1}),
\end{eqnarray}

where $x\_{1:n}$ is the seating arrangement of customers 1 through $n$.

We can modify the IS to exploit this structure. This is useful for example in an online setting, where new datapoints come in dynamically.

To do this, we need to consider not only one target distribution $\pi$, but instead several target distributions given by

\begin{eqnarray}
\bar \pi\_n = \frac{\pi\_n}{Z\_n},
\end{eqnarray}

where $\bar \pi\_n$ can be interpreted as the posterior distribution on seating arrangement given datapoints 1 through $n$.

With this in hand, we now simply use a telescoping product to rewrite the weight as a sequence of *weight updates*:

\begin{eqnarray}
\pi\_n &=& \left( \frac{\pi\_1}{\pi\_0} \right) \left( \frac{\pi\_2}{\pi\_1} \right) \dots \left( \frac{\pi\_n}{\pi\_{n-1}} \right) \\\\
w\_n(x\_{1:n}) &=& \frac{\pi\_n(x\_{1:n})}{\pi\_{n-1}(x\_{1:n-1})} \frac{1}{q(x\_{1:n}|x\_{1:n-1})} \\\\
w &=& \prod\_n w\_n \\\\
\end{eqnarray}

where $\pi\_0 \equiv 1$. 

---

**Example:** weight update for the CRP posterior.

A simple proposal consists in using the CRP prior probabilities, $q(x\_{1:n}|x\_{1:n-1}) = \crp(\part(x\_{1:n})|\part(x\_{1:n-1}))$. If we let $\part = \part(x\_{1:n-1})$ and $\part' = \part(x\_{1:n})$, we get:

\begin{eqnarray}
w\_n(x\_{1:n}) &=& \frac{\CRP(\part') \left( \prod\_{B' \in \part'} m(y\_B') \right)}{\crp(\part) \left( \prod\_{B \in \part} m(y\_B) \right)} \frac{1}{\crp(\part(x\_{1:n})|\part(x\_{1:n-1}))} \\\\
&=& \frac{m(y\_{B'\_0})}{m(y\_{B\_0})}.
\end{eqnarray}

Here $B\_0$ and $B'\_0$ is the block of customers before and after insertion of the new customer. Since $B'\_0 = B\_0 \cup \\{n\\}$, we get: $w\_n(x\_{1:n}) = \ell(y\_n | y\_{B\_0})$.

---

Another view of this algorithm is that we have a sequence of approximating measures $\pi^{(K)}\_n$, where approximation $\pi^{(K)}\_n$ is obtained recursively from $\pi^{(K)}\_{n-1}$ by randomly propagating each the atoms via $q(\cdot|\cdot)$, and updating the weights via $w\_n$. Since we have just rewritten the algorithm in an equivalent form it clearly keeps the asymptotic properties (in $K$) from last section.

One can view this as mass propagation/transportation on a tree organized by layers.

**Resampling:** Since we are often less interested in the previous measures $\pi\_{1}, \pi\_{2}, \dots$ and more interested in the latest one, $\pi\_{n}$ it is useful to trade a loss of information in the far past for a more accurate approximation near the present iteration. This is done via *resampling*.

The goal of resampling is to discard the particles of low weight while preserving the asymptotic guarantees outlined in the previous section.

In its simplest form (*multinomial resampling*), it consists in two steps:

1. Normalize the weights $w^{(1)}\_n, \dots, w^{(K)}\_n$ into $\bar w^{(1)}\_n, \dots, \bar w^{(K)}\_n$
2. Resample $K$ times from the categorical distribution with parameters $\bar w^{(1)}\_n, \dots, \bar w^{(K)}\_n$.

Let $\tilde x^{(1)}, \dots, \tilde x^{(K)}$ denote the list of particles sampled from step 2 above (note that there can be repeats). We construct the following new approximation from the $\tilde x^{(k)}$:

\begin{eqnarray}
\tilde \pi^{(K)}\_n = \sum\_k \frac{1}{K} \delta\_{\tilde x^{(k)}}.
\end{eqnarray}

Note that the weights have been all reset to $1/K$. Informally, this is because the random number of repeats now encodes $\bar w^{(K)}\_n$. Formally, it is a good exercise  to show that for all fixed $A$, 

\begin{eqnarray}
\E\left[\tilde \pi^{(K)}(A) | \pi^{(K)}\right] = \pi^{(K)}(A).
\end{eqnarray}

Using this and Minkowsky's inequality, it is not too hard to show that the resampling preserves the asymptotic guarantees as well (hint: prove convergence in L2). In particular, the normalization can now be consistently approximated by:

\begin{eqnarray}
Z^{(K)} = \prod\_n \frac{1}{K} w_n^{(k)}.
\end{eqnarray}

#### Particle MCMC (PMCMC)

Suppose we want to put a prior on and resampling $\alpha\_0$ and/or likelihood hyper-parameters $h$. Different values can induce potentially very different clusterings (especially for $h$). How can we change our SMC algorithm to handle this? This is where PMCMC comes in.

**Simplest PMCMC algorithm: PMMH**

<!-- PMMH (and other PMCMC methods) are essentially MCMC algorithms in which SMC -->

**In construction**

### What is next?

Here is an overview of some other active topics (each with a review/representative/recent paper) in the Bayesian non-parametric literature. We may spend more time on one or two of these, depending on the level of interest. Let me know if you think one of these will be particularly useful for your research.

- [Infinite Hidden Markov Models](http://www.cs.berkeley.edu/~jordan/papers/stickyHDPHMM_LIDS_TR.pdf)
- [Dependent Dirichlet Processes](http://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/steel/steel_homepage/techrep/ddp.pdf)
- [Indian Buffet Process](http://ai.stanford.edu/~tadayuki/papers/miller-phd-dissertation11.pdf)
- [Normalized Random Measures](http://www.stats.ox.ac.uk/~teh/research/npbayes/FavTeh2013a.pdf)
- [Gamma-exponential process](http://books.nips.cc/papers/files/nips24/NIPS2011_1162.pdf)
- [Fragmentation-Coagulation](http://www.stats.ox.ac.uk/~teh/research/npbayes/TehEllBlu2013a.pdf)
- [Random graphs](http://arxiv.org/abs/1401.1137)


### Supplementary references and notes

**Under construction**
