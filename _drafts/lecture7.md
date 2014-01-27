---
layout: post
title: "Lecture 7: Related topics in computational statistics"
category: 'Lecture'
---

Instructor: Alexandre Bouchard-C&ocirc;t&eacute;

Editor: TBA

You should all be familiar with hard inference problems of the form
\begin{eqnarray}\label{eq:problem}
\pi(x) = \frac{\tilde \pi(x)}{Z},
\end{eqnarray}
where the goal is to compute an expectation with respect to $\pi$ and/or estimate $Z$.

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
- HMC if time permits





### Supplementary references and notes

**Under construction**
