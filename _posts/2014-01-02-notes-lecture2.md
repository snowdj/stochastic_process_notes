---
layout: post
title: "Lecture 2: Forward sampling continued"
category: 'Lecture'
---


### More detailed overview via *Forward Simulation* (continued)

**Recall:** Forward simulation is the process of generating a realization (or FDDs of a realization) of a stochastic process from scratch (without conditioning on any event). Synonym: generating a synthetic dataset. In contrast to: 



#### Random measures: motivation and Dirichlet process forward simulation


Based on lecture notes created the first time I taught this course: [PDF](http://www.stat.ubc.ca/~bouchard/courses/stat547-sp2011/notes-part2.pdf).

**Motivation:** Mixture models. Example: modelling the weight of walruses. Observation: weights of specimens $x\_1, x\_2, $. Inferential question: what is the most atypical weight among the samples?

**Method 1:** 

- Find a normal density $\phi\_{\mu, \sigma^2}$ that best fits the data. 
- Bayesian way: treat the unknown quantity $\phi\_{\mu, \sigma^2}$ as random.
- Equivalently: treat the parameters $\theta = (\mu, \sigma^2)$ as random, $(\sigma^2 > 0)$. Let us call the prior on this, $p(\theta)$ (for example, another normal times a gamma, but this will not be important in this first discussion).

Limitation: fails to model sexual dimorphism. Solution: <img src="{{ site.url }}/images/height-hist.jpg" alt="Drawing" style="width: 200px; float: right"/>

**Method 2:** Use a *mixture models*, with two mixture components, each one assumed to be normal.

- In terms of density, this means our modelled density has the form:
\\begin{eqnarray}
\phi = \pi \phi\_{\theta\_1} + (1 - \pi) \phi\_{\theta\_2},
\\end{eqnarray}
for $\pi \in [0, 1]$.
- Equivalently, we can write it in terms of auxiliary random variables $Z\_i$, one of each  associated to each measurement $X\_i$: 
\\begin{eqnarray}
Z\_i & \sim & \Cat(\pi) \\\\
X\_i | Z\_i, \theta\_1, \theta\_2 & \sim & \Norm(\theta\_{Z\_i}).
\\end{eqnarray}
- Each $Z\_i$ can be interpreted as the sex of animal $i$.

Unfortunately, we did not recorded the male/female information when we collected the data!

- Expensive fix: Do the survey again, collecting the male/female information
- Cheaper fix: Let the model guess, for each datapoint, from which cluster (group, mixture component) it comes from.

Since the $Z\_i$ are unknown we need to model a new parameter $\pi$. Equivalently, two numbers $\pi\_1, \pi\_2$ with the constraint that they should be nonnegative and sum to one. Interpretation: the population frequency of each sex. We need to introduce a prior on these unknown quantities.

- Simplest choice: $\pi \sim \Uni(0,1)$. But fails to model our prior belief that male and female frequencies should be close to $50:50$.
- To encourage this, pick a prior density proportional to: <img src="{{ site.url }}/images/beta.jpg" alt="Drawing" style="width: 200px; float: right"/>
\\begin{eqnarray}
p(\pi\_1, \pi\_2) \propto \pi\_1^{\alpha\_1 - 1} \pi\_2^{\alpha\_2 - 1},
\\end{eqnarray}  
where $\alpha\_1 > 0, \alpha\_2 > 0$ are fixed numbers (called hyper-parameters). The $-1$ are just for convenience, to have a simpler restrictions on the hyper-parameters (but this will become very convenient soon!). 
- The hyper-parameters are sometimes denoted $\alpha = \alpha\_1$ and $\beta = \alpha\_2$. To encourage values of $\pi$ close to $1/2$, pick $\alpha = \beta = 2$. To encourage this even more strongly, pick $\alpha = \beta = 20$. To encourage a ratio $r$ different that $1/2$, make $\alpha$ and $\beta$ grow at different rates, with $\alpha/(\alpha+\beta) = r$. 

**Dirichlet distribution:** This generalizes to more than two mixture components easily. If there are $K$ components, the density of the Dirichlet distribution is proportional to:
\\begin{eqnarray}
p(\pi\_1, \dots, \pi\_K) \propto \prod\_{k=1}^{K} \pi\_k^{\alpha\_k - 1},
\\end{eqnarray} 
Note that the normalization of the Dirichlet distribution is analytically tractable using Gamma functions $\Gamma(\cdot)$ (a generalization of $n!$). See the wikipedia [entry](http://en.wikipedia.org/wiki/Dirichlet_distribution) for details.

**Important take home message:** each component $k\in \{1, \dots, K\}$ is associated with a probability $\pi\_k$, and a value $\theta\_k$. This is therefore a discrete measure with atoms at $(\theta\_1, \dots, \theta\_K)$ and probabilities given by the components of the vector $\pi = (\pi\_1, \dots, \pi\_k)$. Formally, if we are given any set $C = A\times B \subset \RR \times (0, \infty)$ (the space $\RR \times (0, \infty)$ comes from the fact that we want $\theta$ to have two component, a mean and a variance, where the latter has to be positive), we can define the following probability mesure
\\begin{eqnarray}
G(C) = \sum\_{k=1}^{K} \pi\_k \1[\theta\_k \in C].
\\end{eqnarray} 

Since the location of the atoms and their probabilities are random, we can say that $G$ is a random measure. The Dirichlet distribution (together with the prior $p(\theta)$), define the distribution of these random discrete distributions.

Another view of a random measure: a collection of real  random variables indexed by sets in a sigma-algebra: $G\_{C} = G(C)$, $C \in \sa$.

Limitation: we have to specify the number of components $K$. The Dirichlet *process* is also a distribution on atomic distributions, but where the number of atoms can be countably infinite.

**Definition:**  A Dirichlet process (DP), $G\_C(\omega) \in [0, 1]$, $C \in \sa$, is specified by:

1. A *base measure* $G\_0$ (this corresponds to the density $p(\theta)$ in the previous example).
2. A *concentration parameter* $\alpha\_0$ (this plays the same role as $\alpha + \beta$ in the simpler Dirichlet distribution).

To do forward simulation of a DP, do the following:

1. Start with a current stick length of one: $s = 1$
2. Initialize an empty list $r = ()$, which will contain atom-probability pairs.
3. Repeat:
   1. Generate a new independent beta random variable: $\beta \sim \Beta(1, \alpha\_0)$.
   2. Create a new stick length: $\pi = s \times \beta$.
   3. Sample a new independent atom from the base distribution: $\theta \sim G_0$.
   4. Add the new atom and its probability to the result: $r \gets r \circ (\theta, \pi)$
   5. Compute the remaining stick length: $s \gets s \times (1 - \beta)$
   
This is not quite an algorithm (it never terminates), but you will see in the first exercise how to transform it into a valide algorithm.

