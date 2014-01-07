---
layout: post
title: Lecture 1
category: 'Lecture'
---


### Brief Overview

Based on slides I created the first time I taught this course: [PDF](http://www.stat.ubc.ca/~bouchard/courses/stat547-sp2011/lecture1.pdf).

**What is a stochastic process?** A collection of random variables indexed by an arbitrary  set $S$. This gets interesting when $S$ is uncountable. 

**Topic of this course:** Why/when is it useful to have $S$ uncountable? How can we still do inference using our finite brains and computers?

**Note:** $S$ is not necessarily a subset of the real line! 

#### Examples of stochastic processes: 

**Time series:** $S = \\RR$. <img src="{{ site.url }}/images/brownian.jpg" alt="Drawing" style="width: 200px; float: right"/>
 
- $Y_s(\omega) \in \\RR$: value of the random variable with index $s\in S$. The variable $\omega$ denotes an outcome.
- Example: Brownian motion.
- Alternative view: random function. (**Terminology:** Path space of a stochastic process.)
- flexible classes of diffusions and jump processes (we should be able to cover a bit of both in this course) 
 - Continuous-time Markov chains, 
 - stochastic PDEs. 
- Applications: 
 - economic/financial indicators, 
 - frequency of a population having a certain genetic mutation, 
 - periodic phenomena: climate, crop statistics, etc.

**Natural generalizations:**

- $S = \\RR^2$: many applications in spatial statistics. <img src="{{ site.url }}/images/gp.jpg" alt="Drawing" style="width: 200px; float: right"/>
- $S = $ a space with a nice inner product structure (Hilbert space). Example we will cover: Gaussian processes. 

**Less natural but very useful:** random measures: $S = \\sa$, a sigma-algebra.

- Recall: a measure is a function that takes a set and returns a non-negative number, subject to countable additivity of disjoint sets constraints.
- Examples: Dirichlet process, Pitman-Yor process.
- Applications: mixture models (see below).
- Special case: point processes (integer-valued measure, the number of points that fall in a given set).

#### Motivation 

We always have finite observations, why do we need uncountable spaces?

- We may not know yet where predictions will be needed.
- Theory predicts that we may be forced to use stochastic processes as latent variables: *de Finetti theorem*.

**Terminology:** Finite Dimensional Distributions (FDD)

<!-- **What is de Finetti theorem?** TODO -->


### More detailed overview via *Forward Simulation*

**Forward simulation:** generating a realization (or FDDs of a realization) of a stochastic process from scratch (without conditioning on any event). Synonym: generating a synthetic dataset. In contrast to: 

**Posterior simulation:** obtain samples from the posterior distribution conditioning on some observed parts of the process.

Typically, *forward simulation* is much, much easier. This gives us another motivation for using stochastic processes: even when the inference method is very complicated and cannot be explained to non-experts, it is relatively simple to explain the model using forward simulation. 

#### A simple jump process: a continuous-time Markov chain (CTMC)

**Definition:** A CTMC $Y\_t(\omega)\in\Yscr, t\in\RR$ is specified by: 

1. $\Yscr$, a finite set (the state space, not to be confused with $S$), 
2. a transition probability matrix $M(y, y')$ over $\Yscr$ (we will interpret it as a jump probability, so let's assume there are no self transitions, $M(y,y) = 0$);
3. a positive number $\lambda\_y$ associated to each state $y$, called a rate;
4. a starting point $y_0$.

To do forward simulation of a CTMC, do the following:

1. Initialize $y = y_0$, $p = ()$ to an empty list.
2. Repeat:
   1. Simulate a waiting time: use an exponential distribution with the rate corresponding to the current state, $\Delta t \sim \Exp(\lambda\_y)$.
   2. Add a new segment to the piecewise constant path: $p \gets p \circ (y, \Delta t)$, where $\circ$ denotes concatenation (i.e. we represent a path with a list of pair, where each pair has the length of a segment, $\Delta t$, and its state $y$).
   3. Simulate a transition: use a categorical distribution with the probabilities given by row $y$ of the transition matrix, $y' \sim \Cat(M(y,\cdot))$.
   4. Update the current state: $y \gets y'$

**Note:** you may have seen different definitions. They are equivalent. Sometimes the one above is more convenient (for forward simulation, for example), sometimes, other definition/representations are more convenient (for example, the rate matrix representation is more convenient for obtaining FDDs). This will be an important idea in this course, not just with CTMC but also with all the other processes.

**Applications:** change point detection ($t$ is time), segmentation ($t$ is a position on the genome), survival analysis (add an absorbing state $y\_{\textrm{death}}$ with corresponding rate $\lambda = 0$), molecular biology $\Yscr = \\{\textrm{A,C,G,T}\\}$. We may also cover more advanced examples where $\Yscr$ is countably infinite (examples: strings, counts in birth-death processes, RNA folding dynamics, chemistry ($\Yscr =$ count of different molecules undergoing chemical reactions).


#### A simple diffusion: Brownian motion

**Diffusion:** a continuous time, continuous space stochastic process with continuous paths.

**Important example:** Brownian motion (also known as Wiener process).

**Definition:** A collection of random variables $Y_s(\omega)\in\RR, s\in\RR$ such that:

1. For all,  $s,t \in \RR$, the **increments** are normally distributed:
\\begin{eqnarray}
(Y\_{s+t} - Y\_s) \sim N(0,\sigma^2 t)
\\end{eqnarray}
2. For any two disjoint intervals, $(t\_1, t\_2]$ and $(t\_3, t\_4]$, the increments $(Y\_{t\_2} - Y\_{t\_1})$ and $(Y\_{t\_4} - Y\_{t\_3})$ are independent.
3. The process starts at zero ($Y\_0 = 0$) and its paths are continuous functions of $s$.

**Theoretical question:** Why is this a valid definition? In other words, what guarantees existence and uniqueness? Answer: *Kolmogorov consistency theorem*.

**Practical question:** How to sample from a Brownian motion? Example: using Monte Carlo, approximate the probability that a Brownian motion crosses the line $y=+1$ at least once in the interval $[0, 1/2]$.

**Important idea:** Lazy computing. Often, the question we try to answer does not require the knowledge of all points $Y\_s$ (or at least, can be approximated using a finite number of points). Only instantiate the value at those points $s\_i$ that matter.

**Note:** we may not know a priori what points matter. For example, with the above question may want to dynamically refine the approximation in a region very close to the line $y = +1$.

**Refining the Brownian motion FDDs**: Let's say we have only sampled $Y\_{0.1}$ and $Y\_{0.3}$. We want to add a sample $Y\_{0.2}$ after the fact. How can we do this?

- Note: it can be shown that by characterization 1-3 above, $(Y\_{0.1}, Y\_{0.2}, Y\_{0.3})$ has to be Multivariate Normal.
- Mean: $(0,0,0)$.
- Covariance: say $0 < s < t$, we have:
\\begin{eqnarray}
\\Cov(Y\_s, Y\_t) & = & \E[Y\_s Y\_t] \\\\
             & = & \E[Y\_s (Y\_t - Y\_s + Y\_s)] \\\\
             & = & \E[(Y\_s)^2] + \E[Y\_s - Y\_0] \E[Y\_t - Y\_s] \\\\
             & = & \sigma^2 s + 0. 
\\end{eqnarray}
- Therefore, we can sample from the conditional distribution, $Y\_{0.2} | (Y\_{0.1}, Y\_{0.3})$, which is also Normal&mdash;see Rasmussen's book on Gaussian processes, available online [here](http://www.gaussianprocess.org/gpml/chapters/RWA.pdf).



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


### Logistics

Refer to the website. Some points emphasized in class:

- Warning: the material covered in this course is fairly advanced. Stay on top! But you can get *lots* of help at the labs/office hours. 
- Attendance at labs optional, but strongly recommended.
- Discussion forum: To encourage richer open exchanges, Seong and I will *only* use this platform to answer course-related questions (unless for personal matters). See ``Contact`` link in the top menu.
- Languages used in the exercises: Java, Bugs, R.
 - R is great to use existing statistical methodologies, but not to develop new ones.
 - We will organize remedial tutorials for Java/Bugs as we go along. Please try to drop by if you are not familiar with the covered topic.
- Environments supported: Mac OS X, linux.