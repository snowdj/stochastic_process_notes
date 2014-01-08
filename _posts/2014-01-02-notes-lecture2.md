---
layout: post
title: "Lecture 2: Forward sampling continued"
category: 'Lecture'
---


### More detailed overview via *Forward Simulation* (continued)

**Recall:** Forward simulation is the process of generating a realization (or FDDs of a realization) of a stochastic process from scratch (without conditioning on any event). Synonym: generating a synthetic dataset. In contrast to: 

**Another synonym:** generative process/story.

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
where $\alpha\_1 > 0, \alpha\_2 > 0$ are fixed numbers (called hyper-parameters). The $-1$'s give us a simpler restrictions on the hyper-parameters required to ensure finite normalization of the above expression). 
- The hyper-parameters are sometimes denoted $\alpha = \alpha\_1$ and $\beta = \alpha\_2$. 
  - To encourage values of $\pi$ close to $1/2$, pick $\alpha = \beta = 2$. 
  - To encourage this even more strongly, pick $\alpha = \beta = 20$. (and vice versa, one can take value close to zero to encourage realizations with one point mass larger than the other.)
  - To encourage a ratio $r$ different that $1/2$, make $\alpha$ and $\beta$ grow at different rates, with $\alpha/(\alpha+\beta) = r$. But in order to be able to build a bridge between Dirichlet distribution and Dirichlet processes, we will ask that $\alpha = \beta$ <sub>(comment that should be ignored in first reading: one way to obtain the DP from the Dirichlet, which we will not cover in this lecture, see wiki [page](http://en.wikipedia.org/wiki/Dirichlet_process),  is as a weak limit of *symmetric* Dirichlet distributions, i.e. a generalization of the statement $\alpha = \beta$. It may be interesting to ask what happens in the non-symmetric case. This may seem non-exchangeable a priori, but is this addressed when the atoms are sampled iid.).</sub>

**Dirichlet distribution:** This generalizes to more than two mixture components easily. If there are $K$ components, the density of the Dirichlet distribution is proportional to:
\\begin{eqnarray}
p(\pi\_1, \dots, \pi\_K) \propto \prod\_{k=1}^{K} \pi\_k^{\alpha\_k - 1},
\\end{eqnarray} 
Note that the normalization of the Dirichlet distribution is analytically tractable using Gamma functions $\Gamma(\cdot)$ (a generalization of $n!$). See the wikipedia [entry](http://en.wikipedia.org/wiki/Dirichlet_distribution) for details.

**Important points:** 

- Each component $k\in \{1, \dots, K\}$ is associated with a probability $\pi\_k$, and a value $\theta\_k$. 
- This is therefore a discrete measure with atoms at $(\theta\_1, \dots, \theta\_K)$ and probabilities given by the components of the vector $\pi = (\pi\_1, \dots, \pi\_k)$.

**Formally:**

- let us say we are given a set  $A \subset \RR \times (0, \infty)$  <img src="{{ site.url }}/images/setA.jpg" alt="Drawing" style="width: 200px; float: right"/>
  - the space $\RR \times (0, \infty)$ comes from the fact that we want $\theta$ to have two component, a mean and a variance, where the latter has to be positive
  - let us call this space $\Omega = \RR \times (0, \infty)$
- we can define the following probability mesure <img src="{{ site.url }}/images/dirichletRealization.jpg" alt="Drawing" style="width: 200px; float: right"/>
\\begin{eqnarray}
G(A) = \sum\_{k=1}^{K} \pi\_k \1(\theta\_k \in A),
\\end{eqnarray} 
where the notation $\1(\textrm{some logical statement})$ denotes the indicator variable on the event that the boolean logical statement is true, in other words, $\1$ is a random variable that is equal to one if the input outcome satisfies the statement, and zero otherwise.

**Random measure:** Since the location of the atoms and their probabilities are random, we can say that $G$ is a random measure. The Dirichlet distribution (together with the prior $p(\theta)$), define the distribution of these random discrete distributions.

**Connection to previous lecture:** Another view of a random measure: a collection of real  random variables indexed by sets in a sigma-algebra: $G\_{A} = G(A)$, $A \in \sa$.

**Limitation of our setup so far:** we have to specify the number of components/atoms $K$. The Dirichlet *process* is also a distribution on atomic distributions, but where the number of atoms can be countably infinite.

**Definition:**  A Dirichlet process (DP), $G\_A(\omega) \in [0, 1]$, $A \in \sa$, is specified by:

1. A *base measure* $G\_0$ (this corresponds to the density $p(\theta)$ in the previous example).
2. A *concentration parameter* $\alpha\_0$ (this plays the same role as $\alpha + \beta$ in the simpler Dirichlet distribution).

To do forward simulation of a DP, do the following: <img src="{{ site.url }}/images/dpSimulation.jpg" alt="Drawing" style="width: 500px; float: right"/>

1. Start with a current stick length of one: $s = 1$
2. Initialize an empty list $r = ()$, which will contain atom-probability pairs.
3. Repeat, for $k = 1, 2, \dots$
   1. Generate a new independent beta random variable: $\beta \sim \Beta(1, \alpha\_0)$.
   2. Create a new stick length: $\pi\_k = s \times \beta$.
   3. Compute the remaining stick length: $s \gets s \times (1 - \beta)$
   4. Sample a new independent atom from the base distribution: $\theta\_k \sim G_0$.
   5. Add the new atom and its probability to the result: $r \gets r \circ (\theta\_k, \pi\_k)$

Finally, return:
\\begin{eqnarray}
G(A) = \sum\_{k=1}^{\infty} \pi\_k \1(\theta\_k \in A),
\\end{eqnarray} 

This is not quite an algorithm (it never terminates), but you will see in the first exercise how to transform it into a valid algorithm.

#### Sampling from a sampled measure

Since the realization of a DP is a probability distribution, we can sample from this realization! Let us call these samples from the DP sample $\underline{\theta\_1}, \underline{\theta\_2}, \dots$

\\begin{eqnarray}
G & \sim & \DP \\\\
\underline{\theta\_i} | G & \sim & G
\\end{eqnarray}

**Preliminary observation:** If I sample twice from an atomic distribution, there is positive probability that I get two identical copies of the same point (an atomic measure $\mu$ is one that assign a positive mass to a point, i.e. there is a $x$ such that $\mu(\{x\}) > 0$). This phenomenon does not occur with non-atomic distribution (with probability one).

The realizations from the random variables $\underline{\theta\_i}$ live in the same space $\Omega$ as those from $\theta\_i$, but $\underline{\theta}$ and $\theta$ have important differences:

- The list $\underline{\theta\_1}, \underline{\theta\_2}, \dots$ will have duplicates with probability one (why?), while $\theta\_1, \theta\_2, \dots$ generally does not have duplicates (as long as $G\_0$ is non-atomic, **which we will assume today**). 
- Each value taken by a $\underline{\theta\_i}$ correspond to a value taken by a $\theta$, but not vice-versa. 

To differentiate the two, I will use the following terminology:

- *Type:* the variables $\theta\_i$
- *Token:* the variables $\underline{\theta\_i}$

**Question:** How to simulate four tokens, $\underline{\theta\_1}, \dots,  \underline{\theta\_4}$?

1. Hard way: use the stick breaking representation to simulate all the types and their probabilities. This gives us a realization of the random probability measure $G$. Then, sample four times from $G$.
2. Easier way: sample $\underline{\theta\_1}, \dots, \underline{\theta\_4}$ from their marginal distribution directly. This is done via the *Chinese Restaurant process* (CRP).

Let us look at this second method. Again, we will focus on the algorithmic picture, and cover the theory after you have implemented these algorithms.

The first idea is to break the marginal $(\underline{\theta\_1}, \dots, \underline{\theta\_4})$ into sequential decisions. This is done using chain rule: 

1. sample $\underline{\theta\_1}$, then, 
2. sample $\underline{\theta\_2}|\underline{\theta\_1}$, 
3. etc, 
4. until $\underline{\theta\_4}|(\underline{\theta\_1}, \underline{\theta\_2}, \underline{\theta\_3})$. <img src="{{ site.url }}/images/marginalization.jpg" alt="Drawing" style="width: 100px; float: right"/>

The first step is easy: just sample from $G\_0$! For the other ones, we will need to keep track of *token multiplicities*: organize the tokens $\underline{\theta\_i}$ into groups taking so that two tokens are in the same group if and only they picked the same type. Each group is called a *table*, the points in this group are called *customers*, and the value shared by each table, its *dish*.

Once we have these datastructures, the conditional steps 2-4 can be sampled easily using the following decision diagram: 

<img src="{{ site.url }}/images/crp-decisions.jpg" alt="Drawing" style="width: 300px"/>

Formally:

\\begin{eqnarray}
\P(\underline{\theta\_{n+1}} \in A | \underline{\theta\_{1}}, \dots, \underline{\theta\_{n}}) = \frac{\alpha\_0}{\alpha\_0 + n} G\_0(A) + \frac{1}{\alpha\_0 + n} \sum\_{i = 1}^n \1(\underline{\theta\_{i}} \in A).
\\end{eqnarray}

**Clustering/partition view:** This generative story suggests another way of simulating the tokens:

1. Generate a partition of the data points (each block is a cluster)
2. Once this is done, sample one dish for each table independently from $G\_0$.

By a slight abuse of notation, we will also call and denote the distribution on the partition the CRP. It will be clear from the context if the input of $\CRP$ is a partition or a product space $\Omega \times \dots \times \Omega$.

**Important property: exchangeability.** We claim that for any permutation $\sigma : \\{1, \dots, n\\} \to \\{1, \dots, n\\}$, we have:

\\begin{eqnarray}
(\underline{\theta\_1}, \dots, \underline{\theta\_n}) \deq (\underline{\theta\_{\sigma(1)}}, \dots, \underline{\theta\_{\sigma(n)}})
\\end{eqnarray}

You will prove this property in the first exercise. In the meanwhile, you can convince yourself with the following small example, where we compute a joint probability of observing the partition $\{\{1,2,3\},\{4,5\}\}$ with two different orders. First, with the order $1 \to 2 \to 3 \to 4 \to 5$, and $\alpha\_0 = 1$:

<img src="{{ site.url }}/images/CRP-order0.jpg" alt="Drawing" style="width: 450px"/>

Then, with the order $4 \to 5 \to 3 \to 2 \to 1$, we get:

<img src="{{ site.url }}/images/CRP-order1.jpg" alt="Drawing" style="width: 450px"/>

#### FDDs of a DP

Consider the following thought experiment:

1. Let us fix a partition of $\Omega$, $A\_1, A\_2, A\_3$. <img src="{{ site.url }}/images/partition.jpg" alt="Drawing" style="width: 100px; float: right"/> 
2. Imagine that we simulate a realization of $G$. 
3. Let us look at $G(A\_1)$. What is that? Just the sum of the heights of the sticks that fall in $A\_1$. 
4. Do the same for $A\_2, A\_3$. We get a triplet of positive real numbers $(A\_1, A\_2, A\_3)$.
5. Now let's repeat the steps 2-4. We get a list of triplets. Let's ask the following question: what is the multivariate distribution of these triplets? 

**Equivalent definition of a Dirichlet Process:** We say that $G : \sa\_{\Omega'} \to (\sa\_{\Omega} \to [0, 1])$ is distributed according to the Dirichlet process distribution, denoted $G \sim \DP(\alpha\_0, G\_0)$, if for all measurable partitions of $\Omega$, $(A\_1, \dots, A\_m)$, the FDDs are Dirichlet distributed:
\\begin{eqnarray}
(G(A\_1), \dots, G(A\_m)) \sim \Dir(\alpha\_0 G\_0(A\_1), \dots, \alpha\_0 G\_0(A\_m)).
\\end{eqnarray}




#### What's next?

Before going deeper into Bayesian non-parametrics, we will first review some important Bayesian concepts in the parametric setup. This will be the program for next week.


### Supplementary references and notes

**Under construction**


