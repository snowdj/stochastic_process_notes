---
layout: post
title: Lecture 1 notes
category: 'Announcement'
---


Brief Overview
--------------

Based on slides I created the first time I taught this course: [PDF](http://www.stat.ubc.ca/~bouchard/courses/stat547-sp2011/lecture1.pdf).

**What is a stochastic process?** A collection of random variables indexed by an arbitrary  set $S$. This gets interesting when $S$ is uncountable. 

**Topic of this course:** Why/when is it useful to have $S$ uncountable? How can we still do inference using our finite brains and computers?

**Note:** $S$ is not necessarily a subset of the real line! 

**Examples of stochastic processes:**

- Time series: $S = \\RR$. <img src="/images/brownian.jpg" alt="Drawing" style="width: 200px; float: right"/>
 - $Y_s(\omega) \in \\RR$: value of the random variable with index $s\in S$. The variable $\omega$ denotes an outcome.
 - Example: Brownian motion.
 - Alternative view: random function.
 - Applications: 
  - economic/financial indicators, 
  - frequency of a population having a certain genetic mutation, 
  - periodic phenomena: climate, crop statistics, etc.
- Generalizations:
  - $S = \\RR^2$: many applications in spatial statistics. <img src="/images/gp.jpg" alt="Drawing" style="width: 200px; float: right"/>
  - relaxing assumptions on index set: e.g. taking $S$ to be a space with a nice inner product structure (Hilbert space). Example we will cover: Gaussian processes. 
  - flexible classes of diffusions and jump processes: Continuous-time Markov chains, stochastic PDEs. We should be able to cover a bit of both in this course. 
- Random measures: $S = \\sa$, a sigma-algebra.
 - Recall: a measure is a function that takes a set and returns a non-negative number, subject to countable additivity of disjoint sets constraints.
 - Examples: Dirichlet process, Pitman-Yor process.
 - Applications: mixture models (see below).
 - Special case: point processes (integer-valued measure, the number of points that fall in a given set).


**We always have finite observations, why do we need this?** 

- We may not know yet where predictions will be needed.
- Theory predicts that we may be forced to use stochastic processes as latent variables: de Finetti theorem

**What is de Finetti theorem?** TODO


Time series
-----------


Mixture models
--------------

Based on lecture notes created the first time I taught this course: [PDF](http://www.stat.ubc.ca/~bouchard/courses/stat547-sp2011/notes-part2.pdf).

TODO

Logistics
---------

Refer to the website. Some points emphasized in class:

- Warning: the material covered in this course is fairly advanced. Stay on top! But you can get *lots* of help at the labs/office hours. 
- Attendance at labs optional, but strongly recommended.
- Discussion forum: To encourage richer open exchanges, Seong and I will *only* use this platform to answer course-related questions (unless for personal matters). See ``Contact`` link in the top menu.
- Languages used in the exercises: Java, Bugs, R.
 - R is great to use existing statistical methodologies, but not to develop new ones.
 - We will organize remedial tutorials as we go along. Please consider attending if you are not familiar with the covered topic.
- Environments supported: Mac OS X, linux.