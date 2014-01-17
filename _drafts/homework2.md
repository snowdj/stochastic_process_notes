---
layout: post
title: "Exercise 2: Bayesian background"
category: 'Homework'
---

While you start familiarizing yourself with java in preparation for the third exercise, in this exercise we will explore the basics of a simple language specialized to Bayesian inference called JAGS.

JAGS (or similar languages, such as Stan and WinBugs) and Java (or similar languages, such as C++, C# or Go) form a good combo for approaching Bayesian inference problems: the former for small problem and prototyping, and the latter for large scale and non-standard models.

### Optional readings

Start by familiarizing yourself with the main features of the JAGS modelling language. See the [manual](http://people.math.aau.dk/~kkb/Undervisning/Bayes13/sorenh/docs/jags_user_manual.pdf), especially 2.1, 2.2, 2.4, 2.5, 3, and have a quick look at 5, 6. If there are any questions, I will devote some time on Monday for Q&A. Make sure to read the JAGS resource section below as well.

Also, in preparation for the next exercise, it might be useful to start reading the following (but be assured that we will still provide more help during the next lab):

- [Learning the Java Language](http://docs.oracle.com/javase/tutorial/java/index.html) (you can skip the following topics on first reading: nested classes, annotations, generics) 
- [Collections](http://docs.oracle.com/javase/tutorial/collections/index.html) (you can skip the following topics on first reading: Aggregate Operations, custom implementation, Interoperability)

### JAGS resources

#### Logistics

You should hand-in the CODA files produced by the MCMC run (see below), archived at the root of a zip file. Some questions may also require a pdf containing a plot and/or some comments and analysis. Use the same file naming convention as for the previous assignment.

You will have two options for solving this exercise:

- Use a web-based tool designed for this course, where you edit the model online and the code is executed on Amazon EC2 (see instructions below).
- Install JAGS locally (open source software available at [http://mcmc-jags.sourceforge.net](http://mcmc-jags.sourceforge.net)). Even if you choose this option, have a quick look at the web-based tool as it contains some skeleton files that will help you get started.

#### Web-based JAGS

1. Access this [URL](http://54.200.129.218/public_models/1). Use chrome with a large window size if possible. ***TODO: CHANGE THIS*** 
   - Username: ``testUser1``
   - Password: ``password``
2. Click on the file on the left that you want to edit.
3. Write down the model in ``clustering.jags`` and the script controlling the number of iterations, variables to monitor, etc., in ``clustering.txt`` **WARNING:** make sure to press ``commit`` regularly and before switching files to edit.
4. To run:
   - commit your changes
   - you need to select the jags script you want to run. 
   - click on ``clustering.jags`` and then:
   - click on the star located besides the filename. 
   - click on ``Run``; this may take some time
   - once the execution is complete, you can press ``Download`` to get the output in a zip file.
   
#### Reading CODA files

Whether you use the web-based tool or a local install, JAGS uses the CODA format to store samples extracted from an MCMC run. You can use the R libraries ``coda`` and ``mcmcplots`` to interpret a CODA file. For example, you can do the following to construct density estimates from an MCMC output stored in ``CODAchain1.txt`` and ``CODAindex.txt`` of the working directory:

```r
library(coda)

pdf('analysis1.pdf')
res = read.coda("CODAchain1.txt", "CODAindex.txt")
plot(res)

dev.off()
```

#### Misc 

Some potential traps:

- JAGS can have problems with observing deterministic nodes. You can always put a peaked normal around the value to sidestep that issue.
- While the syntax of R and JAGS are generally fairly similar, note one important difference: JAGS uses a *precision* parameter for the spread of normal distributions (the precision is the inverse of the variance), while R uses the standard deviation parameterization.

### Question 1: Truncated DP posterior simulation in JAGS

#### Context

<img src="{{ site.url }}/data/geyser/geyser-data.jpg" alt="Drawing" style="width: 300px; float: right"/> 
In this question, you will use a truncated DP model to analyze the classical [Old Faithful Geyser Data](http://www.stat.cmu.edu/~larry/all-of-statistics/=data/faithful.dat). The normalized version can be downloaded [here]({{ site.url }}/data/geyser/question1.data).  (side question: why do you think normalization useful in this context?)

#### Goal

The dataset contains a collection of pairs. Each pair contains an eruption time in minutes, a a waiting time to the next eruption. Note that we have held-out some values. The goal of this exercise will be to find the posterior distributions of the missing values via a truncated DP (lecture of Jan 15). Your answers should therefore consists in two file: a zip with the CODA samples, and a pdf containing density estimates of the approximate posteriors.

The model should have the following characteristics:

- a truncation $K = 5$
- $\alpha\_0 = 10$, and 
- a normal likelihood model with priors of your choice (note that conjugacy is not required here)

### Question 2: A simple dead reckoning problem

#### Context

This question is inspired by a mammal tracking application recently presented a department seminar. The problem consists in a vastly simplified dead reckoning problem. Suppose we have a marine mammal tagged with a tracking device measuring GPS locations and acceleration. The GPS can give accurate, absolute positions, but when the animal dives, only the accelerometer is operational, giving noisy acceleration estimates. <img src="{{ site.url }}/data/geyser/trajectory.jpg" alt="Drawing" style="width: 300px; float: right"/> 

The objective is to reconstruct a trajectory given only its end-points and noisy intermediate accelerations. At the same time, we want to assess how reliable these estimates are.

#### Goal

The true trajectory (number of time steps=20) is shown in the figure on the right  (where the $z$ axis is the depth), but note that your algorithm should only use the endpoints of this strategy. Your goal is to approximate the posterior density of the depth at time step 3, 10, and 19.

We will use the following simple kinetics model (see the the file located [here]({{ site.url }}/data/trajectory/question2.data) for actual values of the variables used below):

- we ignore friction
- we treat acceleration as instantaneous pulses, one at each point in the discretization, independently and normally distributed on each axis (with mean and sd given by ``true.accel.mean`` and ``true.accel.sd``)
- we only observed a noisy version of each of the accelerations, assumed to be normally distributed around the truth, with sd given by ``accel.error.sd``
- the time interval between time discretization grid points is assumed to be ``deltat``
- only the initial position, ``true.pos[1,]``, initial velocity, ``true.velo[1,]``, and final position, ``true.pos[len,]``

Your answers should consists in two file: a zip with the CODA samples, and a pdf containing density estimates of the posterior depths at time step 3, 10 and 19. Interpret and comment on whether these seem reasonable.


