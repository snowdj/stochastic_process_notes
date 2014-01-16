---
layout: post
title: "Exercise 2: Bayesian background"
category: 'Homework'
---

While you start familiarizing yourself with java in preparation for the third exercise, in this exercise we will explore the basics of a simple language specialized to Bayesian inference called JAGS.

JAGS (or similar languages, such as Stan and WinBugs) and Java (or similar languages, such as C++, C# or Go) form a good combo for approaching Bayesian inference problems: the former for small problem and prototyping, and the latter for large scale and non-standard models.

### Optional readings

Start by familiarizing yourself with the main features of the JAGS modelling language. See the [manual](http://people.math.aau.dk/~kkb/Undervisning/Bayes13/sorenh/docs/jags_user_manual.pdf), especially 2.1, 2.2, 2.4, 2.5, 3, and have a quick look at 5, 6. If there are any questions, I will devote some time on Monday for Q&A. 

In preparation for the next exercise, it might be useful to read the following (but be assured that we will still provide more extensive help at the next lab):

- [Learning the Java Language](http://docs.oracle.com/javase/tutorial/java/index.html) (you can skip the following topics on first reading: nested classes, annotations, generics) 
- [Collections](http://docs.oracle.com/javase/tutorial/collections/index.html) (you can skip the following topics on first reading: Aggregate Operations, custom implementation, Interoperability)

### Question 1: Truncated DP posterior simulation in JAGS

#### Context

<img src="{{ site.url }}/data/geyser/geyser-data.jpg" alt="Drawing" style="width: 200px; float: right"/> 
In this question, you will use a truncated DP model to analyze the classical [Old Faithful Geyser Data](http://www.stat.cmu.edu/~larry/all-of-statistics/=data/faithful.dat). The normalized version (why is normalization useful in this context?) can be downloaded [here]({{ site.url }}/data/geyser/question1.data).

#### Objective

The dataset contains a collection of pairs. Each pair contains an eruption time in minutes, a a waiting time to the next eruption. Note that we have held-out some values. The goal of this exercise will be to find the posterior distributions of the missing values via a truncated DP (lecture of Jan 15).

The model should have the following characteristics:

- a truncation $K = 5$
- $\alpha\_0 = 10$, and 
- a normal likelihood model with priors of your choice (note that conjugacy is not required here)

#### Logistics

You should hand-in the CODA files produced by the MCMC run (see JAGS manual for explanations of CODA).

You will have two options for solving this exercise:

- Use a web-based tool designed for this course, where you edit the model online and the code is executed on Amazon EC2 (see instructions below).
- Install JAGS locally (open source software available at [http://mcmc-jags.sourceforge.net](http://mcmc-jags.sourceforge.net)).

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
   