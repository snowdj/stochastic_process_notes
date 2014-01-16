---
layout: post
title: "Lab 2: Bugs and java basics"
category: 'Lab'
---

### Computer architecture basics

#### Structure of the execution of a program

When you start a program, the memory is divided into two parts:

- stack(s) : roughly, for local variable
- head : roughly, for global variables

#### Organization of the heap:

Organized into blocks. These blocks can be either:

1. an array, or,
2. an object.

- Each block has an address.
- Each block can contain:
  - **primitives**: int (integers), double (real), char (letters), boolean (booleans)
  - **references** or **pointers** to other blocks (via addresses)
- For array, it's a potentially long list of things of the same kind (retrieved with an integer), for object, it's a small number of things that can be of different kinds (retrieved with a name).
 

#### Organization of the stack:

Also organized into blocks, where each block can contains a mix of primitives and references, that can be accessed by their names. 

Important differences of blocks from the stack (vs. those from the heap):

- In contrast to the heap, the blocks in the stack are organized not only internally, but also relative to each other. They form a list, with the following semantics: last in, first out.
- each block in the stack keeps three additional pieces of info:
	- one is the current line number in your code
	- the local variable needed to compute that line of code
	- the last one is a *scope*. Informally this is a context that will help us access the global variable from the heap that we need. Will come back to that later.

#### Main operations done by a program

- **function call:** Creating a new block in the stack, pausing the execution of the current line of code until the code in the block above is completed.
- **new:** Creating a new block in the heap and getting back its address
   - For arrays, we only need to specify the number of things, and the type of these things. For example, create a new array of 10 integers, and assign it to a variable via: ``int[] myArray = new int[5];``
   - For arrays, we need to specify the name and types of the things the block will contain. This specification is called a **class**.

	Example of a class that contains a pair of real numbers (it should be in a file called ``Coordinate.java``):
	
	```java
	class Coordinate {
		double x;
		double y;
	}
	```
	You can then create a new block containing a real called x and a real called y via ``Coordinate myCoordinate = new Coordinate();``.

#### Scope

- When calling a function, information (mentioned above) is placed in the stack directly above the calling function. 
- When a variable is declared, each variable has associated **scope**. 
- The local variables are accessible only inside the currently executing function and are stored in the stack
- The global variables are stored in the heap and are accessible from the executing function depending on the scope of the function

***Some more topics:***
- Differences with R
- Call by value
- multithreading
- generics
- interfaces
- collection api

### Setting up Java environment

  1. Download [Eclipse](http://www.eclipse.org/downloads/)
  2. Install [Java](http://java.com/en/download/manual.jsp) on your machine if it is not installed already
  3. Start Eclipse, create a new **Java Project** via File->New (or keyboard cmd-N or ctrl-N depending on your O/S) -- name the project ``lab2``
  4. Under the src folder, create a new **Java Class** via File->New. Fill out the form as follows:
  
  	 ```
  	 Package:  noise
  	 Name:     Main
  	 Checkbox: public static void main(String[] args)
  	 ```
  5. Type the following code inside the ``main`` function:
  
	  ``` java
	  System.out.println("Hello, World!")
	  ```
	  
	  Run the program.
  6. Download and setup [Gradle](http://www.gradle.org/)
  
