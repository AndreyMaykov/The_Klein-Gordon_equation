# The Klein-Gordon equation
In this project, a code for numerically solving the initial-boundary problem
![ ](https://github.com/AndreiMaikov/The_Klein-Gordon_equation-1/blob/main/img/ibp_2.png)
was developed. 

This problem is related to the Klein-Gordon equation (https://en.wikipedia.org/wiki/Numerical_analysis) and is important for testing the so-called open boundary conditions (also known as transparent or non-reflecting boundary conditions) for equations that describe electromagnetic or acoustic waves in various media.

The code implements a family of explicit and implicit finite difference methods (https://en.wikipedia.org/wiki/Finite_difference_method). For the latter, an iterative algorithm was utilized to obtain the numeric solution. 

An executable compiled with a Visual C++ compiler was integrated via MathLink with a program complex built using Wolfram Mathematica for studying open boundary conditions.

**Note.** The code was (and still is) documented very scarcely, as it was only meant to be used by the author within a short time. I am going to add more detailed documentation later on.
