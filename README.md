Download Link: https://assignmentchef.com/product/solved-fys4150-project-3-2
<br>



<h1></h1>

The task of this project is to integrate first in a brute force manner a sixdimensional integral which is used to determine the ground state correlation energy between two electrons in a helium atom. The integral appears in many quantum mechanical applications. However, if you are not too familiar with quantum mechanics, you can simply look at the mathematical details. We will employ both Gauss-Legendre and Gauss-Laguerre quadrature and Monte-Carlo integration. Furthermore, you will need to parallelize your codes.

We assume that the wave function of each electron can be modelled like the single-particle wave function of an electron in the hydrogen atom. The single-particle wave function for an electron <em>i </em>in the 1<em>s </em>state is given in terms of a dimensionless variable (the wave function is not properly normalized)

<strong>r</strong><em>i </em>= <em>x</em><em>i</em><strong>e</strong><em>x </em>+ <em>y</em><em>i</em><strong>e</strong><em>y </em>+ <em>z</em><em>i</em><strong>e</strong><em>z,</em>

as

<em>ψ</em>1<em>s</em>(<strong>r</strong><em>i</em>) = <em>e</em><em>−αr</em><em>i,</em>

where <em>α </em>is a parameter and

q <em>r</em><em>i </em>=            <em>x</em>2<em>i </em>+ <em>y</em><em>i</em>2 + <em>z</em><em>i</em>2<em>.</em>

We will fix <em>α </em>= 2, which should correspond to the charge of the helium atom <em>Z </em>= 2.

The ansatz for the wave function for two electrons is then given by the product of two so-called 1<em>s </em>wave functions as

Ψ(<strong>r</strong>1<em>,</em><strong>r</strong>2) = <em>e</em><em>−α</em>(<em>r</em>1+<em>r</em>2)<em>.</em>

Note that it is not possible to find a closed-form or analytical solution to Schrödinger’s equation for two interacting electrons in the helium atom.

<em> </em>c 1999-2019, “Computational Physics I

FYS3150/FYS4150″:”http://www.uio.no/studier/emner/matnat/fys/FYS3150/indexeng.html”. Released under CC Attribution-NonCommercial 4.0

license

The integral we need to solve is the quantum mechanical expectation value of the correlation energy between two electrons which repel each other via the classical Coulomb interaction, namely

<em>h </em>1 <em>i </em>= Z <em>∞ d</em><strong>r</strong><sub>1</sub><em>d</em><strong>r</strong><sub>2</sub><em>e</em><em>−</em>2<em>α</em>(<em>r</em>1+<em>r</em>2) 1 <em>. </em>(1) <em>|</em><strong>r</strong>1 <em>− </em><strong>r</strong>2<em>| </em><em>−∞ |</em><strong>r</strong>1 <em>− </em><strong>r</strong>2<em>|</em>

Note that our wave function is not normalized. There is a normalization factor missing, but for this project we don’t need to worry about that.

This integral can be solved in closed form and the answer is 5<em>π</em><sup>2</sup><em>/</em>16<sup>2</sup>. Can you derive this value?

<strong>For this project you can hand in collaborative reports and programs.</strong>

<strong>Project 3a): Gauss-Legendre Quadrature. </strong>Use Gauss-Legendre quadrature and compute the integral by integrating for each variable <em>x</em><sub>1</sub><em>,y</em><sub>1</sub><em>,z</em><sub>1</sub><em>,x</em><sub>2</sub><em>,y</em><sub>2</sub><em>,z</em><sub>2 </sub>from <em>−∞ </em>to <em>∞</em>. How many mesh points do you need before the results converges at the level of the third leading digit? Hint: the single-particle wave function <em>e<sup>−αr</sup></em><em><sup>i </sup></em>is more or less zero at <em>r<sub>i </sub>≈ λ </em>(find the appropriate limit). You can therefore replace the integration limits <em>−∞ </em>and <em>∞ </em>with <em>−λ </em>and +<em>λ</em>, respectively. You need to check that this approximation is satisfactory, that is, make a plot of the function and check if the abovementioned limits are appropriate. You need also to account for the potential problems which may arise when <em>|</em><strong>r</strong><sub>1 </sub><em>− </em><strong>r</strong><sub>2</sub><em>| </em>= 0.

<strong>Project 3b): Improved Gauss-Quadrature. </strong>The Legendre polynomials are defined for <em>x ∈ </em>[<em>−</em>1<em>,</em>1]. The previous exercise gave a very unsatisfactory ad hoc procedure. We wish to improve our results. It can therefore be useful to change to another coordinate frame and employ the Laguerre polynomials. The Laguerre polynomials are defined for <em>x ∈ </em>[0<em>,∞</em>) and if we change to spherical coordinates

<em>d</em><strong>r</strong><sub>1</sub><em>d</em><strong>r</strong><sub>2 </sub>= <em>r</em><sub>1</sub><sup>2</sup><em>dr</em><sub>1</sub><em>r</em><sub>2</sub><sup>2</sup><em>dr</em><sub>2</sub><em>dcos</em>(<em>θ</em><sub>1</sub>)<em>dcos</em>(<em>θ</em><sub>2</sub>)<em>dφ</em><sub>1</sub><em>dφ</em><sub>2</sub><em>,</em>

with

1                               1

=

<em>r</em>12 <sup>p</sup><em>r</em><sub>1</sub><sup>2 </sup>+ <em>r</em><sub>2</sub><sup>2 </sup><em>− </em>2<em>r</em><sub>1</sub><em>r</em><sub>2</sub><em>cos</em>(<em>β</em>)

and <em>cos</em>(<em>β</em>) = <em>cos</em>(<em>θ</em><sub>1</sub>)<em>cos</em>(<em>θ</em><sub>2</sub>) + <em>sin</em>(<em>θ</em><sub>1</sub>)<em>sin</em>(<em>θ</em><sub>2</sub>)<em>cos</em>(<em>φ</em><sub>1 </sub><em>− φ</em><sub>2</sub>))

we can rewrite the above integral with different integration limits. Find these

limits and replace the Gauss-Legendre approach in a) with Laguerre polynomials. The function gauss-laguerre.cpp in the <a href="https://github.com/CompPhysics/ComputationalPhysics/tree/master/doc/Projects/2019/Project3/CodeExamples">CodeExamples</a> folded can be used. Do your results improve? Compare with the results from a).

<strong>Important notice for c++ programmers</strong>: the function which computes the Gauss-Laguerre integration points and weights returns arrays which start at 1 and end <em>n </em>instead of the default values 0 and <em>n − </em>1. You need to declare an array of length <em>n </em>+ 1.

<strong>Project 3c): Monte Carlo Integration. </strong>Compute the same integral but now with brute force Monte Carlo and compare your results with those from the previous points. Discuss the differences. With bruce force we mean that you should use the uniform distribution.

<strong>Project 3d): Improved Monte Carlo Integration. </strong>Improve your brute force Monte Carlo calculation by using importance sampling. Hint: use the exponential distribution and transform to spherical coordinates. Does the variance decrease? Does the CPU time used compared with the brute force Monte Carlo decrease in order to achieve the same accuracy? Comment your results and make a list over the time each method uses. Compare the results also.

<strong>Project 3e): Improved Monte Carlo Integration and Parallization.</strong>

Finally, for the last exercise you should parallelize your code using openMP or MPI. Time your program with various compiler flags and comment these results as well. In particular, we want to see whether you achieve an optimal speed-up or not.

<h1>Introduction to numerical projects</h1>

Here follows a brief recipe and recommendation on how to write a report for each project.

<ul>

 <li>Give a short description of the nature of the problem and the eventual numerical methods you have used.</li>

 <li>Describe the algorithm you have used and/or developed. Here you may find it convenient to use pseudocoding. In many cases you can describe the algorithm in the program itself.</li>

 <li>Include the source code of your program. Comment your program properly.</li>

 <li>If possible, try to find analytic solutions, or known limits in order to test your program when developing the code.</li>

 <li>Include your results either in figure form or in a table. Remember to label your results. All tables and figures should have relevant captions and labels on the axes.</li>

 <li>Try to evaluate the reliabilty and numerical stability/precision of your results. If possible, include a qualitative and/or quantitative discussion of the numerical stability, eventual loss of precision etc.</li>

 <li>Try to give an interpretation of you results in your answers to the problems.</li>

 <li>Critique: if possible include your comments and reflections about the exercise, whether you felt you learnt something, ideas for improvements and other thoughts you’ve made when solving the exercise. We wish to keep this course at the interactive level and your comments can help us improve it.</li>

 <li>Try to establish a practice where you log your work at the computerlab. You may find such a logbook very handy at later stages in your work, especially when you don’t properly remember what a previous test version of your program did. Here you could also record the time spent on solving the exercise, various algorithms you may have tested or other topics which you feel worthy of mentioning.</li>

</ul>

<h1>Format for electronic delivery of report and programs</h1>

The preferred format for the report is a PDF file. You can also use DOC or postscript formats or as an ipython notebook file. As programming language we prefer that you choose between C/C++, Fortran2008 or Python. The following prescription should be followed when preparing the report:

<ul>

 <li>Use Devilry to hand in your projects, log in at <a href="http://devilry.ifi.uio.no/">http://devilry.ifi. </a><a href="http://devilry.ifi.uio.no/">no</a> with your normal UiO username and password and choose either</li>

</ul>

’fys3150’ or ’fys4150’. There you can load up the files within the deadline.

<ul>

 <li>Upload <strong>only </strong>the report file! For the source code file(s) you have developed please provide us with your link to your github domain. The report file should include all of your discussions and a list of the codes you have developed. Do not include library files which are available at the course homepage, unless you have made specific changes to them.</li>

 <li>In your git repository, please include a folder which contains selected results. These can be in the form of output from your code for a selected set of runs and input parametxers.</li>

 <li>In this and all later projects, you should include tests (for example unit tests) of your code(s).</li>

 <li>Comments from us on your projects, approval or not, corrections to be made etc can be found under your Devilry domain and are only visible to you and the teachers of the course.</li>

</ul>

Finally, we encourage you to work two and two together. Optimal working groups consist of 2-3 students. You can then hand in a common report.