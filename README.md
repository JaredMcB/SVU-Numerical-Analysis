# Southern Virginia University Numerical Analysis
## Class Project

This repository contains a grab bag of resources for the learning and practice of numerical analysis at the undergraduate level. It is contributed to by the students MAT 410: _Introduction to Numerical Analysis_ at Southern Virginia University, Buena Vista, Virginia, spring 2024. The students will be adding to the course throughout the semester.  

Soon we will be starting Chapter 2 of Burden, which covers solving equations of a single variable. In preparation for this, I have started a Julia module called EquSolvers.jl. All the code projects for chapter two will be held in this module. 

## Root-finding methods: EquSolvers.jl

### Inverse Quadratic Interpolation

Inverse Quadratic Interpolation improves on the secant method by fitting a sideways parabola (i.e. a parabola is $y$) to three successive points of the iteration. The recursion boils down to the following. 

$$ 
\begin{align}
p_{n+1} = & \frac{f_{n-1}f_{n}}{(f_{n-2} - f_{n-1})(f_{n-2}-f_{n})}p_{n-2} + \frac{f_{n-2}f_{n}}{(f_{n-1} - f_{n-2})(f_{n-1}-f_{n})}p_{n-1} + \frac{f_{n-2}f_{n-1}}{(f_{n} - f_{n-2})(f_{n}-f_{n-1})}p_{n}
 \end{align}
$$

### Newton's Method
One of the benefits of using Newton's Method is the upside of it's intuitiveness. Much as the same way we look at functions and consider them conceptually as plots on some n-planes, 2-space, 3-space, etc. we can consider Newton's Method as being able to do the same thing when our goal is the synthesis of an algorithm. Newton's Method takes our function, some f, and we can either manually, as was done in my case, or we can have some program find the derivative of our function f, to give us an f'. Before we can start and numerical solvers, we need to define a maximum number of iterations of Newton's Method so that we do not encounter and runtime errors, should we input or code something improperly. Then once we have our max. iterations, we can start on actually using Newton's Method to find our zero's which for our case are the numercial solutions we are finding for. We already have our function f, and through some method we found f'. No we need to start Newton's Method with some initial guess, x_0, which we hope to be as close to our actual zero or root or solution, x. We then find the value of f at x_0 which gives us some value and we find the derivative at our x_0. We then subtract the quotient of \frac{f(x)}{f'(x)} which we then define as x_0 again and repeat for as long as the absolute value of that value is greater than our tolerance we defined earlier.

### Fixed Point Method

The fixed point method finds values of our function where the input doesn't change the value of the function entire. This method works by iterating outputs and looping them into out input as long as specific conditions are met, such as ensuring that the derivative at the eventual fixed point is less than 1. The fixed point method determines where g(x) is equal to x by taking an initial value, our p_0, and plugging it into our function g. Whatever value is outputted is then set as our new p_0 and a new value is determined. While the fixed point method is not an explicit root finding method, it is a powerful tool for deriving them. 
