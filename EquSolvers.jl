module EquSolver

function MC_log(x; N = 1000)
    count = 1

    for i = 1: N
        Y = rand()
        X = (x-1)*rand() + 1

        if Y < 1/X
            count = count +1
        end
    end
    I = count/N * (x-1)
end


function hal(Dave)
    println("I can't do that, $Dave.")
end

## Bisection method

"""
# Bisection method
by Caden Eskelsen (Mar 7, 2024)

An implementation of the bisection method [Numerical Analysis 9e, Burden, Sec 2.1]. Using the intermediate value theorem, the bisection method is a binary search method. For more information on bisection see [Wikipedia](https://en.wikipedia.org/wiki/Bisection_method)

### Inputs:
f - the function whose roots we are looking for 

a, b - the interval in which we are looking for a root

maxit - maximum number of iterations (if not specified the function will bisect 50 times)

tol - error tolorance (set at 1e-10)

### Outputs:
p - the root found within the interval provided

j - number of iterations 

ps - a vector with the midpoints after each iteration of bisection

## Example:
```
f = x -> x*cos(x)
a = 1
b = 2

bisection(f,a,b)

```

"""
function bisection(f,a,b; maxit = 50, tol = 1e-10)

    ps = zeros(maxit)
    j = 1               # Initialize counter
    fa = f(a)           # Function value at left end
    p = 0      # Form is prefered to (a+b)/2

    while j < maxit  # iterates the loop while less than the specified maximum iterations
        p = a + (b-a)/2  # finds midpoint of inital guesses a and b
        fp = f(p)  
        if fp == 0 || (b-a)/2 < tol    # check p is the zero or error is with in tolerance
            break
        end

        if fa * fp > 0       # Check if fa and fp are different signs
            a = p            # replace a with p for bisection in next iteration
            fa = fp
        else
            b = p            # replace b with p for bisection in next iteration
        end
        ps[j] = p            # adds to vector of midpoints
        j = j + 1           
    end

    return p, j, ps[1:j-1]
end


## Fixed point iterations

"""
fixed_point by Aaron Miller

The fixed point method determines where g(x) is equal to x for some value x. This is done by taking an initial value x (written as p0 
    in the method code) and plugging it into our user-defined function written as f in the function. Whatever value is outputted is 
    then set as our new p0 and a new value is determined. This new value is p in the while statement. This continues until the fixed 
    point is found, where (x,g(x)) where x = g(x).

Notes on choosing a function  f and an initial point p0: when choosing a function, it is important to note that the derivative of the
    function at the initial point f(p0) and the derivative at the eventual fixed point must be less than 1. In addition, it is important
    to note that, for some function g(x), the values g(a) and g(b) must be within a and b, where a and b are values of x. 
    
"""

function fixed_point(f, p0; tol=1e-10, maxit=50)  # Our inputs for the function are the function itself and our initial point.
    i=1     # We start i, our iteration count, at 1. It will end the 'while' statement once it gets reaches (maxit)+1
    p=0     # Defines p used in the while statement in the global
    while i<=maxit  # The function will continue until this condition is breached
        p=f(p0)     # Defining p (formerly 0) as f(p0)
        if abs(p-p0)<=tol # If they are essentially the same, the process was successful and we move to return
            break
        end
        i=i+1 # Change the iteration count 
        p0=p  # We get our new p0, whatever the output of f(p0) was above
    end
    if i == maxit+1
        println("Method Failed")
    end
    return p, i # Return the fixed point (as determined by the if abs(p-p0)<= tol statement) and i, the number of iterations
end

## Newtons method
using Plots
f(x)= x^3+2x+17     # first function aka f(x)
f_prime(x)= 3x^2+2  # second function aka f'(x)

"""
Newton's Method
by Patrick Whitehouse (3-08-2024, Buena Vista, VA)

This is an application of Newton's Method as illustrated by Numerical Analysis 9e, Burden

Inputs:
f - the originial function
f' - the deriative of the original function
x0 - our initial x value
tolerance - the range of acceptable differences from our f and f'
max_iterations - the maximum number of iterations that the process will undergo, this is implemented to prevent runtime errors

Outputs:
Root found: The actual root that we found using Newton's Method
Residual: The difference between the "true" root and the calculated root that falls within an acceptable range of the tolerance
Number of iterations: this is the number of iterations that the process underwent before falling within an acceptable tolerance range
"""
function newtons_method(f, f_prime, x0, tolerance, max_iterations)  #= setting up and instatiating the guidelines and necessary bounds for NM=#
    x=x0
    iterations = 0

    while abs(f(x)) > tolerance && iterations < max_iterations
        x-= f(x) / f_prime(x)
        iterations += 1
    end

    return x, f(x), iterations
end

x0 = 1.0
tolerance = 1e-10       # increased tolerance thorugh some trial and error for optimizing residuals
max_iterations = 100    #added runtime error safeguard#


root, residual, iterations = newtons_method(f, f_prime, x0, tolerance, max_iterations)


println("Root found: $root")
println("Residual: $residual")
println("Number of iterations: $iterations")


x_values = range(-2, stop=2, length=100)
y_values = f.(x_values)

plot(x_values,y_values,label="Function")
scatter!([root], [0], label="Root", color="red")
title!("Newton's Method Convergence of function, x^3+2x+17")
xlabel!("x")
ylabel!("y")
##############################################################
## Inverse Quadratic interpolation ###########################
##############################################################
""" 
# Inverse quadratic interpolation
By Jared McBride (2-20-2025, Buena Vista, VA)

This is an extension of the secant method method. But rather than using the secant line from two points it uses an arc from a sideways parabola, interpolated by three points. For more information see [Wikipedia](https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation)

Inputs:
f - the function whose root we seek
p0, p1, p2 - three initial points.
N - max number of iterations
tol - error tolerance

outputs:
p - the approximated roots
n - the number of iterations it took
Ps - a history of approximate roots

## Example:
```
f = x -> cos(x) - x
p0 = .5
p1 = .6

p = QuadNewton(f, p0, p1).p
```
"""
function QuadNewton(f, p0, p1; 
    # Use secant to get a third point
    p2  = p0 - (p1 - p0) / (f(p1) - f(p0)) * f(p0), 
    N   = 30,          # Default max number of iterations
    tol = 1e-10        # Defauls tolerance
    )

    
    P = [p0; p1; p2]    # last three approximations
    F = f.(P)           # f evaluated at P
    c = zeros(3)        # initializing vector c
    p = 0
    n = 1
    P_all = zeros(N)
    
    while n < N && abs(P[3] - P[2]) > tol
        # Compute coefficients for the quadratic
        c[1] = (F[2] / (F[1] - F[2])) * (F[3] / (F[1] - F[3]))
        c[2] = (F[1] / (F[2] - F[1])) * (F[3] / (F[2] - F[3]))
        c[3] = (F[1] / (F[3] - F[1])) * (F[2] / (F[3] - F[2]))

        # Get new point
        p = c'P 
        P_all[n] = p       # Save p

        # Roll in the new point
        P = [P[2:3]; p]
        F = [F[2:3]; f(p)]
        
        # increment counter
        n += 1
    end

    return (p = p,
            n = n,
            Ps = P_all)
end





        







end # module
