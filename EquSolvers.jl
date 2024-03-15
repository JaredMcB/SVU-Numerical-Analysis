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

function fixed_point(f, p; tol = 1e-6, maxit = 50)
    x_prev = p
    for i in 1:maxit
        x_next = f(x_prev)
        if abs(x_next - x_prev) < tol
            println("Root found at: ", x_next, " with tolerance: ", tol, " in ", i, " iterations.")
            return x_next
        end
        x_prev = x_next
    end
    return x_prev
end

## Newtons method

function newton_method(f, f_p, p0; tol=1e-6, maxit=100)
    x = p0
    i = 0
    while abs(f(x)) > tol && i < maxit
        x -= f(x) / f_p(x)
        i += 1
    end
    return x, i
end

## Secant method

function secant_method(f, x0, x1; tol=1e-6, maxit=100)
    i = 0
    while abs(f(x1)) > tol && i < maxit
        x_next = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        x0 = x1
        x1 = x_next
        i += 1
    end
    return x1, i
end

## False Position

function false_position(f, a, b; tol=1e-6, maxit=100)
    if f(a) * f(b) > 0
        error("The function must have opposite signs at the endpoints.")
    end
    
    i = 0
    while abs(f(b) - f(a)) > tol && i < maxit
        c = (a * f(b) - b * f(a)) / (f(b) - f(a))
        if f(c) * f(a) < 0
            b = c
        else
            a = c
        end
        i += 1
    end
    return (a + b) / 2, i
end

end # module

