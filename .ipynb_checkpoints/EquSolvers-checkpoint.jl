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
Bisection method
by Jared McBride (Feb 7, 2024)

An implementation of the bisection method [Numerical Analysis 9e, Burden, Sec 2.1].

"""
function bisect_method(f,a,b; maxit = 50, tol = 1e-10)

    i = 1               # Initialize counter
    fa = f(a)           # Function value at left end
    p = a +(b-a)/2      # Form is prefered to (a+b)/2

    for i = 1:maxit

        fp = f(p)
        if fp == 0 || (b-a)/2 < tol     # check p is the zero or error is with in tolerance
            break
        end

        if fa * fp > 0          # Check if fa and fp are different signs
            a = p 
            fa = fp
        else
            b = p
        end
    end

    return p , i
end


























## Fixed point iterations


## Newtons method


end # module
