module NumCalc

"""
Newton-Cotes by undetermined coefficients
By Jared McBride (3-13-2024, Buena Vista, Virginia)

## Example
````
f = x -> sin(x)^2 - 2x*sin(x) + 1

b = 1.3
a = 0.75

newton_cotes(f,a,b, n = 4).S
````
"""
function newton_cotes(f, a, b; n = 5, closed = true)
    # Choose h and x-grid according to n and open/closed
    if closed
        h = (b-a)/n
        x = a : h : b
    else
        h = (b-a)/(n+2)
        x = a+h : h : b-h/2
    end

    
    # Vector of Data
    A = [(b^(i+1) - a^(i+1))/(i+1) for i = 0:n]
    
    # Design matrix B_ij = [x^(i-1)_(j-1)]
    B = cat(dims = 2, [[x^i for i = 0:n] for x in x]...) 

    as = B \ A
    return (I = as'f.(x), # The accual quadrature
            S = as/h)     # Quadrature scheme
end

"""
Composite Newton-Cotes
by Jared McBride (3-25-2024, Buena Vista, VA)

Uses the function `newton_cotes` to form a composite quadrature rule. 

inputs:
    f - integrand,
    a, b - left and right end points of he whole interval of integration,
    num_blks - the number of subintervals,
    n - degree of quadrature, as in we take n+1 points from each subinterval resulting in `n*num_blks + 1` total points,

outputs:

"""
function composite_newton_cotes(f, 
    a,                     # Left limit of integration
    b; 
    num_blks   = 25,       # Total number of points 
    n          = 3,        # degree of quadrature
    tot_points = num_blks * n - (num_blks -1), # Total number of points used
    x          = a : (b - a) / (tot_points - 1) : b,
    closed     = true
)
h = x[2] - x[1]
S = zeros(tot_points)
for m = 1:num_blks
    S[(n-1)*(m-1)+1 : (n-1)*m+1] += 
        nc.newton_cotes(f, x[(n-1)*(m-1)+1], x[(n-1)*m+1]; n = n-1, closed).S
end

x = a : (b - a) / (tot_points - 1) : b

return (I = h*S'f.(x), # The accual quadrature
        S = S,         # Quadrature scheme
        x = x)
end

end # module