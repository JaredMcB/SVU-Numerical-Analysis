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

end # module