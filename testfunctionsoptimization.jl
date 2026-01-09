#! test functions for optimization
# https://en.wikipedia.org/wiki/Test_functions_for_optimization

##
# functions

# Rastrigin function
# rastrigin(0) == 0
# search domain: -5.12 <= x[i] <= 5.12
function rastrigin(x, A=10)
    n = length(x)
    A*n + sum(x[i]^2 - A*cos(2π*x[i]) for i in 1:n)
end

# Ackley function
# ackley(0) == 0
# search domain: -5 <= x[i] <= 5
function ackley(x, params=[20, 0.2, 2π, 2])
    n = length(x)
    a, b, c, d = (params...,)
    -a * exp(-b * sqrt(1/d * sum(x[i]^2 for i in 1:n))) - exp(1/d * sum(cos(c * x[i]) for i in 1:n)) + a + exp(1)
end

# Sphere function
# sphere(0) == 0
# search domain: -inf <= x[i] <= inf
function sphere(x)
    n = length(x)
    sum(x[i]^2 for i in 1:n)
end

# Rosenbrock function
# rosenbrock(1) == 0
# search domain: -inf <= x[i] <= inf
function rosenbrock_2(xx, params)
    a, b = (params...,)
    x = xx[1]
    y = xx[2]
    return (a - x)^2 + b*(y - x^2)^2 
end

function rosenbrock_n2(x)
    n = length(x)
    isodd(n) && return throw(error("length of x must be even"))
    n2 = div(n, 2)
    sum(100(x[2i-1]^2 - x[2i])^2 + (x[2i-1] - 1)^2 for i in 1:n2)
end

function rosenbrock_n(x)
    n1 = length(x) - 1
    sum(100(x[i+1] - x[i]^2)^2 + (1 - x[i])^2 for i in 1:n1)
end

# Beale function
# beale([3, 0.5]) == 0
# search domain: -4.5 <= x[i] <= 4.5
function beale(xx)
    x = xx[1]
    y = xx[2]
    (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2
end

# Goldstein-Price function
# goldsteinprice([0, -1]) == 3
# search domain: -2 <= x[i] <= 2
function goldsteinprice(xx)
    x = xx[1]
    y = xx[2]
    (1 + (x + y + 1)^2 * (19 - 14x + 3x^2 - 14y + 6x*y + 3y^2))*
    (30 + (2x - 3y)^2 * (18 - 32x + 12x^2 + 48y - 36x*y + 27y^2))
end

# Booth function
# booth([1, 3]) == 0
# search domain: -10 <= x[i] <= 10
function booth(xx)
    x = xx[1]
    y = xx[2]
    (x + 2y - 7)^2 + (2x + y - 5)^2
end

# Burkin function no. 6
# burkin6([-10, 1]) == 0
# search domain: -15 <= x <= -5, -3 <= y <= 3
function burkin6(xx)
    x = xx[1]
    y = xx[2]
    100sqrt(abs(y - 0.01x^2)) + 0.01abs(x + 10)
end

# Matyas function
# matyas(0) == 0
# search domain: -10 <= x[i] <= 10
function matyas(xx)
    x = xx[1]
    y = xx[2]
    0.26(x^2 + y^2) - 0.48x*y
end

# Lévi function no. 13
# levi13(1) == 0
# search domain: -10 <= x[i] <= 10
function levi13(xx)
    x = xx[1]
    y = xx[2]
    sin(3π*x)^2 + (x - 1)^2 * (1 + sin(3π*y)^2) + (y - 1)^2 * (1 + sin(2π*y)^2)
end

# Griewank function
# griewank(0) == 0
# search domain: -inf <= x[i] <= inf
function griewank(x)
    n = length(x)
    p(x, i) = cos(x[i]/sqrt(i))
    1 + 1/4000 * sum(x[i]^2 for i in 1:n) - prod(p(x[i], i) for i in 1:n)
end

function griewank_ns(x)
    n = length(x)
    p(x, i) = abs(cos(x[i]/(2sqrt(i)))) - abs(sin(x[i]/(2sqrt(i))))
    1 + 1/4000 * sum(x[i]^2 for i in 1:n) - prod(p(x[i], i) for i in 1:n)
end

# Himmelblau function
# himmelblau([3, 2]) == himmelblau([-2.805118, 3.131312]) == himmelblau([-3.779310, -3.283186]) == himmelblau([3.584428, -1.848126]) == 0
# search domain: -5 <= x[i] <= 5
function himmelblau(xx)
    x = xx[1]
    y = xx[2]
    (x^2 + y - 11)^2 + (x + y^2 - 7)^2
end

# Three-hump camel function
# camel(0) == 0
# search domain: -5 <= x[i] <= 5
function camel(xx)
    x = xx[1]
    y = xx[2]
    2x^2 - 1.05x^4 + x^6/6 + x*y + y^2
end

# Easom function
# easom([π, π]) == -1
# search domain: -100 <= x[i] <= 100
function easom(xx)
    x = xx[1]
    y = xx[2]
    -cos(x) * cos(y) * exp(-((x - π)^2 + (y - π)^2))
end

# Cross-in-tray function
# crossintray([1.34941, -1.34941]) == crossintray([1.34941, 1.34941]) == crossintray([-1.34941, 1.34941]) == crossintray([-1.34941, -1.34941]) == -2.06261
# search domain: -10 <= x[i] <= 10
function crossintray(xx)
    x = xx[1]
    y = xx[2]
    -0.0001(abs(sin(x) * sin(y) * exp(abs(100 - sqrt(x^2 + y^2)/π))) + 1)^0.1
end

# Eggholder function
# eggholder([512, 404.2319]) == -959.6407
# search domain: -512 <= x[i] <= 512
function eggholder(xx)
    x = xx[1]
    y = xx[2]
    -(y + 47) * sin(sqrt(abs(x/2 + (y + 47)))) - x * sin(sqrt(abs(x - (y + 47))))
end

# Hölder table function
# holder([8.05502, 9.66459]) == holder([-8.05502, 9.66459]) == holder([8.05502, -9.66459]) == holder([-8.05502, -9.66459]) == -19.2085
# search domain: -10 <= x[i] <= 10
function holder(xx)
    x = xx[1]
    y = xx[2]
    -abs(sin(x) * cos(y) * exp(abs(1 - sqrt(x^2 + y^2)/π)))
end

# McCormick function
# mccormick([-0.54719, -1.54719]) == -1.9133
# search domain: -1.5 <= x <= 4, -3 <= y <= 4
function mccormick(xx)
    x = xx[1]
    y = xx[2]
    sin(x + y) + (x - y)^2 - 1.5x + 2.5y + 1
end

# Schaffer function no. 2
# schaffer2(0) == 0
# search domain: -100 <= x[i] <= 100
function schaffer2(xx)
    x = xx[1]
    y = xx[2]
    0.5 + (sin(x^2 - y^2)^2 - 0.5)/(1 + 0.001(x^2 + y^2))^2
end

# Schaffer function no. 4
# schaffer4([0, 1.25313]) == schaffer4([0, -1.25313]) == schaffer4([1.25313, 0]) == schaffer4([-1.25313, 0]) == 0.292579
# search domain: -100 <= x[i] <= 100
function schaffer4(xx)
    x = xx[1]
    y = xx[2]
    0.5 + (cos(sin(abs(x^2 - y^2)))^2 - 0.5)/(1 + 0.001(x^2 + y^2))^2
end

# Styblinski-Tang function
# -39.16617n < styblinskitang(-2.903534) < -39.16616n
# search domain: -5 <= x[i] <= 5
function styblinskitang(x)
    n = length(x)
    sum(x[i]^4 - 16x[i]^2 + 5x[i] in 1:n)/2
end

# Shekel function
#
# search domain: -inf <= x[i] <= inf
function shekel(x, c=repeat([1.2], 30), a=repeat([3.5], length(x), length(c)))
    n = length(x)
    m = length(c)
    sum(inv(c[i] + sum((x[j] - a[j, i]))^2 for j in 1:n) for i in 1:m)
end
