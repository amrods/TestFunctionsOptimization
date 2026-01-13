# test functions for optimization
# https://en.wikipedia.org/wiki/Test_functions_for_optimization

# Rastrigin function
# rastrigin(0) == 0
# search domain: -5.12 <= x[i] <= 5.12
rastrigin <- function(x, A = 10) {
  n <- length(x)
  A * n + sum(x^2 - A * cos(2 * pi * x))
}

# Ackley function
# ackley(0) == 0
# search domain: -5 <= x[i] <= 5
ackley <- function(x, params = c(20, 0.2, 2 * pi, 2)) {
  n <- length(x)
  a <- params[1]; b <- params[2]; c_ <- params[3]; d <- params[4]
  -a * exp(-b * sqrt((1 / d) * sum(x^2))) - exp((1 / d) * sum(cos(c_ * x))) + a + exp(1)
}

# Sphere function
# sphere(0) == 0
# search domain: -inf <= x[i] <= inf
sphere <- function(x) {
  sum(x^2)
}

# Rosenbrock function
# rosenbrock(1) == 0
# search domain: -inf <= x[i] <= inf
rosenbrock_2d <- function(xx, params) {
  a <- params[1]; b <- params[2]
  x <- xx[1]; y <- xx[2]
  (a - x)^2 + b * (y - x^2)^2
}

rosenbrock_n2 <- function(x) {
  n <- length(x)
  if (n %% 2 == 1) stop("length of x must be even")
  n2 <- n %/% 2
  i  <- 1:n2
  x1 <- x[2 * i - 1]
  x2 <- x[2 * i]
  sum(100 * (x1^2 - x2)^2 + (x1 - 1)^2)
}

rosenbrock_n <- function(x) {
  n1 <- length(x) - 1
  i  <- 1:n1
  sum(100 * (x[i + 1] - x[i]^2)^2 + (1 - x[i])^2)
}

# Beale function
# beale(c(3, 0.5)) == 0
# search domain: -4.5 <= x[i] <= 4.5
beale <- function(xx) {
  x <- xx[1]; y <- xx[2]
  (1.5   - x + x * y)^2 + (2.25  - x + x * y^2)^2 + (2.625 - x + x * y^3)^2
}

# Goldstein-Price function
# goldstein_price(c(0, -1)) == 3
# search domain: -2 <= x[i] <= 2
goldstein_price <- function(xx) {
  x <- xx[1]; y <- xx[2]
  (1 + (x + y + 1)^2 * (19 - 14 * x + 3 * x^2 - 14 * y + 6 * x * y + 3 * y^2)) *
    (30 + (2 * x - 3 * y)^2 * (18 - 32 * x + 12 * x^2 + 48 * y - 36 * x * y + 27 * y^2))
}

# Booth function
# booth(c(1, 3)) == 0
# search domain: -10 <= x[i] <= 10
booth <- function(xx) {
  x <- xx[1]; y <- xx[2]
  (x + 2 * y - 7)^2 + (2 * x + y - 5)^2
}

# Bukin function no. 6
# bukin6(c(-10, 1)) == 0
# search domain: -15 <= x <= -5, -3 <= y <= 3
bukin6 <- function(xx) {
  x <- xx[1]; y <- xx[2]
  100 * sqrt(abs(y - 0.01 * x^2)) + 0.01 * abs(x + 10)
}

# Matyas function
# matyas(0) == 0
# search domain: -10 <= x[i] <= 10
matyas <- function(xx) {
  x <- xx[1]; y <- xx[2]
  0.26 * (x^2 + y^2) - 0.48 * x * y
}

# Lévi function no. 13
# levi13(1) == 0
# search domain: -10 <= x[i] <= 10
levi13 <- function(xx) {
  x <- xx[1]; y <- xx[2]
  sin(3 * pi * x)^2 +
    (x - 1)^2 * (1 + sin(3 * pi * y)^2) +
    (y - 1)^2 * (1 + sin(2 * pi * y)^2)
}

# Griewank function
# griewank(0) == 0
# search domain: -inf <= x[i] <= inf
griewank <- function(x) {
  n <- length(x)
  1 + (1 / 4000) * sum(x^2) - prod(cos(x / sqrt(1:n)))
}

griewank_ns <- function(x) {
  n <- length(x)
  p <- abs(cos(x / (2 * sqrt(1:n)))) - abs(sin(x / (2 * sqrt(1:n))))
  1 + (1 / 4000) * sum(x^2) - prod(p)
}

# Himmelblau function
# himmelblau([3, 2]) == himmelblau(c(-2.805118, 3.131312)) == himmelblau(c(-3.779310, -3.283186)) == himmelblau(c(3.584428, -1.848126)) == 0
# search domain: -5 <= x[i] <= 5
himmelblau <- function(xx) {
  x <- xx[1]; y <- xx[2]
  (x^2 + y - 11)^2 + (x + y^2 - 7)^2
}

# Three-hump camel function
# camel(0) == 0
# search domain: -5 <= x[i] <= 5
camel <- function(xx) {
  x <- xx[1]; y <- xx[2]
  2 * x^2 - 1.05 * x^4 + x^6 / 6 + x * y + y^2
}

# Easom function
# easom(c(pi, pi)) == -1
# search domain: -100 <= x[i] <= 100
easom <- function(xx) {
  x <- xx[1]; y <- xx[2]
  -cos(x) * cos(y) * exp(-((x - pi)^2 + (y - pi)^2))
}

# Cross-in-tray function
# cross_in_tray(c(1.34941, -1.34941)) == cross_in_tray(c(1.34941, 1.34941)) == cross_in_tray(c(-1.34941, 1.34941)) == cross_in_tray(c(-1.34941, -1.34941)) == -2.06261
# search domain: -10 <= x[i] <= 10
cross_in_tray <- function(xx) {
  x <- xx[1]; y <- xx[2]
  -0.0001(abs(sin(x) * sin(y) * exp(abs(100 - sqrt(x^2 + y^2)/pi))) + 1)^0.1
}

# Eggholder function
# eggholder(c(512, 404.2319)) == -959.6407
# search domain: -512 <= x[i] <= 512
eggholder <- function(xx) {
  x <- xx[1]; y <- xx[2]
  -(y + 47) * sin(sqrt(abs(x / 2 + (y + 47)))) -
    x * sin(sqrt(abs(x - (y + 47))))
}

# Hölder table function
# holder(c(8.05502, 9.66459)) == -19.2085 etc.
# search domain: -10 <= x[i] <= 10
holder <- function(xx) {
  x <- xx[1]; y <- xx[2]
  -abs(sin(x) * cos(y) * exp(abs(1 - sqrt(x^2 + y^2) / pi)))
}

# McCormick function
# mccormick(c(-0.54719, -1.54719)) == -1.9133
# search domain: -1.5 <= x <= 4, -3 <= y <= 4
mccormick <- function(xx) {
  x <- xx[1]; y <- xx[2]
  sin(x + y) + (x - y)^2 - 1.5 * x + 2.5 * y + 1
}

# Schaffer function no. 2
# schaffer2(c(0, 0)) == 0
# search domain: -100 <= x[i] <= 100
schaffer2 <- function(xx) {
  x <- xx[1]; y <- xx[2]
  0.5 + (sin(x^2 - y^2)^2 - 0.5) / (1 + 0.001 * (x^2 + y^2))^2
}

# Schaffer function no. 4
# schaffer4(c(0, 1.25313)) == 0.292579 etc.
# search domain: -100 <= x[i] <= 100
schaffer4 <- function(xx) {
  x <- xx[1]; y <- xx[2]
  0.5 + (cos(sin(abs(x^2 - y^2)))^2 - 0.5) / (1 + 0.001 * (x^2 + y^2))^2
}

# Styblinski-Tang function
# -39.16617n < styblinski_tang(-2.903534) < -39.16616n
# search domain: -5 <= x[i] <= 5
styblinski_tang <- function(x) {
  sum(x^4 - 16 * x^2 + 5 * x) / 2
}

# Shekel function
# 
# search domain: -inf <= x[i] <= inf
shekel <- function(x, c = rep(1.2, 30), a = matrix(3.5, nrow = length(x), ncol = length(c))) {
  n <- length(x)
  m <- length(c)
  s <- 0
  for (i in 1:m) {
    s <- s + 1 / (c[i] + sum((x - a[, i])^2))
  }
  s
}
