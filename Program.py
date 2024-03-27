# Project 2
# Francis Williamson, Isaac Toapanta & Lindsay Wakley
#-----------------------------------------------------------------------------------------------
#Sphere Volume based on the radius
print ("Following code shows how we can calculate volume based on the radius")
print (" ")
import math

def calculate_sphere_volume(radius):
    """
    Calculates the volume of a sphere based on its radius.

    :param radius: The radius of the sphere.
    :return: The volume of the sphere.
    """
    volume = 4/3 * math.pi * radius**3
    return volume
# Demonstrate the function with radii 2, 4, 5, and 10
radii = [2, 4, 5, 10]
for r in radii:
    volume = calculate_sphere_volume(r)
    print(f"The volume of a sphere with radius {r} is {volume:.2f}")
#-----------------------------------------------------------------------------------------------
# Volume of a piramid
print (" ")
print ("Now, this shows how to calculate the Volume of a Piramid")
print (" ")

def calculate_pyramid_volume(base_length, height):
    """
    Calculates the volume of a pyramid with a square base.

    :param base_length: The length of one side of the base.
    :param height: The height of the pyramid.
    :return: The volume of the pyramid.
    """
    volume = (base_length**2 * height) / 3
    return volume

# Demonstrate the function with the following pairs
pairs = [(4, 4), (3, 5), (5, 2)]
for s, h in pairs:
    volume = calculate_pyramid_volume(s, h)
    print(f"The volume of a pyramid with base length {s} and height {h} is {volume:.2f}")

#-----------------------------------------------------------------------------------------------
# Slope of an Equation
print (" ")
print ("This is showing how to calculate the slope of an equation")
print (" ")

def calculate_slope(f, x0, h=1e-9):
    """
    Calculates the slope of a function f(x) at a point x0 using the limit definition of a derivative.

    :param f: The function as a callable.
    :param x0: The point at which to calculate the slope.
    :param h: The small value used to approximate the derivative. Defaults to 1e-9.
    :return: The slope of the function at the point x0.
    """
    return (f(x0 + h) - f(x0)) / h
# Example usage:
def f(x):
    return 3 * x**2 - 6 * x + 1
x0 = 1
slope = calculate_slope(f, x0)
print(f"The slope of f(x) = 3x^2-6x+1 at x = {x0} is {slope:.9f}.")
#-----------------------------------------------------------------------------------------------
# Quadratic equation Program

print (" ")
print ("This is the program for the Quadratic Equation")
print (" ")

import cmath

def solve_quadratic_equation(a, b, c):
    """
    Solves the quadratic equation ax^2 + bx + c = 0.

    :param a: The coefficient of x^2.
    :param b: The coefficient of x.
    :param c: The constant term.
    :return: A list of tuples containing the roots of the equation.
    """
    discriminant = b**2 - 4*a*c
    if discriminant > 0:
        root1 = (-b + cmath.sqrt(discriminant)) / (2*a)
        root2 = (-b - cmath.sqrt(discriminant)) / (2*a)
        return [(root1.real, root2.real)]
    elif discriminant == 0:
        root = -b / (2*a)
        return [(root,)]
    else:
        real = (-b) / (2*a)
        imaginary = cmath.sqrt(discriminant) / (2*a)
        return [(real, imaginary), (-real, -imaginary)]

# Example usage:
coefficients = (3, -4, 1)
roots = solve_quadratic_equation(*coefficients)
print(f"The roots of the equation {coefficients[0]}x^2 + {coefficients[1]}x + {coefficients[2]} = 0 are {roots}.")
#-----------------------------------------------------------------------------------------------
#Function that applies Newtons Method

print (" ")
print ("Newton’s method function")
print (" ")

def f(x):
    return 3*x**2 - 4*x + 1

def f_prime(x):
    return 6*x - 4

def newtons_method(f, f_prime, x0, tol=1e-6, max_iter=100):
    """
    Newton's method for finding a zero of a function.
    
    :param f: The function.
    :param f_prime: The derivative of the function.
    :param x0: Initial guess.
    :param tol: Tolerance for stopping criteria.
    :param max_iter: Maximum number of iterations.
    :return: Approximation of the zero.
    """
    x = x0
    for i in range(max_iter):
        x_new = x - f(x) / f_prime(x)
        if abs(x_new - x) < tol:
            return x_new
        x = x_new
    return None

# Finding zeros of f(x) = 3x^2 - 4x + 1
zero1 = newtons_method(f, f_prime, 0.5)
zero2 = newtons_method(f, f_prime, 1.5)

# Using the quadratic formula for comparison
import cmath
a = 3
b = -4
c = 1
quadratic_zeros = [(-b + cmath.sqrt(b**2 - 4*a*c)) / (2*a), (-b - cmath.sqrt(b**2 - 4*a*c)) / (2*a)]

print("Newton's method zeros:", zero1, zero2)
print("Quadratic formula zeros:", quadratic_zeros)

# Evaluating square root and cubed root of 15 using Newton's method
sqrt_15 = newtons_method(lambda x: x**2 - 15, lambda x: 2*x, 4)
cbrt_15 = newtons_method(lambda x: x**3 - 15, lambda x: 3*x**2, 2)

print("Square root of 15:", sqrt_15)
print("Cubed root of 15:", cbrt_15)

#-----------------------------------------------------------------------------------------------
# Riemann summation
print (" ")
print ("Riemann summation to stimate the areas")
print (" ")
def riemann_sum(f, a, b, n, method ='left'):
    """
    Estimates the area under a curve using a Riemann summation.

    :param f: The function as a callable.
    :param a: The left endpoint of the interval.
    :param b: The right endpoint of the interval.
    :param n: The number of rectangles.
    :param method: The method to use for the Riemann summation. Can be either 'left', 'right', or 'middle'.
    :return: The estimated area under the curve.
    """
    if method == 'left':
        x_values = [a + i * (b - a) / n for i in range(n)]
    elif method == 'right':
        x_values = [a + (i + 1) * (b - a) / n for i in range(n)]
    elif method == 'middle':
        x_values = [a + (i + 0.5) * (b - a) / n for i in range(n)]
    else:
        raise ValueError("Invalid method for the Riemann summation.")
    y_values = [f(x) for x in x_values]
    rectangles = [y_values[i] * (b - a) / n for i in range(n)]
    return sum(rectangles)

print ("This are the requirements solved from the left")
# Example usage by left:
def f(x):
    return x**3 + 3*x**2 + 4*x - 1

a = 1
b = 4
n = 100
area = riemann_sum(f, a, b, n)
method = 'left'
print(f"The estimated area under the curve y = x^3 + 3x^2 + 4x - 1 from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by left:
import math

def f(x):
    return math.log(x)

a = 5
b = 10
n = 100
area = riemann_sum(f, a, b, n)
method = 'left'
print(f"The estimated area under the curve y = ln(x) from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by left:
def f(x):
    return math.sin(x)

a = 0
b = math.pi
n = 100
area = riemann_sum(f, a, b, n)
method = 'left'
print(f"The estimated area under the curve y = sin(x) from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by left:
def f(x):
    return math.exp(x)

a = 0
b = 2
n = 100
area = riemann_sum(f, a, b, n)
method = 'left'
print(f"The estimated area under the curve y = e^x from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

print ("This are the requirements solved from the right")
# Example usage by right:
def f(x):
    return x**3 + 3*x**2 + 4*x - 1

a = 1
b = 4
n = 100
area = riemann_sum(f, a, b, n)
method = 'right'
print(f"The estimated area under the curve y = x^3 + 3x^2 + 4x - 1 from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by right:
import math

def f(x):
    return math.log(x)

a = 5
b = 10
n = 100
area = riemann_sum(f, a, b, n)
method = 'right'
print(f"The estimated area under the curve y = ln(x) from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by right:
def f(x):
    return math.sin(x)

a = 0
b = math.pi
n = 100
area = riemann_sum(f, a, b, n)
method = 'right'
print(f"The estimated area under the curve y = sin(x) from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by right:
def f(x):
    return math.exp(x)

a = 0
b = 2
n = 100
area = riemann_sum(f, a, b, n)
method = 'right'
print(f"The estimated area under the curve y = e^x from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

print ("This are the requirements solved from the middle")
# Example usage by middle:
def f(x):
    return x**3 + 3*x**2 + 4*x - 1

a = 1
b = 4
n = 100
area = riemann_sum(f, a, b, n)
method = 'middle'
print(f"The estimated area under the curve y = x^3 + 3x^2 + 4x - 1 from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by middle:
import math

def f(x):
    return math.log(x)

a = 5
b = 10
n = 100
area = riemann_sum(f, a, b, n)
method = 'middle'
print(f"The estimated area under the curve y = ln(x) from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by middle:
def f(x):
    return math.sin(x)

a = 0
b = math.pi
n = 100
area = riemann_sum(f, a, b, n)
method = 'middle'
print(f"The estimated area under the curve y = sin(x) from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")

# Example usage by middle:
def f(x):
    return math.exp(x)

a = 0
b = 2
n = 100
area = riemann_sum(f, a, b, n)
method = 'middle'
print(f"The estimated area under the curve y = e^x from x = {a} to x = {b} is {area:.2f} (using {method} Riemann summation with {n} rectangles).")