from math import sqrt, sin, cos, tan

e = 1e-5
A = 0.982
B = -1.298
a = 0.374
b = 1.575
c = -0.546
d = 0.476

# First


def y_k_1(x):
    return (-sin(x + A) + c)/b


def x_k_1(y):
    return d - cos(y + B)


def func_1(x, y):
    return sin(x+A) + b*y - c


def func_2(x, y):
    return x + cos(y + B) - d


def norma(x, y):
    return sqrt(x**2 + y**2)


def simple_iteration(x_k, y_k):
    i = 1
    fun1 = func_1(x_k, y_k)
    fun2 = func_2(x_k, y_k)
    print("i =", i)
    print("(x{0}, y{0}) = ({1} , {2})".format(i, x_k, y_k))
    print("(f1(x{0}, y{0}), f2(x{0}), y{0})) = ({1}, {2})".format(i, fun1, fun2))
    print("||F|| =", norma(fun1, fun2))

    while True:
        i += 1
        x_k1 = x_k_1(y_k)
        y_k1 = y_k_1(x_k)
        dx = x_k1 - x_k
        dy = y_k1 - y_k
        fun1 = func_1(x_k1, y_k1)
        fun2 = func_2(x_k1, y_k1)
        print("\ni =", i)
        print("(x{0}, y{0}) = ({1} , {2})".format(i, x_k1, y_k1))
        print("(f1(x{0}, y{0}), f2(x{0}), y{0})) = ({1}, {2})".format(i, fun1, fun2))
        print("||F|| =", norma(fun1, fun2))
        x_k = x_k1
        y_k = y_k1
        if norma(fun1, fun2) < e and norma(dx, dy) < e:
            break
    print("\n\nSolution: ({0}, {1})".format(x_k1, y_k1))

# SECOND


def det(x, y):
    return (2*b*(y**2) - 2*a*(x**2))/((cos(x*y + A))**2) - 4*b*x*y


def func_12(x, y):
    return tan(x*y + A) - x**2


def func_22(x, y):
    return a*(x**2) + b*(y**2) - 1


def xn_k_1(x, y):
    return x - (2*b*y*func_12(x, y) - x*func_22(x, y)/((cos(x*y + A))**2))/det(x, y)


def yn_k_1(x, y):
    return y - (-2*a*x*func_12(x, y) + (y/((cos(x*y + A))**2) - 2*x)*func_22(x, y))/det(x, y)


def newton(x_k, y_k):
    i = 1
    f1 = func_12(x_k, y_k)
    f2 = func_22(x_k, y_k)
    print("i =", i)
    print("(x{0}, y{0}) = ({1}, {2})".format(i, x_k, y_k))
    print("(f1(x{0} , y{0}), f2(x{0} , y{0}) = {1}, {2}".format(i, f1, f2))
    print("||F|| =", norma(f1, f2))
    while True:
        i += 1
        x_k1 = xn_k_1(x_k, y_k)
        y_k1 = yn_k_1(x_k, y_k)
        dx = x_k1 - x_k
        dy = y_k1 - y_k
        f1 = func_12(x_k1, y_k1)
        f2 = func_22(x_k1, y_k1)
        print("\ni =", i)
        print("(x{0}, y{0}) = ({1}, {2}".format(i, x_k1, y_k1))
        print("(f1(x{0} , y{0}), f2(x{0} , y{0}) = {1}, {2}".format(i, f1, f2))
        print("||F|| =", norma(f1, f2))
        x_k = x_k1
        y_k = y_k1
        if norma(f1, f2) < e and norma(dx, dy) < e:
            break
    print("\n\nSolution: ({0}, {1})".format(x_k1, y_k1))


if __name__ == "__main__":
    print("_______________________SIMPLE_________________\n")
    simple_iteration(1, -1)
    print("\n_______________________NEWTON_1_________________\n")
    newton(-1.6, -0.1)
    print("\n_______________________NEWTON_2_________________\n")
    newton(-0.7, 0.7)
    print("\n_______________________NEWTON_3_________________\n")
    newton(0.7, -0.7)
    print("\n_______________________NEWTON_4_________________\n")
    newton(1.6, 0.1)
