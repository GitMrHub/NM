import math

file = open("Result.txt", "w")

function_global = lambda x: 2*(x**3)-4*(x**2)-x+1
function_first  = lambda x: 6*(x**2)-8*x-1

a1, b1 = (-1, 0)
a2, b2 = (0, 1)
a3, b3 = (2, 3)
e = 0.00001
iterations = 0


def bisector_method(a, b, f, i):
    file.write('-------------Bisector method---------------------\n\n')
    file.write('Interval [' + str(a) + '; ' + str(b) + ']\n')
    file.write(str(e) + '\n')

    file.write('[' + str(a) + '; ' + str(b) + ']  f(' + str(a) + ')=' + str(f(a)) + '  f(' + str(b) + ')=' + str(f(b)) + '\n')
    x = (a+b)/2
    i +=1

    while (math.fabs(f(x)) or math.fabs(b-a)) >= e:
        x = (a+b)/2
        i += 1
        a, b = (a, x) if f(a)*f(x) < 0 else (x, b)
        file.write('[' + str(a) + '; ' + str(b) + ']  f(' + str(a) + ')=' + str(f(a)) + '  f(' + str(b) + ')=' + str(f(b))+ '\n')

    file.write('\nIterations = ' + str(i) + '\n')
    file.write('x = ' + str(x) + '\nf(x) = ' + str(f(x)) + '\n\n')
    return (a+b)/2


def newtons_method(a, b, f, f1, i):
    file.write('-------------Newton method---------------------\n\n')
    file.write('Interval [' + str(a) + '; ' + str(b) + ']\n')
    file.write(str(e) + '\n')

    x0 = (a + b) / 2
    x1 = x0 - (f(x0) / f1(x0))
    i += 1

    file.write('x= ' + str(x1) + '\n')
    file.write('f(x) = ' + str(f(x1)) + '\n')

    while True:
        if (math.fabs(x1 - x0) or math.fabs(f(x1))) < e:
            file.write('\nIterations = ' + str(i) + '\n')
            file.write('x = ' + str(x1) + '\nf(x) = ' + str(function_global(x1)) + '\n\n')
            return x1

        i += 1
        x0 = x1
        x1 = x0 - (f(x0) / f1(x0))
        file.write('x= ' + str(x1) + '\n')
        file.write('f(x) = ' + str(f(x1)) + '\n')


def chord_method(a, b, f, i):
    file.write('-------------Chord method---------------------\n')
    file.write('Interval [' + str(a) + '; ' + str(b) + ']\n')
    file.write(str(e) + '\n')
    x = (a*f(b)-b*f(a))/(f(b)-f(a))
    file.write('x = ' + str(x) + '\nf(x) = ' + str(function_global(x)) + '\n\n')
    i += 1

    while math.fabs(f(x)) >= e:
        x = (a * f(b) - b * f(a)) / (f(b) - f(a))
        i += 1
        a, b = (a, x) if f(a) * f(x) < 0 else (x, b)
        file.write('[' + str(a) + '; ' + str(b) + ']  f(' + str(a) + ')=' + str(f(a)) + '  f(' + str(b) + ')=' + str(f(b))+ '\n')

    file.write('\nIterations = ' + str(i) + '\n')
    file.write('x = ' + str(x) + '\nf(x) = ' + str(f(x)) + '\n\n')
    return (a * f(b) - b * f(a)) / (f(b) - f(a))


bisector_method(a1, b1, function_global, iterations)
newtons_method(a1, b1, function_global, function_first, iterations)
chord_method(a1, b1, function_global, iterations)

bisector_method(a2, b2, function_global, iterations)
newtons_method(a2, b2, function_global, function_first, iterations)
chord_method(a2, b2, function_global, iterations)

bisector_method(a3, b3, function_global, iterations)
newtons_method(a3, b3, function_global, function_first, iterations)
chord_method(a3, b3, function_global, iterations)

file.close()
