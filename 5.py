from math import sqrt
import matplotlib.pyplot as plt

a = -0.7
b = 0.7


def polinom(x):
    return (x ** 2) / (2 * sqrt(1 - 3 * (x ** 4)))


def deriv_2(x):
    return 3 * x ** 5 / (1 - 3 * x ** 4) ** (3 / 2) + x / sqrt(1 - 3 * x ** 4)


def deriv_6(x):
    return 540 * (-3780 * x ** 12 / (3 * x ** 4 - 1) ** 3 + 1620 * x ** 8 / (3 * x ** 4 - 1) ** 2 - 153 * x ** 4 / (
                3 * x ** 4 - 1) - 36 * x ** 4 * (2268 * x ** 12 / (3 * x ** 4 - 1) ** 3 - 1260 * x ** 8 / (
                3 * x ** 4 - 1) ** 2 + 195 * x ** 4 / (3 * x ** 4 - 1) - 7) / (1 - 3 * x ** 4) + 63 * x ** 4 * (
                              7128 * x ** 16 / (3 * x ** 4 - 1) ** 4 - 4860 * x ** 12 / (
                                  3 * x ** 4 - 1) ** 3 + 1050 * x ** 8 / (3 * x ** 4 - 1) ** 2 - 75 * x ** 4 / (
                                          3 * x ** 4 - 1) + 1) / (1 - 3 * x ** 4) + 1) / (1 - 3 * x ** 4) ** (3 / 2)


def newton(x, xi):
    s = polinom(xi[0])

    for i in range(1, len(xi)):
        F = 0
        for j in range(i + 1):
            z = 1.0
            for k in range(i + 1):
                if k != j:
                    z = z * (xi[j] - xi[k])
            F += polinom(xi[j]) / z

        for k in range(i):
            F *= (x - xi[k])

        s = s + F
    return s*10


def lagrange(x, xi):
    s = 0
    for i in range(len(xi)):
        w = 1
        for j in range(len(xi)):
            if i != j:
                w = w * (x - xi[j]) / (xi[i] - xi[j])
            else:
                w = 0
        s = s + w * polinom(xi[i])
    return s


def majoranta(t, a, b, xi):
    if (a >= b):
        a, b = b, a
    if (a > 0):
        sup = abs(deriv_6(a))
    elif (a < 0):
        sup = abs(deriv_6(b))
    w = 1
    for i in range(len(xi)):
        w *= (t - xi[i])
    return (sup)*abs(w)

def spline(t, x, d2f):
    h = []
    l = []
    dlt = []
    lmb = []
    a = []
    b = []
    c = []
    d = []
    y = []
    for i in range(len(x)):
        y.append(polinom(x[i]))
        h.append(0)
        l.append(0)
        dlt.append(0)
        lmb.append(0)
        a.append(y[i])
        b.append(0)
        c.append(0)
        d.append(0)

    for k in range(1, len(x)):
        h[k] = x[k] - x[k - 1]
        l[k] = (y[k] - y[k - 1]) / h[k]

    dlt[1] = -h[2] / (2 * (h[1] + h[2]))
    lmb[1] = 1.5 * (l[2] - l[1]) / (h[1] + h[2])

    for k in range(3, len(x)):
        dlt[k - 1] = -h[k] / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * dlt[k - 2])
        lmb[k - 1] = (3 * l[k] - 3 * l[k - 1] - h[k - 1] * lmb[k - 2]) / (
                2 * h[k - 1] + 2 * h[k] + h[k - 1] * dlt[k - 2])

    c[0] = d2f(x[0]) / 2
    c[len(x) - 1] = d2f(x[-1]) / 2
    for k in range(len(x) - 1, 1, -1):
        c[k - 1] = dlt[k - 1] * c[k] + lmb[k - 1]
    for k in range(1, len(x)):
        d[k] = (c[k] - c[k - 1]) / (3 * h[k])
        b[k] = l[k] + (2 * c[k] * h[k] + h[k] * c[k - 1]) / 3

    res = 0.0
    for k in range(1, len(x)):
        if (x[k - 1] <= t <= x[k]):
            res = a[k] + b[k] * (t - x[k]) + c[k] * (t - x[k]) ** 2 + d[k] * (t - x[k]) ** 3
    return res*10


xi = []
xr = []
yi = []

i = a
while i <= b:
    xi.append(i)
    i += 0.1

xi.reverse()
xr.extend(xi)
xi.reverse()

x = xi[0]
dx = 0.01
domain = []

while not (x >= b + dx):
    domain.append(x)
    x = round(x + dx, 2)

for i in range(len(domain)):
    yi.append(polinom(domain[i]))

y_l, y_n, y_nb, y_m, y_s = [], [], [], [], []

for i in range(len(domain)):
    y_l.append(lagrange(domain[i], xi))
    y_n.append(newton(domain[i], xi))
    y_nb.append(newton(domain[i], xr))
    y_m.append(majoranta(domain[i], a, b, xi))
    y_s.append(spline(domain[i], xi, deriv_2))

# ______abs(error)____#
y_l_f, y_n_f, y_nb_f, y_s_f, = [], [], [], []
for i in range(len(domain)):
    y_l_f.append(abs(y_l[i] - yi[i]))
    y_n_f.append(abs(y_n[i] - yi[i]))
    y_nb_f.append(abs(y_nb[i] - yi[i]))
    y_s_f.append(abs(y_s[i] - yi[i]))

"""____________Graphics____________"""
# grid1 = plt.grid(True)
# graph1 = plt.plot(domain, y_l_f, label = 'Lagrange')
# graph2 = plt.plot(domain, y_n_f, label = 'Newton')
# graph2_2 = plt.plot(domain, y_nb_f, label = 'NewtonBack')
# graph3 = plt.plot(domain, y_s_f, label = 'Spline')
# graph4 = plt.plot(domain, y_m, label = 'Majoranta')
# plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
#      ncol=2, mode="expand", borderaxespad=0.)

fig, axes = plt.subplots(1, 5, figsize=(10, 5))
axes[0].plot(domain, y_l_f, 'g')
axes[0].grid(True)
axes[0].set_title('Lagrange')
axes[1].plot(domain, y_n_f, 'b')
axes[1].grid(True)
axes[1].set_title('Newton')
axes[2].plot(domain, y_nb_f, 'y')
axes[2].set_title('Newton Back')
axes[2].grid(True)
axes[3].plot(domain, y_s_f, 'm')
axes[3].grid(True)
axes[3].set_title('Spline')
axes[4].plot(domain, y_m, 'r')
axes[4].grid(True)
axes[4].set_title('Majoranta')
fig.tight_layout()

plt.show()
