import numpy as np


def f1(x):
    if 0.2 <= x <= 2:
        return np.log(x+1)/(x+1)
    return None


def f1_der(x):
    if 0.2 <= x <= 2:
        return (1 - np.log(x+1)) / (x+1)**2
    return None


def f1_2der(x):
    if 0.2 <= x <= 2:
        return (2*np.log(x+1)-3) / (x+1)**3
    return None


def f1_osc(x):
    if 0.2 <= x <= 2:
        return np.log(x+1)/(x+1) + np.cos(10*x)
    return None


def f1_osc_der(x):
    if 0.2 <= x <= 2:
        return (1 - np.log(x+1)) / (x+1)**2 - 10*np.sin(10*x)
    return None


def f1_osc_2der(x):
    if 0.2 <= x <= 2:
        return (2*np.log(x+1)-3) / (x+1)**3 - 100*np.cos(10*x)
    return None


def f2(x):
    if 2 <= x <= 4:
        return np.log(x+1)/x
    return None


def f2_der(x):
    if 2 <= x <= 4:
        return (x + (-x+1)*np.log(x+1))/(x**3 + x**2)
    return None


def f2_2der(x):
    if 2 <= x <= 4:
        return (-np.log(x+1)*(x**3 + x**2) - 3*x**3 - 2*x**2 +
                3*np.log(x+1)*x**2 * (x+1) + 2*(x+2)*x*np.log(x+1))/(x**3 + x**2)**2
    return None


def f2_osc(x):
    if 2 <= x <= 4:
        return np.log(x+1)/x + np.cos(10*x)
    return None


def f2_osc_der(x):
    if 2 <= x <= 4:
        return (x + (-x+1)*np.log(x+1))/(x**3 + x**2) - 10*np.sin(10*x)
    return None


def f2_osc_2der(x):
    if 2 <= x <= 4:
        return (-np.log(x+1)*(x**3 + x**2) - 3*x**3 - 2*x**2 +
                3*np.log(x+1)*x**2 * (x+1) + 2*(x+2)*x*np.log(x+1))/(x**3 + x**2)**2 - 100*np.cos(10*x)
    return None


def f3(x):
    if 1 <= x <= np.pi:
        return np.sin(x + 1) / x
    return None


def f3_der(x):
    if 1 <= x <= np.pi:
        return (x*np.cos(x+1) - np.sin(x + 1)) / (x**2)
    return None


def f3_2der(x):
    if 1 <= x <= np.pi:
        return (-x**2 * np.sin(x+1) - 2*x*np.cos(x+1) + 2 * np.sin(x + 1)) / (x**3)
    return None


def f3_osc(x):
    if 1 <= x <= np.pi:
        return np.sin(x + 1) / x + np.cos(10*x)
    return None


def f3_osc_der(x):
    if 1 <= x <= np.pi:
        return (x*np.cos(x+1) - np.sin(x + 1)) / (x**2) - 10*np.sin(10*x)
    return None


def f3_osc_2der(x):
    if 1 <= x <= np.pi:
        return (-x**2 * np.sin(x+1) - 2*x*np.cos(x+1) + 2 * np.sin(x + 1)) / (x**3) - 100*np.cos(10*x)
    return None


def phi(x):
    if -1 <= x <= 0:
        return x**3 + 3 * x**2
    elif 0 < x <= 1:
        return -x**3 + 3 * x**2
    else:
        return None


def phi_derivative(x):
    if -1 <= x <= 0:
        return 3*x**2 + 6 * x
    elif 0 < x <= 1:
        return -3*x**2 + 6 * x
    else:
        return None


def phi_2derivative(x):
    if -1 <= x <= 0:
        return 6*x + 6
    elif 0 < x <= 1:
        return -6*x + 6
    else:
        return None


class CubicalSpline(object):
    a: float
    b: float
    c: float
    d: float
    x: tuple

    def __init__(self, x: tuple[float], a: float, b: float, c: float, d: float) -> None:
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.x = x

    def get_func(self):
        def S(x):
            return self.a + self.b*(x-self.x[1]) + (self.c/2)*(x-self.x[1])**2 + (self.d/6)*(x-self.x[1])**3
        return S

    def get_func_derivative(self):
        def Sder(x):
            return self.b + self.c*(x-self.x[1]) + (self.d/2)*(x-self.x[1])**2
        return Sder

    def get_func_2derivative(self):
        def S2der(x):
            return self.c + self.d*(x-self.x[1])
        return S2der


def TDMA(A, B, C, phi, n, mu1, mu2):
    """ реализация прогонки из прошлой л.р."""

    alpha = np.zeros(n+1)
    beta = np.zeros(n+1)
    for i in range(0, n-1):
        alpha[i+1] = B / (-C - A*alpha[i])
        beta[i+1] = (-phi[i+1] + A*beta[i]) / (-C - A*alpha[i])

    v = np.zeros(n+1)
    v[0] = mu1
    v[n] = mu2

    for i in range(n-1, 0, -1):
        v[i] = alpha[i] * v[i+1] + beta[i]

    return v


def find_splines(n, f, a, b, sa, sb) -> list:
    x = np.linspace(a, b, n + 1)
    h = (x[-1]-x[0])/n
    A = h
    B = h
    C = 4*h
    phi = np.zeros(n+1)
    phi[0] = sa
    phi[-1] = sb
    for i in range(1, n):
        phi[i] = 6 * ((f(x[i+1])-f(x[i])) - (f(x[i])-f(x[i-1]))) / h

    c = TDMA(A, B, C, phi, n, sa, sb)
    a = np.zeros(n)
    d = np.zeros(n)
    b = np.zeros(n)

    splines = []
    for i in range(n):
        a[i] = f(x[i+1])
        d[i] = (c[i+1]-c[i])/h
        b[i] = (f(x[i+1]) - f(x[i]))/h + c[i+1]*(h/3) + c[i]*(h/6)
        splines.append(CubicalSpline((x[i],x[i+1]), a[i], b[i], c[i+1], d[i]))

    return splines
