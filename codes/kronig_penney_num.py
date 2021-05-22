import numpy as np
import matplotlib.pyplot as plt

# Parameters (global)
U_0 = 10.0
a = 10.0
b = 0.01
#k = 0.1
k = np.pi/(a+b)

def my_eq_E(E):
    SQRT2 = np.sqrt(2.0)
    return -np.cos(k*(a + b)) + \
    np.cos(SQRT2*np.sqrt(E)*a)*np.cosh(b*np.sqrt(-2*E + 2*U_0)) + \
    SQRT2*(-4*E + 2*U_0)*np.sinh(b*np.sqrt(-2*E + 2*U_0))/(4*np.sqrt(E)*np.sqrt(-2*E + 2*U_0))

E = np.linspace(1e-3, 2.0, 2000)
f = my_eq_E(E)

plt.clf()
plt.plot(E, f)
plt.grid(True)
plt.savefig("IMG_k_0.pdf")