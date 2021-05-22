from sympy import *

κ, K = symbols("kappa K")
U_0, E = symbols("U_0 E", real=True)
a, b = symbols("a b", real=True, positive=True)
k = symbols("k", real=True)

# Using atomic units hbar**2/m=1
dict_subs_E = {
    K: sqrt(2*E),
    κ: sqrt(2*(U_0 - E))
}

dict_num = {
    U_0: 1.0,
    a: 1.0,
    b: 1.0,
    k: 0.0,
}

eq_bands = (κ**2 - K**2)/(2*κ*K) * sinh(κ*b) * sin(K*a) + cosh(κ*b)*cos(K*a) - cos(k*(a+b))

pprint(eq_bands.subs(dict_subs_E))

print(eq_bands.subs(dict_subs_E))

eq_E_k = eq_bands.subs(dict_subs_E).subs(dict_num)
pprint(eq_E_k)

#pprint( solve(eq_E_k, E) )