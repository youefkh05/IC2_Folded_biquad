import sympy as sy
import matplotlib.pyplot as plt
import control as ct
import math
import numpy as np
import mplcursors as mpl
%matplotlib qt

s_n,s  = sy.symbols('s_n s')

Omega_S = (60-10)/(40-20)
Omega_0 = math.sqrt(40*20) * 2*math.pi*1e6
A_p = -2
As = -20
Q = math.sqrt(40*20)/(40-20)

# Assume Butterworth filter 

epsilon = math.sqrt(10**(2/10) -1)
n_order = math.ceil((20- 20 * math.log10(epsilon))/(20*math.log10(Omega_S)))

poles = [math.pow(1/epsilon,1/n_order)*complex(-1*math.sin((2*i+1)/(2*n_order)*math.pi),math.cos((2*i+1)/(2*n_order)*math.pi)) for i in range(n_order)]
num = np.prod([(pole.real**2 + pole.imag**2) for pole in poles])
num = -np.prod(poles).real
num
#n_order
poles[1] = complex(poles[1].real,0)  # Imaginary number so small 
poles = np.array(poles)
Transfer_function_normalized_control = ct.tf([num], np.poly(poles))
print(Transfer_function_normalized_control)

num_expr = sum(coef * s_n**i for i, coef in enumerate(reversed(Transfer_function_normalized_control.num[0][0])))
den_expr = sum(coef * s_n**i for i, coef in enumerate(reversed(Transfer_function_normalized_control.den[0][0])))

# Create symbolic transfer function
Transfer_function_normalized_sym = num_expr / den_expr
Transfer_function_normalized_sym
sy.init_printing()
Transfer_function_actual_sym = Transfer_function_normalized_sym.evalf(subs={s_n:Q*(s/Omega_0+Omega_0/s)})
sy.simplify(Transfer_function_actual_sym,Rational=True)


n, d = map(lambda p: [float(c) for c in sy.Poly(p, s).all_coeffs()], sy.fraction(sy.simplify(Transfer_function_actual_sym)))
Transfer_function_actual_control = ct.tf(n, d)
print(Transfer_function_actual_control)
Poles_actual = ct.poles(Transfer_function_actual_control)
ct.zeros(Transfer_function_actual_control)
ct.pzmap(Transfer_function_actual_control)
ct.pzmap(Transfer_function_normalized_control)

ct.bode(Transfer_function_actual_control,dB=True,Hz=True)
plt.show()

# Implementation using Sallen-Key 
omega_array,zeta_array,Poles_actual= ct.damp(Transfer_function_actual_control)
Poles_actual
Q_array = 1/(2*zeta_array)
Q_array
omega_array
# Assume K =2 , R1=4*R2 
K = 2;
R1_over_R2 = 4
K_over_R1C5 = np.pow(n[0]/d[0],1/3)
K_over_R1C5
R1C5= K/K_over_R1C5
R1C5
R4C3_array = (1+R1_over_R2)/(omega_array**2 * R1C5)
R4C3_array

R4C5_array = 1/(omega_array/Q_array - 1/R4C3_array + R1_over_R2*(K-1)/R1C5 -1/R1C5)
R4C5_array


# Final Paramters 
# let R4 =10k
C5_array = R4C5_array/10000
C5_array
C3_array = R4C3_array/10000
C3_array
R1_array = R1C5/C5_array
R2_array = R1_array/R1_over_R2
print(R1_array)
R2_array
print(C3_array)
print(C5_array)
Q_array # I re-ordered it to make Lowest Q is the first stage 
