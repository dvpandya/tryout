# -*- coding: utf-8 -*-
"""
@author: Darshan
"""
import numpy as np
from scipy import integrate
from matplotlib.pylab import *
from scipy.optimize import fsolve

R=8.314 # J/mol.K
print("CSTR MODEL")

print("The Reaction Considered is ")

print("NH4OH ----------> NH3 + H2O")

H_ammonia=-46110 # in J/mol

print("Enthalpy of Formation of NH3 at 298 K is ",H_ammonia )

H_water=-285830   # in J/mol

print("Enthalpy of Formation of H2O at 298 K is ",H_water )

H_nh4oh=-80800 # in J/mol

print("Enthalpy of Formation of NH4OH at 298 K is ",H_nh4oh )


F1=input("Enter the Input Flow Rate for the First Tank : ")
Cin=input("Enter the Concentration of A in the Input : " )



def tanks_in_series(t, C):
    """
    Reaction : A-----> B + C
    Dynamic balance for a series of 3 CSTRs
 
    C_A1 = C[0] = CA leaving tank 1 [mol/L]
    C_A2 = C[1] = CA leaving tank 2 [mol/L]
    C_A3 = C[2] = CA leaving tank 3 [mol/L]
 
    Returns dy/dt = [F/V*(C_{A,in} - C_{A,out}) - k*C_{A,out}] 
                    repeated for each tank                    
    """
    F = F1            # L/min
    CA_in_1 = Cin       # mol/L
    V1 = V2 = V3 = 500.0  # L     (changed in part 1 and 2)
    k = 0.6            # 1/min
 
    # Assign some variables for convenience of notation
    CA1 = C[0]
    CA2 = C[1]
    CA3 = C[2]
    
    CA_in_2 = CA1
    CA_in_3 = CA2
    
    # Activate this for time greater than 20 mins(disturbance after 20 mins)
    if t > 20:
        CA_in_1 = Cin + 1.0
     
    # Output from ODE function must be a COLUMN vector, with n rows
    n = len(C)      # 2: implies we have two ODEs
    dCdt = np.zeros((n,1))
    dCdt[0] = F/V1*(CA_in_1 - CA1) - k*CA1
    dCdt[1] = F/V2*(CA_in_2 - CA2) - k*CA2
    dCdt[2] = F/V3*(CA_in_3 - CA3) - k*CA3    
    return dCdt

# Set the integrator    
r = integrate.ode(tanks_in_series).set_integrator('vode', method='bdf')

# Set the time range (in min)
t_start = 0.0
t_final = 100.0
delta_t = 0.1
# Number of time steps: 1 extra for initial condition
num_steps = np.floor((t_final - t_start)/delta_t) + 1

# Set initial condition(s): for integrating variable and time
CA1_t_zero = Cin  # mol/L
CA2_t_zero = 0  # mol/L
CA3_t_zero = 0  # mol/L
r.set_initial_value([CA1_t_zero, CA2_t_zero, CA3_t_zero], t_start)

# Additional Python step: create vectors to store trajectories
t = np.zeros((num_steps, 1))
CA1 = np.zeros((num_steps, 1))
CA2 = np.zeros((num_steps, 1))
CA3 = np.zeros((num_steps, 1))

t[0] = t_start
CA1[0] = CA1_t_zero
CA2[0] = CA2_t_zero
CA3[0] = CA3_t_zero

# Integrate the ODE(s) across each delta_t timestep
k = 1
while r.successful() and k < num_steps:
    r.integrate(r.t + delta_t)

    # Store the results to plot later
    t[k] = r.t
    CA1[k] = r.y[0]
    CA2[k] = r.y[1]
    CA3[k] = r.y[2]
    k += 1

# Plot the trajectories:
fig = figure()
plot(t, CA1, 'r')
plot(t, CA2, 'g')
plot(t, CA3, 'b')
xlabel('Time [minutes]')
ylabel('Concentrations [mol/L]')


k = 0.1                   # 1/min
F = 100.0                    # L/min
CA_in_1 = Cin               # mol/L
V = np.arange(50., 1000, 1) # L
tau = F/V

outlet_conversion1=1-tau/(k+tau)
fig = figure()
plot(V, outlet_conversion1)

xlabel('Volume in each CSTR [L]')
ylabel('Conversion for First Tank')

outlet_conversion2 = 1 - np.power(tau/(k+tau),2)
fig = figure()
plot(V, outlet_conversion2)

xlabel('Volume in each CSTR [L]')
ylabel('Conversion for Second Tank')

outlet_conversion3 = 1 - np.power(tau/(k+tau),3)
fig = figure()
plot(V, outlet_conversion3)

plot(V[np.where(V==500)], 
     outlet_conversion3[np.where(V==500)], 'ro', ms=10)
     
xlabel('Volume in each CSTR [L]')
ylabel('Overall Conversion for Third Tank')


H_rxn=H_ammonia+H_water-H_nh4oh

tau_ass=F/500.0

X1=1-tau_ass/(k+tau_ass)

X2=1-np.power(tau_ass/(k+tau_ass),2)

X3=1 - np.power(tau_ass/(k+tau_ass),3)

print("No. of moles of NH4OH decomposed in first tank = ",X1*Cin*F)

print("No. of moles of NH4OH decomposed in second tank = ",X2*CA2[100]*F)

print("No. of moles of NH4OH decomposed in third tank = ",X3*CA3[100]*F)

conv_1=X1*Cin*F

conv_2=X2*CA2[100]*F

conv_3=X3*CA3[100]*F

H_in_1=Cin*F*H_nh4oh
H_out_1=(Cin*F-conv_1)*H_nh4oh+conv_1*H_ammonia+conv_1*H_water

Delta_H_1=H_out_1-H_in_1

print("Change in Enthalpy ",Delta_H_1)

#Cp/R=A+BT+CT^2+DT^-2 for gases
#Cp/R=A+BT+CT^2 for liquids

A_amm=3.578
B_amm=3.020*10**-3
D_amm=-0.186*10**5

A_nh4oh=22.626
B_nh4oh=-100.75*10**-3
C_nh4oh=192.71*10**-6

A_water=8.712
B_water=1.25*10**-3
C_water=-0.18*10**-6

def enthalpy(T):
    f=Delta_H_1+R*((Cin*F-conv_1)*(A_nh4oh*(T-298.15)+B_nh4oh*(T**2-298.15**2)/2+C_nh4oh*(T**3-298.15**3)/3)+ (conv_1)*(A_water*(T-298.15)+B_water*(T**2-298.15**2)/2+C_water*(T**3-298.15**3)/3)+(conv_1)*(A_amm*(T-298.15)+B_amm*(T**2-298.15**2)/2-D_amm*(T**-1-298.15**-1)/2))
    return f

T_soln=fsolve(enthalpy,400)

print("The temp. is " ,T_soln)

H_R=np.zeros(100)
Temp_soln=np.zeros(100)

for i in range(0,100):
    j=i/5
    H_R[i]=j*H_rxn
    def graph(Temp):
        F=H_R[i]+R*(j*(A_nh4oh*(Temp-298.15)+B_nh4oh*(Temp**2-298.15**2)/2+C_nh4oh*(Temp**3-298.15**3)/3)+ i*(A_water*(Temp-298.15)+B_water*(Temp**2-298.15**2)/2+C_water*(Temp**3-298.15**3)/3)+i*(A_amm*(Temp-298.15)+B_amm*(Temp**2-298.15**2)/2-D_amm*(Temp**-1-298.15**-1)/2))
        return F
    
    Temp_soln[i]=fsolve(graph,400)
    i=i+1

fig = figure()
plot(-1*0.001*H_R,Temp_soln)

xlabel('Enthalpy Change in kJ')
ylabel('Temperature of the system in K')