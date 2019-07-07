"""
1D PFR: concentration profile of one component system
"""
# Importing the package for computations in Python
import numpy as np
# Importing the package for plotting in Python
import matplotlib.pyplot as plt

def main():
    """
    1. Defining the following types of variables:
        - Reactor variables
        - Variables for concentration profile
        - Variables for temperature profile
        - Computational stencil parameters
    """
    # Reactor variables
    L_reactor = 1.0     	    # Length of reactor (m)
    velocity_inlet = 1.0        # Velocity of the entering reactant gas mixture (m/s)
    
    # Variables for concentration profile
    c_A_in = 1.0                # Inlet feed concentration of A (mol/m3)
    c_B_in = 0.0                # Inlet feed concentration of B (mol/m3)
    
    D_f_A = 5e-4                # Diffusion coefficient of reactant A in gas (m2/s)
    D_f_B = 5e-4                # Diffusion coefficient of reactant B in gas (m2/s)

    T_reactor = 500.0           # Reactor temperature (K)

    # Kinetic parameters
    Ea = 1.6e5                  # Activation energy (J/mol)
    R = 8.314                   # Gas constant (J/(K*mol))
    k0 = 3.98e17                # Arrhenius pre-factor (1/s)3.98e16
    
    # Combining the variables
    ingoing = [c_A_in, c_B_in, T_reactor]
    D_f = [D_f_A, D_f_B]
    
    # Computational stencil parameters
    delta_t = 1e-3              # Time-step value (s)
    delta_x = 5e-3              # Grid size for reactor transport (m)
    
    """
    2. Discretization of space and determination of the total time duration
    """
    time = 1  
    timesteps = int(time/delta_t)

    spacesteps = int(L_reactor/delta_x)         # Number of steps in space
    x = np.linspace(0,L_reactor,spacesteps+1)   # Vector with the steps in space
    
    """
    3. Calculating the numerical solution
    """
    # Initial concentration profile (at t=0, only inert gas is present in the reactor)
    conc_A_current = np.ones(spacesteps+1)
    conc_B_current = np.zeros(spacesteps+1)

    current_time = 0
    
    molfraction_A_0 = conc_A_current / (conc_A_current + conc_B_current)
    molfraction_B_0 = conc_B_current / (conc_A_current + conc_B_current)
    
    # Numerical solutions at several points in time
    while current_time < 0.01:
        y, conc_A_current, conc_B_current =  solver(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, conc_A_current, conc_B_current, T_reactor, ingoing, k0, Ea, R)
        current_time = current_time + delta_t
        molfraction_A_0_01 = conc_A_current / (conc_A_current + conc_B_current)
        molfraction_B_0_01 = conc_B_current / (conc_A_current + conc_B_current)
    
    while current_time < 0.1:
        y, conc_A_current, conc_B_current =  solver(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, conc_A_current, conc_B_current, T_reactor, ingoing, k0, Ea, R)
        current_time = current_time + delta_t
        molfraction_A_0_1 = conc_A_current / (conc_A_current + conc_B_current)
        molfraction_B_0_1 = conc_B_current / (conc_A_current + conc_B_current)
    
    while current_time < 0.5:
        y, conc_A_current, conc_B_current =  solver(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, conc_A_current, conc_B_current, T_reactor, ingoing, k0, Ea, R)
        current_time = current_time + delta_t
        molfraction_A_0_5 = conc_A_current / (conc_A_current + conc_B_current)
        molfraction_B_0_5 = conc_B_current / (conc_A_current + conc_B_current)
        
    while current_time < time:
        y, conc_A_current, conc_B_current =  solver(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, conc_A_current, conc_B_current, T_reactor, ingoing, k0, Ea, R)
        current_time = current_time + delta_t
        molfraction_A_1 = conc_A_current / (conc_A_current + conc_B_current)
        molfraction_B_1 = conc_B_current / (conc_A_current + conc_B_current)
    
    # Analytical solution for concentration A
    analytical_conc_A = np.exp(-(k0 * np.exp(-Ea/(R*T_reactor))*x))
    
    """
    4. Plotting the numerical and analytical solution
    """
    fig = plt.subplots(2,1)
    fig = plt.subplots_adjust(hspace=1.0)
    plt.grid()
    plt.suptitle('1D PFR: molfraction of reactants in 1D PFR with one component system')
    plt.subplot(2,1,1)
    plt.plot(x,molfraction_A_0, label='Numerical solution of concentration A at t = 0 second')
    plt.plot(x,molfraction_A_0_01, label='Numerical solution of concentration A at t = 0.01 second')
    plt.plot(x,molfraction_A_0_1, label='Numerical solution of concentration A at t = 0.1 second')
    plt.plot(x,molfraction_A_0_5, label='Numerical solution of concentration A at t = 0.5 second')
    plt.plot(x,molfraction_A_1, label='Numerical solution of concentration A at t = 1 second')
    plt.plot(x,analytical_conc_A, label='Analytical steady state solution of concentration A')
    plt.xlabel('Position in reactor [m]')
    plt.ylabel('Molfraction [-]')
    plt.title('Reactant A')
    plt.legend(loc = 'center left', bbox_to_anchor = (1.0, 0.5))
    plt.subplot(2,1,2)
    plt.plot(x,molfraction_B_0, label='Numerical solution of concentration B at t = 0 second')
    plt.plot(x,molfraction_B_0_01, label='Numerical solution of concentration B at t = 0.01 second')
    plt.plot(x,molfraction_B_0_1, label='Numerical solution of concentration B at t = 0.1 second')
    plt.plot(x,molfraction_B_0_5, label='Numerical solution of concentration B at t = 0.5 second')
    plt.plot(x,molfraction_B_1, label='Numerical solution of concentration B at t = 1 second')
    plt.xlabel('Position in reactor [m]')
    plt.ylabel('Molfraction [-]')
    plt.title('Reactant B')
    plt.legend(loc = 'center left', bbox_to_anchor = (1.0, 0.5))
    filename = '1D_PFR_conc_profile_one_component.png'
    plt.savefig(filename)
    
    """
    6. Functions needed for solving the concentration profile numerical
    """
    
def solver(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, conc_A_current, conc_B_current, T_reactor, ingoing, k0, Ea, R):
    # Generate matrices
    matrix_A = generate_matrix_concentration(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, conc_A_current, conc_B_current, T_reactor, ingoing, 0)
    matrix_B = generate_matrix_concentration(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, conc_A_current, conc_B_current, T_reactor, ingoing, 1)
#   matrix_T, rho_Cp_T = generate_matrix_temperature(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, k_f, rho, Cp, conc_A_current, conc_B_current, T_current, ingoing, 2)
  
    # Generate right-hand side vector
    vector_A = np.copy(conc_A_current)
    vector_B = np.copy(conc_B_current)
#   vector_T = np.copy(T_current)
        
    # Process the reaction term within the right-hand side vector
    for i in range(1,spacesteps):
#           heat_transfer = - (U * a * (vector_T[i]-T_wall) * delta_t )/(rho_Cp_T[i])
            reaction = (vector_A[i]*delta_t) * -(k0 * np.exp(-Ea/(R*T_reactor)))
            vector_A[i] = vector_A[i] + reaction
            vector_B[i] = vector_B[i] - reaction
#           vector_T[i] = vector_T[i] + (reaction*delta_H*(delta_t/(rho_Cp_T[i]))) + heat_transfer
    # Implement the boundary condition within the right-hand side vector
    vector_A[0] = ingoing[0]
    vector_B[0] = ingoing[1]
#   vector_T[0] = ingoing[2]
    
    vector_A[spacesteps] = 0
    vector_B[spacesteps] = 0
#   vector_T[spacesteps] = 0
    
    # Solve the equation Ax=B with x being the concentration profile 
    # at the next timestep
    A = np.linalg.solve(matrix_A,vector_A)
    B = np.linalg.solve(matrix_B,vector_B)
#   T = np.linalg.solve(matrix_T,vector_T)    
    
    return np.linspace(0,L_reactor,spacesteps+1), A, B
    
def generate_matrix_concentration(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, conc_A_current, conc_B_current, T_reactor, ingoing, situation):
    # Define step
    step = L_reactor/(spacesteps+1)

    # Initialization of matrix
    matrix = np.zeros((spacesteps+1,spacesteps+1))
        
    # Loop over position
    for i in range(1,spacesteps):
        # Calculating the current velocity
        velocity_current = velocity(velocity_inlet, L_reactor, conc_A_current[i], conc_B_current[i], T_reactor, ingoing)
        # Combining the variables
        alpha = -(delta_t*D_f[situation])/(step**2)
        beta = (velocity_current*delta_t)/(2.0*step)
        # Generating the matrix
        matrix[i,i-1] = alpha - beta
        matrix[i,i+1] = alpha + beta
        matrix[i,i] = 1.0 - 2.0*alpha
        
    # Boundary conditions
    matrix[0,0] = 1.0
    matrix[spacesteps,spacesteps-3] = -2.0 / (6.0*step)
    matrix[spacesteps,spacesteps-2] = 9.0 / (6.0*step)
    matrix[spacesteps,spacesteps-1] = -18.0 / (6.0*step)
    matrix[spacesteps,spacesteps] = 11.0 / (6.0*step)

    return matrix

#def generate_matrix_temperature(L_reactor, velocity_inlet, spacesteps, delta_t, D_f, k_f, rho, Cp, conc_A_current, conc_B_current, T_current, ingoing, situation):
#    # Define step
#    step = L_reactor/(spacesteps+1)
#
#    # Initialization of matrix
#    matrix = np.zeros((spacesteps+1,spacesteps+1))
#    rho_Cp_T = np.zeros(spacesteps)
#
#    # Loop over position
#    for i in range(1,spacesteps):
#        molfraction = [0,0,0,0]
#        # Determing average density and specific heat capacity depending on the molfraction of the reactant
#        molfraction[0] = conc_A_current[i]/(conc_A_current[i] + conc_B_current[i])
#        molfraction[1] = conc_B_current[i]/(conc_A_current[i] + conc_B_current[i])
#        rho_Cp_T[i] = ((molfraction[0]*rho[0]*Cp[0])+(molfraction[1]*rho[1]*Cp[1]))
#        # Calculating the current velocity
#        velocity_current = velocity(velocity_inlet, L_reactor, conc_A_current[i], conc_B_current[i], T_current[i], ingoing, D_f, k_f)
#        # Combining the variables
#        alpha = -(delta_t * k_f)/(rho_Cp_T[i]*step**2)
#        beta = (velocity_current*delta_t)/(2.0*step) 
#        # Generating the matrix
#        matrix[i,i-1] = alpha - beta
#        matrix[i,i+1] = alpha + beta
#        matrix[i,i] = (1 - 2*alpha)
#        
#    # Boundary conditions
#    matrix[0,0] = 1.0
#    matrix[spacesteps,spacesteps-3] = -2.0 / (6.0*step)
#    matrix[spacesteps,spacesteps-2] = 9.0 / (6.0*step)
#    matrix[spacesteps,spacesteps-1] = -18.0 / (6.0*step)
#    matrix[spacesteps,spacesteps] = 11.0 / (6.0*step)
#
#    return matrix, rho_Cp_T
 
def velocity(velocity_inlet, L_reactor, conc_A, conc_B, temp, ingoing):
    # Defining the variables that aren't used before
    pressure_ingoing = 20.0
    pressure = pressure_ingoing
    
    temp_ingoing = temp
    
    # Calculating the current velocity in the reactor
    velocity_current = velocity_inlet * (conc_A+conc_B)/(ingoing[0]+ingoing[1]) * (pressure_ingoing)/(pressure) * (temp)/(temp_ingoing)
    
    return velocity_current
 
if __name__ == '__main__':
    main()