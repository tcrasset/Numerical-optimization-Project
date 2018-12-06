include("./data_struct.jl")

using Data
using JuMP

Pkg.add("Ipopt")
using Ipopt
# Pkg.add("CPLEX")
# using CPLEX

###################################################

# struct Measurements
#     V_NaturalGas::Array{Float64}
#     V_Air::Array{Float64}
#     V_HotFumes::Array{Float64}
#     wi_NaturalGas::Array{Array{Float64}} # order : CH4, C2H6, C3H8
#     wi_Fumes::Array{Array{Float64}} # order : CO2, H2O, N2
# end


#------------------------ CONSTANTS --------------------------------------------------------

#Masses molaires
M_CH4 = 16.04246
M_O2 = 31.9988
M_N2 = 28.0134
M_Air = 28.850334
M_CO2 = 44.0095
M_H2O = 18.01528

#Pressions
P_Air = 1
P_CH4 = 1
P_HotFumes = 1

#Temperature
T_CH4 = 25
T_Air = 25
T_HotFumes = 1600

#Coefficients réactions
coeff_CH4 = 1
coeff_O2 = 2
coeff_N2 = 2 * (79/21)
coeff_Air = coeff_O2 + coeff_N2
coeff_CO2 = 1
coeff_H2O = 2

# #Fractions molaires réactifs
# x_CH4 = coeff_CH4/n_reactifs_tots
# x_Air = coeff_Air/ n_reactifs_tots
# x_O2 = coeff_02/ n_reactifs_tots
# x_N2 = coeff_N2/n_reactifs_tots

# #Fractions molaires produits
# x_C02 = coeff_C02/n_produits_tots
# x_H20 = coeff_H2O/ n_produits_tots
# x_N2 = coeff_N2/n_produits_tots






#------------------------ END OF CONSTANTS --------------------------------------------------------


#------------------------FUNCTIONS --------------------------------------------------------------

function V_to_N_gaslaw(P::Int64, T::Int64,V::Float64)
    """Params 
        P - pressure in atm
        T - temperature in Celsius
        V - volume in m3
    """
    return (P * 1000 * V)/(0.08205 * (273.15 + T))
end


#------------------------ MODEL ---------------------------------------------------------------------
measurements = loadDataFromFile("q1_easy")

print(measurements.V_NaturalGas)
print(measurements.V_Air)
print(measurements.V_HotFumes)
print(measurements.wi_NaturalGas)
print(measurements.wi_Fumes)

length(measurements.V_Air)

#Solve for wi_Air




m = Model(solver=IpoptSolver())
n_Obs = length(measurements.V_Air)
@objective(min (sum(erreur)))
# for t = 1:n_Obs


    @variable(m, err_CH4)
    @variable(m, err_Air)
    @variable(m, V_CH4_measured >= 0)
    @variable(m, V_Air_measured >= 0)
    @variable(m, V_HotFumes_measured >= 0)

    #The measure of V is imprecise
    @constraint(V_CH4_measured == measurements.V_NaturalGas[1] + err_CH4)
    @constraint(V_Air_measured == measurements.V_Air[1] + err_Air)
    @constraint(V_HotFumes_measured == measurements.V_HotFumes[1] + err_Air)


    
    n_CH4_measured = V_to_N_gaslaw(P_CH4,T_CH4,V_CH4_measured)
    n_Air_measured = V_to_N_gaslaw(P_Air,T_Air,V_Air_measured)
    
    
    n_O2_measured = 21/100 * n_Air_measured
    n_N2_measured = 79/100 * n_Air_measured
    
    n_reactifs_tots = n_Air_measured + n_CH4_measured
    
    x_O2 = n_O2_measured/n_reactifs_tots
    x_N2 = n_N2_measured/n_reactifs_tots
    x_CH4 = n_CH4_measured/n_reactifs_tots
    
    
    n_produits_tots = n_CH4_measured * (coeff_CO2 + coeff_H2O + coeff_N2)
    
    x_CO2 = coeff_CO2/n_produits_tots
    x_H2O = coeff_H2O/ n_produits_tots
    x_N2 = coeff_N2/n_produits_tots
    
    # #Masse molaire moyenne
    M_reactifs = x_CH4 * M_CH4 + x_N2 * M_N2 + x_O2 * M_O2
    M_produits = x_C02 * M_CO2 + x_H20 * M_H20 + x_N2 * M_N2
    
    @contraint(measurements.wi_Fumes[1][t] == x_CO2 * (M_CO2/M_produits))
    @contraint(measurements.wi_Fumes[2][t] == x_H2O * (M_H2O/M_produits))
    @contraint(measurements.wi_Fumes[3][t] == x_N2 * (M_N2/M_produits))

    @contraint(measurements.wi_NaturalGas[1][1] == x_CH4 * (M_CH4/M_reactifs))

    
#     solve(m)
# end