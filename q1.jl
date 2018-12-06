include("./data_struct.jl")

using Data
using JuMP

# Pkg.add("Ipopt")
# using Ipopt
Pkg.add("CPLEX")
using CPLEX
# Pkg.add("Gurobi")
# using Gurobi
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
T_CH4 = 25 + 273.15
T_Air = 25 + 273.15
T_HotFumes = 1600 + 273.15

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

#------------------------ MODEL ---------------------------------------------------------------------
measurements = loadDataFromFile("q1_easy")

m = Model(solver=CplexSolver())
n_Obs = length(measurements.V_Air)



@variable(m, err_CH4_bound >= 0.0)
@variable(m, err_Air_bound >= 0.0)
@variable(m, err_Hot_bound >= 0.0)

@variable(m, V_CH4 >= 0.0)
@variable(m, V_Air >= 0.0)
@variable(m, V_HotFumes >= 0.0)
@variable(m, V_O2 >= 0.0)
@variable(m, V_CO2 >= 0.0)
@variable(m, V_H2O >= 0.0)

@objective(m, Min, err_CH4_bound + err_Air_bound + err_Hot_bound)

#Transformer les volumes d'air et de hotfumes en leur composants
#N2 ne participe pas à la réaction
@constraint(m, V_O2 == V_Air/(1 + (79/21)))
@constraint(m, V_CO2 == V_HotFumes/(1 + 2 + 2*(79/21)))
@constraint(m, V_H2O == (V_HotFumes * 2)/(1 + 2 + 2*(79/21)))

#Mettre l'eqaution en contraintes
@constraint(m, V_CH4/T_CH4 == V_CO2/T_HotFumes)
@constraint(m, V_CH4/T_CH4 == 0.5 * V_H2O/T_HotFumes)
@constraint(m, V_CH4/T_CH4 == 0.5 * V_O2/T_HotFumes)

#Linearisation
@constraint(m, -err_CH4_bound <= measurements.V_NaturalGas[1] - V_CH4 )
@constraint(m, measurements.V_NaturalGas[1] - V_CH4  <= err_CH4_bound)
@constraint(m, -err_Air_bound <= measurements.V_Air[1] - V_Air)
@constraint(m, measurements.V_Air[1] - V_Air <= err_Air_bound)
@constraint(m, -err_Hot_bound <= measurements.V_HotFumes[1] - V_HotFumes)
@constraint(m, measurements.V_HotFumes[1] - V_HotFumes <= err_Hot_bound)

println("The optimization problem to be solved is:")
print(m)

solve(m)

println("Objective value: ", getObjectiveValue(m))
println("err_CH4_bound = ", getValue(err_CH4_bound))
println("err_Air_bound = ", getValue(err_Air_bound))
println("err_Hot_bound = ", getValue(err_Hot_bound))


