include("./data_struct.jl")

using Data
using JuMP

Pkg.add("Clp")
using Clp

###################################################

# struct Measurements
#     V_NaturalGas::Array{Float64}
#     V_Air::Array{Float64}
#     V_HotFumes::Array{Float64}
#     wi_NaturalGas::Array{Array{Float64}} # order : CH4, C2H6, C3H8
#     wi_Fumes::Array{Array{Float64}} # order : CO2, H2O, N2
# end
###################################################


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

#------------------------ END OF CONSTANTS --------------------------------------------------------


#------------------------FUNCTIONS --------------------------------------------------------------

#------------------------ MODEL ---------------------------------------------------------------------
measurements = loadDataFromFile("q1_easy")
m = Model(solver=ClpSolver())
n_Obs = length(measurements.V_NaturalGas)
time = 1:n_Obs
@variable(m, err_Air_bound[time] >= 0.0)
@variable(m, err_Hot_bound[time] >= 0.0)
@variable(m, err_CH4_bound[time] >= 0.0)

@variable(m,  V_CH4[time] >= 0.0)
@variable(m,  V_Air[time] >= 0.0)
@variable(m,  V_HotFumes[time] >= 0.0)
@variable(m,  V_O2[time] >= 0.0)
@variable(m,  V_CO2[time] >= 0.0)
@variable(m,  V_H2O[time] >= 0.0)

@objective(m, Min, sum(err_CH4_bound) + sum(err_Air_bound) + sum(err_Hot_bound))

#Transformer les volumes d'air et de hotfumes en leur composants
#N2 ne participe pas à la réaction
@constraint(m, V_O2[time] .== V_Air[time]/(1 + (79/21)))
@constraint(m, V_CO2[time] .== V_HotFumes[time]/(1 + 2 + 2*(79/21)))
@constraint(m, V_H2O[time] .== (V_HotFumes[time] * 2)/(1 + 2 + 2*(79/21)))

# #Mettre l'eqaution en contraintes
@constraint(m, V_CH4[time]/T_CH4 .== V_CO2[time]/T_HotFumes)
@constraint(m, V_CH4[time]/T_CH4 .== 0.5 * V_H2O[time]/T_HotFumes)
@constraint(m, V_CH4[time]/T_CH4 .== 0.5 * V_O2[time]/T_Air)

# @constraint(m,V_CH4/T_CH4 + 2 * V_O2/T_Air .==  V_CO2/T_HotFumes + 2 * V_H2O/T_HotFumes)

#Linearisation
@constraint(m, -err_CH4_bound[time] .<= measurements.V_NaturalGas[time] - V_CH4[time] )
@constraint(m, measurements.V_NaturalGas[time] - V_CH4[time]  .<= err_CH4_bound[time])
@constraint(m, -err_Air_bound[time] .<= measurements.V_Air[time] - V_Air[time])
@constraint(m, measurements.V_Air[time] - V_Air[time] .<= err_Air_bound[time])
@constraint(m, -err_Hot_bound[time] .<= measurements.V_HotFumes[time] - V_HotFumes[time])
@constraint(m, measurements.V_HotFumes[time] - V_HotFumes[time] .<= err_Hot_bound[time])

println("The optimization problem to be solved is:")
print(m)

solve(m)

println("Objective value: ", getobjectivevalue(m))
println("err_CH4_bound = ", getvalue(err_CH4_bound))
println("err_Air_bound = ", getvalue(err_Air_bound))
println("err_Hot_bound = ", getvalue(err_Hot_bound))


