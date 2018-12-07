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
M_C2H6 = 30.06904
M_C3H8 = 44.09562
M_O2 = 31.9988
M_N2 = 28.0134
M_Air = 28.850334
M_CO2 = 44.0095
M_H2O = 18.01528
M_NG = M_CH4 + M_C2H6 + M_C3H8

#Pressions
P_Air = 1
P_Gas = 1
P_HotFumes = 1

#Temperature
T_NG = 25 + 273.15
T_Air = 25 + 273.15
T_HotFumes = 1600 + 273.15

#Coefficients réactions
coeff_CH4 = 1
coeff_CO2_CH4 = 1
coeff_H2O_CH4 = 2

coeff_C2H6 = 2
coeff_CO2_C2H6 = 4
coeff_H2O_C2H6 = 6

coeff_C3H8 = 1
coeff_CO2_C3H8 = 3
coeff_H2O_C3H8 = 4

#------------------------ END OF CONSTANTS --------------------------------------------------------



#------------------------ MODEL ---------------------------------------------------------------------
measurements = loadDataFromFile("q2_easy")

m = Model(solver=ClpSolver())
n_Obs = length(measurements.V_NaturalGas)
time = 1:n_Obs
@variable(m, err_Air_bound[time] >= 0.0)
@variable(m, err_Hot_bound[time] >= 0.0)
@variable(m, err_CH4_bound[time] >= 0.0)

@variable(m,  V_CH4[time] >= 0.0)
@variable(m,  V_C2H6[time] >= 0.0)
@variable(m,  V_C3H8[time] >= 0.0)
@variable(m,  V_HotFumes[time] >= 0.0)
@variable(m,  V_Air[time] >= 0.0)
@variable(m,  V_O2[time] >= 0.0)
@variable(m,  V_CO2[time] >= 0.0)
@variable(m,  V_H2O[time] >= 0.0)


@objective(m, Min, sum(err_CH4_bound) + sum(err_Air_bound) + sum(err_Hot_bound))

#Transformer les volumes d'air et de hotfumes en leur composants
#N2 ne participe pas à la réaction

M_Fumes_Inv = measurements.wi_Fumes[1][time]/M_CO2 + measurements.wi_Fumes[2][time]/M_H2O + measurements.wi_Fumes[3][time]/M_N2

#Proportions fumes POUR 1 MOLE de NG
prop_CO2_CH4 = 1 /(1+2+3)
prop_CO2_C2H6 = 2 /(1+2+3) # 1 mol
prop_CO2_C3H8 = 3 /(1+2+3)
prop_H2O_CH4 =  2/(2+3+4)   
prop_H2O_C2H6 =  3/(2+3+4) # 1 mol
prop_H2O_C3H8 =  4/(2+3+4)

@constraint(m, V_CH4[time]/T_NG .== prop_CO2_CH4 * (coeff_CO2_CH4/coeff_CH4) * measurements.wi_Fumes[1][time] ./(M_CO2 * M_Fumes_Inv[time]) .* V_HotFumes[time]/T_HotFumes)
@constraint(m, V_CH4[time]/T_NG .== prop_H2O_CH4 * (coeff_H2O_CH4/coeff_CH4) * measurements.wi_Fumes[2][time] ./(M_H2O * M_Fumes_Inv[time]) .* V_HotFumes[time]/T_HotFumes)
# 3rd constraint missing (q1 still wrong)

@constraint(m, V_C2H6[time]/T_NG .== prop_CO2_C2H6 * (coeff_CO2_C2H6/coeff_C2H6) * measurements.wi_Fumes[1][time] ./(M_CO2 * M_Fumes_Inv[time]) .* V_HotFumes[time]/T_HotFumes)
@constraint(m, V_C2H6[time]/T_NG .== prop_H2O_C2H6 * (coeff_H2O_C2H6/coeff_C2H6) * measurements.wi_Fumes[2][time] ./(M_H2O * M_Fumes_Inv[time]) .* V_HotFumes[time]/T_HotFumes)
# 3rd constraint missing (q1 still wrong)

@constraint(m, V_C3H8[time]/T_NG .== prop_CO2_C3H8 * (coeff_CO2_C3H8/coeff_C3H8) * measurements.wi_Fumes[1][time] ./(M_CO2 * M_Fumes_Inv[time]) .* V_HotFumes[time]/T_HotFumes)
@constraint(m, V_C3H8[time]/T_NG .== prop_H2O_C3H8 * (coeff_H2O_C3H8/coeff_C3H8) * measurements.wi_Fumes[2][time] ./(M_H2O * M_Fumes_Inv[time]) .* V_HotFumes[time]/T_HotFumes)
# 3rd constraint missing (q1 still wrong)


# @constraint(m, V_CH4[time] .== measurements.wi_NaturalGas[t][1] * (M_CH4/M_NG) * measurements.V_NaturalGas[t])
# @constraint(m, V_C2H6[time] .== measurements.wi_NaturalGas[t][2] * (M_C2H6/M_NG) * measurements.V_NaturalGas[t])
# @constraint(m, V_C3H8[time] .== measurements.wi_NaturalGas[t][3] * (M_C3H8/M_NG) * measurements.V_NaturalGas[t])

# @constraint(m, V_O2[time] .== V_Air[time]/(1 + (79/21)))
# @constraint(m, V_CO2[time] .== V_HotFumes[time]/(1 + 2 + 2*(79/21)))
# @constraint(m, V_H2O[time] .== (V_HotFumes[time] * 2)/(1 + 2 + 2*(79/21)))


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
println(getvalue(err_CH4_bound))
println(getvalue(err_Air_bound))
println(getvalue(err_Hot_bound))
println(getvalue(V_CH4))
println(getvalue(V_Air))
println(getvalue(V_HotFumes))
