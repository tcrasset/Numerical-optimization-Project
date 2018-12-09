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

#Temperature
T_NG = 25 + 273.15
T_Air = 25 + 273.15
T_HotFumes = 1600 + 273.15

#Coefficients réactions
coeff_CH4 = 1
coeff_CO2_CH4 = 1
coeff_H2O_CH4 = 2
coeff_O2_CH4 = 2

coeff_C2H6 = 2
coeff_CO2_C2H6 = 4
coeff_H2O_C2H6 = 6
coeff_O2_C2H6 = 7

coeff_C3H8 = 1
coeff_CO2_C3H8 = 3
coeff_H2O_C3H8 = 4
coeff_O2_C3H8 = 5


#------------------------ END OF CONSTANTS --------------------------------------------------------



#------------------------ MODEL ---------------------------------------------------------------------
measurements = loadDataFromFile("q3")

m = Model(solver=ClpSolver())
n_Obs = length(measurements.V_NaturalGas)
time = 1:n_Obs
@variable(m, err_V_NG_bound[time] >= 0.0)
@variable(m, err_V_Air_bound[time] >= 0.0)
@variable(m, err_V_Hot_bound[time] >= 0.0)

@variable(m, err_w_CH4_bound[time] >= 0.0)
@variable(m, err_w_C2H6_bound[time] >= 0.0)
@variable(m, err_w_C3H8_bound[time] >= 0.0)

@variable(m, err_w_CO2_bound[time] >= 0.0)
@variable(m, err_w_H2O_bound[time] >= 0.0)
@variable(m, err_w_O2_bound[time] >= 0.0)

@variable(m, err_w_CH4[time] >= 0.0)
@variable(m, err_w_C2H6[time] >= 0.0)
@variable(m, err_w_C3H8[time] >= 0.0)

@variable(m, err_w_CO2[time] >= 0.0)
@variable(m, err_w_H2O[time] >= 0.0)
@variable(m, err_w_N2[time] >= 0.0)

@variable(m,  err_V_NG[time] >= 0.0)
@variable(m,  err_V_HotFumes[time] >= 0.0)
@variable(m,  err_V_Air[time] >= 0.0)


@objective(m, Min, sum(err_V_NG_bound) + sum(err_V_Air_bound) + sum(err_V_Hot_bound)
                + sum(err_w_CH4_bound) + sum(err_w_C2H6_bound) + sum(err_w_C3H8_bound)
                + sum(err_w_CO2_bound) + sum(err_w_H2O_bound))

#Transformer les volumes d'air et de hotfumes en leur composants
#N2 ne participe pas à la réaction

# M_Fumes_Inv = measurements.wi_Fumes[1][time]/M_CO2 + measurements.wi_Fumes[2][time]/M_H2O + measurements.wi_Fumes[3][time]/M_N2
# M_NG_Inv = measurements.wi_NaturalGas[1][time]/M_CH4 + measurements.wi_NaturalGas[2][time]/M_C2H6 + measurements.wi_NaturalGas[3][time]/M_C3H8
# M_Air_Inv = 1/(0.21 * M_O2 + 0.79 * M_N2)

#Proportions fumes POUR 1 MOLE de NG
prop_CO2_CH4 = 1
prop_H2O_CH4 = 2    
prop_O2_CH4 = 2

prop_CO2_C2H6 = 2  # 1 mol
prop_H2O_C2H6 = 3 # 1 mol
prop_O2_C2H6 = 3.5

prop_CO2_C3H8 = 3 
prop_H2O_C3H8 = 4
prop_O2_C3H8 = 5


#Constraint de CO2
@constraint(m, measurements.V_NaturalGas[time]/T_NG .* (  (1 + err_w_CO2[time] + err_w_CH4[time] + err_V_NG[time])  .*(prop_CO2_CH4 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[1][time]/M_CO2/M_CH4)
                                                        + (1 + err_w_CO2[time] + err_w_C2H6[time] + err_V_NG[time]) .*(prop_CO2_C2H6 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[2][time]/M_CO2/M_C2H6)      
                                                        + (1 + err_w_CO2[time] + err_w_C3H8[time] + err_V_NG[time]) .*(prop_CO2_C3H8 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[3][time]/M_CO2/M_C3H8)
                                                        + (1 + err_w_H2O[time] + err_w_CH4[time] + err_V_NG[time])  .*(prop_CO2_CH4 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[1][time]/M_H2O/M_CH4)
                                                        + (1 + err_w_H2O[time] + err_w_C2H6[time] + err_V_NG[time]) .*(prop_CO2_C2H6 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[2][time]/M_H2O/M_C2H6)      
                                                        + (1 + err_w_H2O[time] + err_w_C3H8[time] + err_V_NG[time]) .*(prop_CO2_C3H8 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[3][time]/M_H2O/M_C3H8)
                                                                                
                                                        + (1 + err_w_N2[time] + err_w_CH4[time] + err_V_NG[time])  .*(prop_CO2_CH4 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[1][time]/M_N2/M_CH4)
                                                        + (1 + err_w_N2[time] + err_w_C2H6[time] + err_V_NG[time]) .*(prop_CO2_C2H6 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[2][time]/M_N2/M_C2H6)      
                                                        + (1 + err_w_N2[time] + err_w_C3H8[time] + err_V_NG[time]) .*(prop_CO2_C3H8 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[3][time]/M_N2/M_C3H8)
                                                        )
                                                    
                                                    .==
                measurements.V_HotFumes[time]/M_CO2/T_HotFumes .*( (1/M_CH4) * measurements.wi_NaturalGas[1][time] .*(1 + err_w_CH4[time] + err_V_HotFumes[time] + err_w_CO2[time])
                                                                +  (1/M_C2H6) * measurements.wi_NaturalGas[2][time] .*(1 + err_w_C2H6[time] + err_V_HotFumes[time] + err_w_CO2[time])
                                                                +  (1/M_C3H8) * measurements.wi_NaturalGas[3][time] .*(1 + err_w_C3H8[time] + err_V_HotFumes[time] + err_w_CO2[time])
                                                                )
)

#Constraint de H2O
@constraint(m, measurements.V_NaturalGas[time]/T_NG .* (  (1 + err_w_CO2[time] + err_w_CH4[time] + err_V_NG[time])  .*(prop_H2O_CH4 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[1][time]/M_CO2/M_CH4)
                                                        + (1 + err_w_CO2[time] + err_w_C2H6[time] + err_V_NG[time]) .*(prop_H2O_C2H6 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[2][time]/M_CO2/M_C2H6)      
                                                        + (1 + err_w_CO2[time] + err_w_C3H8[time] + err_V_NG[time]) .*(prop_H2O_C3H8 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[3][time]/M_CO2/M_C3H8)
                                                        + (1 + err_w_H2O[time] + err_w_CH4[time] + err_V_NG[time])  .*(prop_H2O_CH4 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[1][time]/M_H2O/M_CH4)
                                                        + (1 + err_w_H2O[time] + err_w_C2H6[time] + err_V_NG[time]) .*(prop_H2O_C2H6 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[2][time]/M_H2O/M_C2H6)      
                                                        + (1 + err_w_H2O[time] + err_w_C3H8[time] + err_V_NG[time]) .*(prop_H2O_C3H8 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[3][time]/M_H2O/M_C3H8)
                                                                                
                                                        + (1 + err_w_N2[time] + err_w_CH4[time] + err_V_NG[time])  .*(prop_H2O_CH4 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[1][time]/M_N2/M_CH4)
                                                        + (1 + err_w_N2[time] + err_w_C2H6[time] + err_V_NG[time]) .*(prop_H2O_C2H6 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[2][time]/M_N2/M_C2H6)      
                                                        + (1 + err_w_N2[time] + err_w_C3H8[time] + err_V_NG[time]) .*(prop_H2O_C3H8 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[3][time]/M_N2/M_C3H8)
                                                        )
                                                    
                                                    .==
                measurements.V_HotFumes[time]/M_H2O/T_HotFumes .*( (1/M_CH4) * measurements.wi_NaturalGas[1][time] .*(1 + err_w_CH4[time] + err_V_HotFumes[time] + err_w_H2O[time])
                                                                +  (1/M_C2H6) * measurements.wi_NaturalGas[2][time] .*(1 + err_w_C2H6[time] + err_V_HotFumes[time] + err_w_H2O[time])
                                                                +  (1/M_C3H8) * measurements.wi_NaturalGas[3][time] .*(1 + err_w_C3H8[time] + err_V_HotFumes[time] + err_w_H2O[time])
                                                                )
)

#Constraint on O2
@constraint(m, measurements.V_NaturalGas[time]/T_NG .* (  (1 +  err_w_CH4[time] + err_V_NG[time])  .*(prop_O2_CH4 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[1][time]/M_CH4)
                                                        + (1 +  err_w_C2H6[time] + err_V_NG[time]) .*(prop_O2_C2H6 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[2][time]/M_C2H6)      
                                                        + (1 +  err_w_C3H8[time] + err_V_NG[time]) .*(prop_O2_C3H8 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[3][time]/M_C3H8)
                                                       )
                                                    .==
                measurements.V_Air[time]* 0.21 * M_Air/T_Air .* (  (1/M_CH4) * measurements.wi_NaturalGas[1][time] .*(1 + err_w_CH4[time] + err_V_Air[time])
                                                                +  (1/M_C2H6) * measurements.wi_NaturalGas[2][time] .*(1 + err_w_C2H6[time] + err_V_Air[time])
                                                                +  (1/M_C3H8) * measurements.wi_NaturalGas[3][time] .*(1 + err_w_C3H8[time] + err_V_Air[time])
                                                                )
                                                                .*(0.21 * M_O2 + 0.79 * M_N2)

)













# @constraint(m, measurements.V_NaturalGas[time]/T_NG .* ( (1 + err_w_CO2[time]) .*( prop_H2O_CH4 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[1][time]/M_CO2/M_CH4
#                                                                                 + prop_H2O_C2H6 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[2][time]/M_CO2/M_C2H6      
#                                                                                 + prop_H2O_C3H8 * measurements.wi_Fumes[1][time] .* measurements.wi_NaturalGas[3][time]/M_CO2/M_C3H8
#                                                                                 )
#                                                         + (1 + err_w_H2O[time]) .* ( prop_H2O_CH4 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[1][time]/M_H2O/M_CH4
#                                                                                     + prop_H2O_C2H6 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[2][time]/M_H2O/M_C2H6      
#                                                                                     + prop_H2O_C3H8 * measurements.wi_Fumes[2][time] .* measurements.wi_NaturalGas[3][time]/M_H2O/M_C3H8
#                                                                                 )
#                                                         + (1 + err_w_N2[time]) .* ( prop_H2O_CH4 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[1][time]/M_N2/M_CH4
#                                                                                 + prop_H2O_C2H6 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[2][time]/M_N2/M_C2H6      
#                                                                                 + prop_H2O_C3H8 * measurements.wi_Fumes[3][time] .* measurements.wi_NaturalGas[3][time]/M_N2/M_C3H8
#                                                                                 )
#                                                     )
#                                                     .==
#                 measurements.V_HotFumes[time]/M_H2O/T_HotFumes .*((1/M_CH4) * measurements.wi_NaturalGas[1][time] .*(1 + err_w_CH4[time] + err_V_HotFumes[time] + err_w_H2O[time])
#                                                                 +(1/M_C2H6) * measurements.wi_NaturalGas[2][time] .*(1 + err_w_C2H6[time] + err_V_HotFumes[time] + err_w_H2O[time])
#                                                                 +(1/M_C3H8) * measurements.wi_NaturalGas[2][time] .*(1 + err_w_C3H8[time] + err_V_HotFumes[time] + err_w_H2O[time])
#                                                                 )
# )


#Constraint de O2

#Linearisation
@constraint(m, -err_V_NG_bound[time] .<= measurements.V_NaturalGas[time] - V_NG[time])
@constraint(m, measurements.V_NaturalGas[time] - V_NG[time]  .<= err_V_NG_bound[time])

@constraint(m, -err_V_Air_bound[time] .<= measurements.V_Air[time] - V_Air[time])
@constraint(m, measurements.V_Air[time] - V_Air[time] .<= err_V_Air_bound[time])

@constraint(m, -err_V_Hot_bound[time] .<= measurements.V_HotFumes[time] - V_HotFumes[time])
@constraint(m, measurements.V_HotFumes[time] - V_HotFumes[time] .<= err_V_Hot_bound[time])

@constraint(m, -err_w_CH4_bound[time] .<= err_w_CH4[time])
@constraint(m, err_w_CH4[time] .<= err_w_CH4_bound[time])
@constraint(m, -err_w_C2H6_bound[time] .<= err_w_C2H6[time])
@constraint(m, err_w_C2H6[time] .<= err_w_C2H6_bound[time])
@constraint(m, -err_w_C3H8_bound[time] .<= err_w_C3H8[time])
@constraint(m, err_w_C3H8[time] .<= err_w_C3H8_bound[time])

@constraint(m, -err_w_CO2_bound[time] .<= err_w_CO2[time])
@constraint(m,  err_w_CO2[time] .<= err_w_CO2_bound[time])
@constraint(m, -err_w_H2O_bound[time] .<= err_w_H2O[time])
@constraint(m, err_w_H2O[time] .<= err_w_H2O_bound[time])

println("The optimization problem to be solved is:")
print(m)

solve(m)

println("Objective value: ", getobjectivevalue(m))
println("===================================================================================")
println(getvalue(err_V_NG_bound))
println("===================================================================================")
println(getvalue(err_V_Air_bound))
println("===================================================================================")
println(getvalue(err_V_Hot_bound))
println("===================================================================================")
println(getvalue(err_w_CH4_bound))
println("===================================================================================")
println(getvalue(err_w_C2H6_bound))
println("===================================================================================")
println(getvalue(err_w_C3H8_bound))
println("===================================================================================")
println(getvalue(err_w_CO2_bound))
println("===================================================================================")
println(getvalue(err_w_H2O_bound))
println("===================================================================================")
println(getvalue(V_NG))
println("===================================================================================")
println(getvalue(V_Air))
println("===================================================================================")
println(getvalue(V_HotFumes))

                
                  
                
               
