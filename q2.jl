include("./data_struct.jl")
using Clp
using Data
using JuMP
using PyCall
using PyPlot


###################################################

# struct Measurements
#     V_NaturalGas::Array{Float64}
#     V_Air::Array{Float64}
#     V_HotFumes::Array{Float64}
#     wi_NaturalGas::Array{Array{Float64}} # order : CH4, C2H6, C3H8
#     wi_Fumes::Array{Array{Float64}} # order : CO2, H2O, N2
# end
###################################################

#------------------------ CONSTANTS -----------------------
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

#------------------------ MODEL ---------------------------------------------------------------------
measurements = loadDataFromFile("q2")

m = Model(solver=ClpSolver())
n_Obs = length(measurements.V_NaturalGas)
time = 1:n_Obs
@variable(m, err_Air_bound[time] >= 0.0)
@variable(m, err_Hot_bound[time] >= 0.0)
@variable(m, err_NG_bound[time] >= 0.0)

@variable(m,  V_NG[time] >= 0.0)
@variable(m,  V_HotFumes[time] >= 0.0)
@variable(m,  V_Air[time] >= 0.0)

@objective(m, Min, sum(err_NG_bound) + sum(err_Air_bound) + sum(err_Hot_bound))

#Transformer les volumes d'air et de hotfumes en leur composants
#N2 ne participe pas à la réaction
M_Fumes_Inv = measurements.wi_Fumes[1][time]/M_CO2 + measurements.wi_Fumes[2][time]/M_H2O + measurements.wi_Fumes[3][time]/M_N2
M_NG_Inv = measurements.wi_NaturalGas[1][time]/M_CH4 + measurements.wi_NaturalGas[2][time]/M_C2H6 + measurements.wi_NaturalGas[3][time]/M_C3H8
M_Air_Inv = 1/(0.21 * M_O2 + 0.79 * M_N2)

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

# Contrainte avec CO2
@constraint(m, V_NG[time]/T_NG .==  (M_NG_Inv./M_Fumes_Inv)   .* (measurements.wi_Fumes[1][time]/M_CO2)
                                                                ./  ( prop_CO2_CH4 * measurements.wi_NaturalGas[1][time] / M_CH4
                                                                    + prop_CO2_C2H6 * measurements.wi_NaturalGas[2][time] / M_C2H6
                                                                    + prop_CO2_C3H8 * measurements.wi_NaturalGas[3][time] / M_C3H8
                                                                    )
                                                                .* V_HotFumes[time]/T_HotFumes)
# Contrainte avec H2O
@constraint(m, V_NG[time]/T_NG .==  (M_NG_Inv./M_Fumes_Inv)   .* (measurements.wi_Fumes[2][time]/M_H2O)
                                                                ./  ( prop_H2O_CH4 * measurements.wi_NaturalGas[1][time] / M_CH4
                                                                    + prop_H2O_C2H6 * measurements.wi_NaturalGas[2][time] / M_C2H6
                                                                    + prop_H2O_C3H8 * measurements.wi_NaturalGas[3][time] / M_C3H8
                                                                    )
                                                                .* V_HotFumes[time]/T_HotFumes)
# Contrainte avec O2 
@constraint(m, V_NG[time]/T_NG .== M_NG_Inv    * 0.21
                                                ./  ( prop_O2_CH4 * measurements.wi_NaturalGas[1][time] / M_CH4
                                                    + prop_O2_C2H6 * measurements.wi_NaturalGas[2][time] / M_C2H6
                                                    + prop_O2_C3H8 * measurements.wi_NaturalGas[3][time] / M_C3H8
                                                    )
                                                .* V_Air[time]/T_Air)
#Linearisation

@constraint(m, -err_NG_bound[time] .<= measurements.V_NaturalGas[time] - V_NG[time])
@constraint(m, measurements.V_NaturalGas[time] - V_NG[time]  .<= err_NG_bound[time])

@constraint(m, -err_Air_bound[time] .<= measurements.V_Air[time] - V_Air[time])
@constraint(m, measurements.V_Air[time] - V_Air[time] .<= err_Air_bound[time])

@constraint(m, -err_Hot_bound[time] .<= measurements.V_HotFumes[time] - V_HotFumes[time])
@constraint(m, measurements.V_HotFumes[time] - V_HotFumes[time] .<= err_Hot_bound[time])


println("The optimization problem to be solved is:")
print(m)

status = solve(m)


if(status == :Optimal)
    
    println("Objective value: ", getobjectivevalue(m))
    println("===================================================================================")
    println(getvalue(err_NG_bound))
    println("===================================================================================")
    println(getvalue(err_Air_bound))
    println("===================================================================================")
    println(getvalue(err_Hot_bound))
    println("===================================================================================")
    println(getvalue(V_NG))
    println("===================================================================================")
    println(getvalue(V_Air))
    println("===================================================================================")
    println(getvalue(V_HotFumes))
        
    figure()
    suptitle("Natural Gas", fontsize=12)
    plot(time, measurements.V_NaturalGas, linestyle=":",linewidth=2, label="Data")
    plot(time, [ getvalue(V_NG[t]) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Natural Gas volume flow")
    legend()


    figure()
    suptitle("Air", fontsize=12)
    plot(time, measurements.V_Air, linestyle=":",linewidth=2, label="Data")
    plot(time, [ getvalue(V_Air[t]) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Air Volume flow")
    legend()

    figure()
    suptitle("Hot Fumes", fontsize=12)
    plot(time, measurements.V_HotFumes, linestyle=":",linewidth=2, label="Data")
    plot(time, [ getvalue(V_HotFumes[t]) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Hot Fumes volume flow")
    legend()
end