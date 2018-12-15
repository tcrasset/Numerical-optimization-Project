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

#------------------------ CONSTANTS --------------------------------------------------------
# Molar masses
M_CH4 = 16.04246
M_C2H6 = 30.06904
M_C3H8 = 44.09562
M_O2 = 31.9988
M_N2 = 28.0134
M_Air = 28.850334
M_CO2 = 44.0095
M_H2O = 18.01528
M_NG = M_CH4 + M_C2H6 + M_C3H8
# Temperature
T_NG = 25 + 273.15
T_Air = 25 + 273.15
T_HotFumes = 1600 + 273.15

#------------------------ END OF CONSTANTS --------------------------------------------------------

#Fumes proportions for 1 MOLE of NG
prop_CO2_CH4 = 1
prop_H2O_CH4 = 2    
prop_O2_CH4 = 2
prop_CO2_C2H6 = 2
prop_H2O_C2H6 = 3
prop_O2_C2H6 = 3.5
prop_CO2_C3H8 = 3 
prop_H2O_C3H8 = 4
prop_O2_C3H8 = 5

#------------------------ MODEL ---------------------------------------------------------------------
measurements = loadDataFromFile("q3_easy")

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



#Constraint on CO2

@constraint(m, (measurements.wi_Fumes[1][time]/M_CO2) .* (measurements.V_HotFumes[time]/T_HotFumes) .* ( (measurements.wi_NaturalGas[1][time]/M_CH4) .* (1 + err_V_HotFumes[time] + err_w_CH4[time] + err_w_CO2[time])
                                                                                                        + (measurements.wi_NaturalGas[2][time]/M_C2H6) .* (1 + err_V_HotFumes[time] + err_w_C2H6[time] + err_w_CO2[time])
                                                                                                        + (measurements.wi_NaturalGas[3][time]/M_C3H8) .* (1 + err_V_HotFumes[time] + err_w_C3H8[time] + err_w_CO2[time])
                                                                                                        )
                                                                                                        .==
                (measurements.V_NaturalGas[time]/T_NG)   .* ( (prop_CO2_CH4 * measurements.wi_NaturalGas[1][time]/M_CH4) .* (   (measurements.wi_Fumes[1][time]/M_CO2 .* (1 + err_w_CO2[time] + err_w_CH4[time] + err_V_NG[time]))
                                                                                                                            + (measurements.wi_Fumes[2][time]/M_H2O .* (1 + err_w_H2O[time] + err_w_CH4[time] + err_V_NG[time]))
                                                                                                                            + (measurements.wi_Fumes[3][time]/M_N2 .* (1 + err_w_N2[time] + err_w_CH4[time] + err_V_NG[time]))
                                                                                                                          ) 
                                                            + (prop_CO2_C2H6 * measurements.wi_NaturalGas[2][time]/M_C2H6) .* ((measurements.wi_Fumes[1][time]/M_CO2 .* (1 + err_w_CO2[time] + err_w_C2H6[time] + err_V_NG[time]))
                                                                                                                            +  (measurements.wi_Fumes[2][time]/M_H2O .* (1 + err_w_H2O[time] + err_w_C2H6[time] + err_V_NG[time]))
                                                                                                                            +  (measurements.wi_Fumes[3][time]/M_N2 .* (1 + err_w_N2[time] + err_w_C2H6[time] + err_V_NG[time]))
                                                                                                                        ) 
                                                            + (prop_CO2_C3H8 * measurements.wi_NaturalGas[3][time]/M_C3H8) .* ((measurements.wi_Fumes[1][time]/M_CO2 .* (1 + err_w_CO2[time] + err_w_C3H8[time] + err_V_NG[time]))
                                                                                                                            +  (measurements.wi_Fumes[2][time]/M_H2O .* (1 + err_w_H2O[time] + err_w_C3H8[time] + err_V_NG[time]))
                                                                                                                            +  (measurements.wi_Fumes[3][time]/M_N2 .* (1 + err_w_N2[time] + err_w_C3H8[time] + err_V_NG[time]))
                                                                                                                        )
                                                            )  
)                                 
         
#Constraint on H2O

@constraint(m, (measurements.wi_Fumes[2][time]/M_H2O) .* (measurements.V_HotFumes[time]/T_HotFumes) .* ( (measurements.wi_NaturalGas[1][time]/M_CH4) .* (1 + err_V_HotFumes[time] + err_w_CH4[time] + err_w_H2O[time])
                                                                                                        + (measurements.wi_NaturalGas[2][time]/M_C2H6) .* (1 + err_V_HotFumes[time] + err_w_C2H6[time] + err_w_H2O[time])
                                                                                                        + (measurements.wi_NaturalGas[3][time]/M_C3H8) .* (1 + err_V_HotFumes[time] + err_w_C3H8[time] + err_w_H2O[time])
                                                                                                        )
                                                                                                        .==
                (measurements.V_NaturalGas[time]/T_NG)   .* ( (prop_H2O_CH4 * measurements.wi_NaturalGas[1][time]/M_CH4) .* (   (measurements.wi_Fumes[1][time]/M_CO2 .* (1 + err_w_CO2[time] + err_w_CH4[time] + err_V_NG[time]))
                                                                                                                            + (measurements.wi_Fumes[2][time]/M_H2O .* (1 + err_w_H2O[time] + err_w_CH4[time] + err_V_NG[time]))
                                                                                                                            + (measurements.wi_Fumes[3][time]/M_N2 .* (1 + err_w_N2[time] + err_w_CH4[time] + err_V_NG[time]))
                                                                                                                          ) 
                                                            + (prop_H2O_C2H6 * measurements.wi_NaturalGas[2][time]/M_C2H6) .* ((measurements.wi_Fumes[1][time]/M_CO2 .* (1 + err_w_CO2[time] + err_w_C2H6[time] + err_V_NG[time]))
                                                                                                                            +  (measurements.wi_Fumes[2][time]/M_H2O .* (1 + err_w_H2O[time] + err_w_C2H6[time] + err_V_NG[time]))
                                                                                                                            +  (measurements.wi_Fumes[3][time]/M_N2 .* (1 + err_w_N2[time] + err_w_C2H6[time] + err_V_NG[time]))
                                                                                                                        ) 
                                                            + (prop_H2O_C3H8 * measurements.wi_NaturalGas[3][time]/M_C3H8) .* ((measurements.wi_Fumes[1][time]/M_CO2 .* (1 + err_w_CO2[time] + err_w_C3H8[time] + err_V_NG[time]))
                                                                                                                            +  (measurements.wi_Fumes[2][time]/M_H2O .* (1 + err_w_H2O[time] + err_w_C3H8[time] + err_V_NG[time]))
                                                                                                                            +  (measurements.wi_Fumes[3][time]/M_N2 .* (1 + err_w_N2[time] + err_w_C3H8[time] + err_V_NG[time]))
                                                                                                                        )
                                                            )  
)        

#Constraint on O2
@constraint(m, 0.21 * measurements.V_Air[time]/T_Air .* ( measurements.wi_NaturalGas[1][time]/M_CH4 .* (1 +  err_w_CH4[time] + err_V_Air[time] )
                                                        + measurements.wi_NaturalGas[2][time]/M_C2H6 .* (1 +  err_w_C2H6[time] + err_V_Air[time] )
                                                        + measurements.wi_NaturalGas[3][time]/M_C3H8 .* (1 +  err_w_C3H8[time] + err_V_Air[time] ))
                                                        .==
                measurements.V_NaturalGas[time]/T_NG .* (  prop_O2_CH4 * measurements.wi_NaturalGas[1][time]/M_CH4 .* (1 +  err_w_CH4[time] + err_V_NG[time]) 
                                                        + prop_O2_C2H6 * measurements.wi_NaturalGas[2][time]/M_C2H6 .* (1 +  err_w_C2H6[time] + err_V_NG[time]) 
                                                        + prop_O2_C3H8 * measurements.wi_NaturalGas[3][time]/M_C3H8 .* (1 +  err_w_C3H8[time] + err_V_NG[time]) 
                                                        )

) 

@constraint(m, measurements.wi_NaturalGas[1][time] .* (1 + err_w_CH4[time])
             + measurements.wi_NaturalGas[2][time] .* (1 + err_w_C2H6[time]) 
             + measurements.wi_NaturalGas[3][time] .* (1 + err_w_C3H8[time]) .== 1)

@constraint(m, measurements.wi_Fumes[1][time] .* (1 + err_w_CO2[time])
             + measurements.wi_Fumes[2][time] .* (1 + err_w_H2O[time]) 
             + measurements.wi_Fumes[3][time] .* (1 + err_w_N2[time]) .== 1)


# Linearisation
@constraint(m, -err_V_NG_bound[time] .<= err_V_NG[time] )
@constraint(m,  err_V_NG[time]  .<= err_V_NG_bound[time])
@constraint(m, -err_V_Air_bound[time] .<=  err_V_Air[time])
@constraint(m, err_V_Air[time] .<= err_V_Air_bound[time])
@constraint(m, -err_V_Hot_bound[time] .<= err_V_HotFumes[time])
@constraint(m,  err_V_HotFumes[time] .<= err_V_Hot_bound[time])

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

status = solve(m)
if(status == :Optimal)
    println("Objective value: ", getobjectivevalue(m))
    println("===================================================================================")
    println(getvalue(err_V_NG_bound))
    println("===================================================================================")
    print("V_NG:")
    println(measurements.V_NaturalGas[time])
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

                  
                
               
