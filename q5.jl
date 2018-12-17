include("./data_struct.jl")
include("./sensitivity_gurobi.jl")
using Gurobi
using Data
using Ipopt
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

#------------------------ CONSTANTS ----------------------------------
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

#------------------------ MODEL -----------------------------------
measurements = loadDataFromFile("q3")

m = Model(solver=GurobiSolver())
n_Obs = length(measurements.V_NaturalGas)
time = 1:n_Obs

@variable(m, err_w_CH4_bound[time] >= 0.0)
@variable(m, err_w_C2H6_bound[time] >= 0.0)
@variable(m, err_w_C3H8_bound[time] >= 0.0)
@variable(m, err_w_CO2_bound[time] >= 0.0)
@variable(m, err_w_H2O_bound[time] >= 0.0)
@variable(m, err_w_N2_bound[time] >= 0.0)
@variable(m,  err_V_NG_bound[time] >= 0.0)
@variable(m,  err_V_HotFumes_bound[time] >= 0.0)
@variable(m,  err_V_Air_bound[time] >= 0.0)

@variable(m, err_w_CH4[time])
@variable(m, err_w_C2H6[time])
@variable(m, err_w_C3H8[time])
@variable(m, err_w_CO2[time])
@variable(m, err_w_H2O[time])
@variable(m, err_w_N2[time])
@variable(m,  err_V_NG[time])
@variable(m,  err_V_HotFumes[time])
@variable(m,  err_V_Air[time])


@objective(m, Min, sum(err_V_NG_bound[t] + err_V_Air_bound[t] + err_V_HotFumes_bound[t]
                + err_w_CH4_bound[t] + err_w_C2H6_bound[t] + err_w_C3H8_bound[t]
                + err_w_CO2_bound[t] + err_w_H2O_bound[t] + err_w_N2_bound[t] for t=time))

#Constraint on CO2
CO2_constr = @constraint(m, (measurements.wi_Fumes[1][time]/M_CO2) .* (measurements.V_HotFumes[time]/T_HotFumes) .* ( (measurements.wi_NaturalGas[1][time]/M_CH4) .* (1 + err_V_HotFumes[time] + err_w_CH4[time] + err_w_CO2[time])
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
H2O_constr = @constraint(m, (measurements.wi_Fumes[2][time]/M_H2O) .* (measurements.V_HotFumes[time]/T_HotFumes) .* ( (measurements.wi_NaturalGas[1][time]/M_CH4) .* (1 + err_V_HotFumes[time] + err_w_CH4[time] + err_w_H2O[time])
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
O2_constr = @constraint(m, 0.21 * measurements.V_Air[time]/T_Air .* ( measurements.wi_NaturalGas[1][time]/M_CH4 .* (1 +  err_w_CH4[time] + err_V_Air[time] )
                                                        + measurements.wi_NaturalGas[2][time]/M_C2H6 .* (1 +  err_w_C2H6[time] + err_V_Air[time] )
                                                        + measurements.wi_NaturalGas[3][time]/M_C3H8 .* (1 +  err_w_C3H8[time] + err_V_Air[time] )
                                                        )
                                                        .==
                measurements.V_NaturalGas[time]/T_NG .* (  prop_O2_CH4 * measurements.wi_NaturalGas[1][time]/M_CH4 .* (1 +  err_w_CH4[time] + err_V_NG[time]) 
                                                        + prop_O2_C2H6 * measurements.wi_NaturalGas[2][time]/M_C2H6 .* (1 +  err_w_C2H6[time] + err_V_NG[time]) 
                                                        + prop_O2_C3H8 * measurements.wi_NaturalGas[3][time]/M_C3H8 .* (1 +  err_w_C3H8[time] + err_V_NG[time]) 
                                                        )

)

# Constraints on the sum of mass fractions
wNG_constr = @constraint(m, measurements.wi_NaturalGas[1][time] .* (1 + err_w_CH4[time])
             + measurements.wi_NaturalGas[2][time] .* (1 + err_w_C2H6[time]) 
             + measurements.wi_NaturalGas[3][time] .* (1 + err_w_C3H8[time]) .== 1)

wFumes_constr = @constraint(m, measurements.wi_Fumes[1][time] .* (1 + err_w_CO2[time])
             + measurements.wi_Fumes[2][time] .* (1 + err_w_H2O[time]) 
             + measurements.wi_Fumes[3][time] .* (1 + err_w_N2[time]) .== 1)


# Constraints on the bounds
errVNG_constr = @constraint(m,[t in time],  norm(err_V_NG[t]) <=  err_V_NG_bound[t])
errVAir_constr = @constraint(m,[t in time], norm(err_V_Air[t])<=  err_V_Air_bound[t])
errVHot_constr = @constraint(m,[t in time],  norm(err_V_HotFumes[t])<=  err_V_HotFumes_bound[t])

errwCH4_constr = @constraint(m,[t in time], norm(err_w_CH4[t]) <= err_w_CH4_bound[t])
errwC2J6_constr = @constraint(m,[t in time], norm(err_w_C2H6[t])<=  err_w_C2H6_bound[t])
errwC3H8_constr = @constraint(m,[t in time], norm(err_w_C3H8[t])<=  err_w_C3H8_bound[t])

errwH2O_constr = @constraint(m,[t in time], norm(err_w_H2O[t])<=  err_w_H2O_bound[t])
errwCO2_constr = @constraint(m,[t in time], norm(err_w_CO2[t]) <= err_w_CO2_bound[t])
errwN2_constr = @constraint(m, [t in time], norm(err_w_N2[t]) <= err_w_N2_bound[t])


status = solve(m)
if(status == :Optimal)
    figure()
    plot(time, measurements.V_NaturalGas, linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.V_NaturalGas[t] * (1 +getvalue(err_V_NG[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Natural Gas volume flow")
    legend()
    savefig("q5_V_NG.svg")

    figure()
    plot(time, measurements.V_Air, linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.V_Air[t] * (1 +getvalue(err_V_Air[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Air Volume flow")
    legend()
    savefig("q5_V_Air.svg")

    figure()
    plot(time, measurements.V_HotFumes, linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.V_HotFumes[t] * (1 +getvalue(err_V_Hot[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Hot Fumes volume flow")
    legend()
    savefig("q5_V_HF.svg")

    figure()
    plot(time, measurements.wi_NaturalGas[1], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_NaturalGas[1][t] * (1 +getvalue(err_w_CH4[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("CH4 mass percent")
    legend()
    savefig("q5_w_CH4.svg")

    figure()
    plot(time, measurements.wi_NaturalGas[2], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_NaturalGas[2][t] * (1 +getvalue(err_w_C2H6[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("C2H6 mass percent")
    legend()
    savefig("q5_w_C2H6.svg")

    figure()
    plot(time, measurements.wi_NaturalGas[3], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_NaturalGas[3][t] * (1 +getvalue(err_w_C3H8[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("C3H8 mass percent")
    legend()
    savefig("q5_w_C3H8.svg")


    figure()
    plot(time, measurements.wi_Fumes[1], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_Fumes[1][t] * (1 +getvalue(err_w_CO2[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("CO2 mass percent")
    legend()   
    savefig("q5_w_CO2.svg")
    
    
    figure()
    plot(time, measurements.wi_Fumes[2], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_Fumes[2][t] * (1 +getvalue(err_w_H2O[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("H2O mass percent")
    legend()
    savefig("q5_w_H2O.svg")

    figure()
    plot(time, measurements.wi_Fumes[3], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_Fumes[3][t] * (1 +getvalue(err_w_N2[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("N2 mass percent")
    legend()
    savefig("q5_w_N2.svg")
end
               
