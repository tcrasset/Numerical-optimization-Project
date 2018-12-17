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


#------------------------ CONSTANTS -----------------------------
# Molar masses
M_CH4 = 16.04246
M_O2 = 31.9988
M_N2 = 28.0134
M_Air = 28.850334
M_CO2 = 44.0095
M_H2O = 18.01528

# Temperature
T_NG = 25 + 273.15
T_Air = 25 + 273.15
T_HotFumes = 1600 + 273.15

# Reaction coefficients
coeff_CH4 = 1
coeff_O2 = 2
coeff_N2 = 2 * (79/21)
coeff_Air = coeff_O2 + coeff_N2
coeff_CO2 = 1
coeff_H2O = 2

#------------------------ MODEL -----------------------------------
measurements = loadDataFromFile("q1")
n_Obs = length(measurements.V_NaturalGas)
time = 1:n_Obs
m = Model(solver=ClpSolver())

@variable(m, err_Air_bound[time] >= 0.0)
@variable(m, err_Hot_bound[time] >= 0.0)
@variable(m, err_NG_bound[time] >= 0.0)

@variable(m,  V_NG[time] >= 0.0)
@variable(m,  V_Air[time] >= 0.0)
@variable(m,  V_HotFumes[time] >= 0.0)
@variable(m,  V_O2[time] >= 0.0)

@objective(m, Min, sum(err_NG_bound) + sum(err_Air_bound) + sum(err_Hot_bound))

M_Fumes_Inv = measurements.wi_Fumes[1][time]/M_CO2 + measurements.wi_Fumes[2][time]/M_H2O + measurements.wi_Fumes[3][time]/M_N2

@constraint(m, V_NG[time]/T_NG .== (coeff_CH4/coeff_CO2) * measurements.wi_Fumes[1][time] ./(M_CO2 * M_Fumes_Inv[time]) .* V_HotFumes[time]/T_HotFumes)
@constraint(m, V_NG[time]/T_NG .== (coeff_CH4/coeff_H2O) * measurements.wi_Fumes[2][time] ./(M_H2O * M_Fumes_Inv[time]) .* V_HotFumes[time]/T_HotFumes)
@constraint(m, V_NG[time]/T_NG .== (coeff_CH4/coeff_O2) * 21/100 * (V_Air[time]/T_Air))

#Linearisation of the absolute difference constraints
@constraint(m, -err_NG_bound[time] .<= measurements.V_NaturalGas[time] - V_NG[time] )
@constraint(m, measurements.V_NaturalGas[time] - V_NG[time]  .<= err_NG_bound[time])
@constraint(m, -err_Air_bound[time] .<= measurements.V_Air[time] - V_Air[time])
@constraint(m, measurements.V_Air[time] - V_Air[time] .<= err_Air_bound[time])
@constraint(m, -err_Hot_bound[time] .<= measurements.V_HotFumes[time] - V_HotFumes[time])
@constraint(m, measurements.V_HotFumes[time] - V_HotFumes[time] .<= err_Hot_bound[time])

status = solve(m)
if(status == :Optimal)

    println("Objective value: ", getobjectivevalue(m))
    println(getvalue(err_NG_bound))
    println(getvalue(err_Air_bound))
    println(getvalue(err_Hot_bound))
    println(getvalue(V_NG))
    println(getvalue(V_Air))
    println(getvalue(V_HotFumes))

    figure()
    plot(time, measurements.V_NaturalGas, linestyle=":",linewidth=2, label="Data")
    plot(time, [ getvalue(V_NG[t]) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Natural Gas volume flow")
    legend()
    savefig("q1_V_NG.svg")

    figure()
    plot(time, measurements.V_Air, linestyle=":",linewidth=2, label="Data")
    plot(time, [ getvalue(V_Air[t]) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Air Volume flow")
    legend()
    savefig("q1_V_Air.svg")
    
    figure()
    plot(time, measurements.V_HotFumes, linestyle=":",linewidth=2, label="Data")
    plot(time, [ getvalue(V_HotFumes[t]) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Hot Fumes volume flow")
    legend()
    savefig("q1_V_HF.svg")
end
