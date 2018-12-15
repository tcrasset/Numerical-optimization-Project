include("./data_struct.jl")
using Clp
using Data
using JuMP
<<<<<<< HEAD
using PyCall
using PyPlot


=======
using Clp
using PyCall
using PyPlot
>>>>>>> 40be3613c7eea66728386a250ca4fd89c27ac5f9

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
<<<<<<< HEAD
    println(getvalue(err_NG_bound))
    println(getvalue(err_Air_bound))
    println(getvalue(err_Hot_bound))
    println(getvalue(V_NG))
    println(getvalue(V_Air))
    println(getvalue(V_HotFumes))

    print("------------ VOLUME NATURAL GAZ -------------------")
    print(V_NG)
    print("---------------- VOLUME AIR -----------------------")
    print(V_Air)
    print("------------- VOLUME HOT FUMES --------------------")
    print(V_HotFumes)

    figure()
    suptitle("Natural Gas", fontsize=12)
    plot(time, measurements.V_NaturalGas, linestyle=":",linewidth=2, label="Data")
    plot(time, [ getvalue(V_NG[t]) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
=======
    println("===================================================================================")
    println("Volume of Natural gases: ", getvalue(V_NG))
    println("===================================================================================")
    println("Volume of Air: ", getvalue(V_Air))
    println("===================================================================================")
    println("Volume of Hot Fumes: ", getvalue(V_HotFumes))
    println("===================================================================================")
    println("Error on Natural Gases volume: ", getvalue(err_NG_bound))
    println("===================================================================================")
    println("Error on Air volume: ", getvalue(err_Air_bound))
    println("===================================================================================")
    println("Error on Hot Fumes volume: ", getvalue(err_Hot_bound))
    
    figure()
    suptitle("Natural Gas", fontsize=12)
    plot(time, measurements.V_NaturalGas, "r", label="Data")
    plot(time, [getvalue(err_NG_bound[t]) for t = time], "b", label="Clean Data")
>>>>>>> 40be3613c7eea66728386a250ca4fd89c27ac5f9
    xlabel("Time period")
    ylabel("Natural Gas volume flow")
    legend()

<<<<<<< HEAD

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
=======
    figure()
    suptitle("Hot Fumes", fontsize=12)
    plot(time, measurements.V_HotFumes, "r", label="Data")
    plot(time, [getvalue(err_Hot_bound[t]) for t = time], "b", label="Clean Data")
    xlabel("Time period")
    ylabel("Hot Fumes volume flow")
    legend()

    figure()
    suptitle("Air", fontsize=12)
    plot(time, measurements.V_Air, "r", label="Data")
    plot(time, [getvalue(err_Air_bound[t]) for t = time], "b", label="Clean Data")
    xlabel("Time period")
    ylabel("Air Volume flow")
>>>>>>> 40be3613c7eea66728386a250ca4fd89c27ac5f9
    legend()
end
