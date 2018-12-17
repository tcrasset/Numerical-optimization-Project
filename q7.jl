include("./data_struct.jl")
using Gurobi
using Juniper
using Ipopt
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

#Polynomial coefficients
A_CO2 = 58.16639
B_CO2 = 2.720074
C_CO2 = -0.492289
D_CO2 = 0.038844
E_CO2 = -6.447293
F_CO2 = -425.9186
G_CO2 = 263.6125
H_CO2 = -393.5224
CO2_poly_coeff = [A_CO2, B_CO2, C_CO2, D_CO2, E_CO2, F_CO2, G_CO2, H_CO2]

A_H2O = 41.96426
B_H2O = 8.62205
C_H2O = -1.499780
D_H2O = 0.098119
E_H2O = -11.15764
F_H2O = -272.1797
G_H2O = 219.7809
H_H2O = -241.8264
H2O_poly_coeff  = [A_H2O, B_H2O, C_H2O, D_H2O, E_H2O, F_H2O, G_H2O, H_H2O]

A_N2 = 19.50583
B_N2 = 19.88705
C_N2 = -8.598535
D_N2 = 1.369784
E_N2 = 0.527601
F_N2 = -4.935202
G_N2 = 212.3900
H_N2 = 0.0
N2_poly_coeff  = [A_N2, B_N2, C_N2, D_N2, E_N2, F_N2, G_N2, H_N2]

CH4_heat_comb = 802.34
C2H6_heat_comb = 1437.2
C3H8_heat_comb = 2044.2
#------------------------ END OF CONSTANTS --------------------------------------------------------
function heat_absorbed(temp, gaz)
    """
    Temperature in K
    Gaz = 1,2 or 3 for [ 'CO2', 'H2O, 'N2' ] respectively
    """
    temp = temp/1000
    if gaz == 1
        coeff = CO2_poly_coeff
        # return A_N2*temp +B_N2/2 * temp^2 +C_N2/3 * temp * temp * temp
        #        D_N2/4 * temp * temp * temp * temp - E_N2/temp + G_N2 -H_N2
    elseif gaz == 2
        coeff = H2O_poly_coeff
    elseif gaz == 3
        coeff = N2_poly_coeff
    end
    return coeff[1]*temp + coeff[2]/2 * temp^2 + coeff[3]/3 * temp * temp * temp
            coeff[4]/4 * temp * temp * temp * temp - coeff[5]/temp + coeff[6] - coeff[7]
end
# #------------------------ MODEL ---------------------------------------------------------------------

measurements = loadDataFromFile("q3_easy")

m = Model(solver = IpoptSolver(max_iter = 5000))
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
@variable(m, err_w_N2_bound[time] >= 0.0)

@variable(m, err_w_CH4[time])
@variable(m, err_w_C2H6[time])
@variable(m, err_w_C3H8[time])

@variable(m, err_w_CO2[time])
@variable(m, err_w_H2O[time])
@variable(m, err_w_N2[time])

@variable(m,  err_V_NG[time])
@variable(m,  err_V_Hot[time])
@variable(m,  err_V_Air[time])

@variable(m, T_HotFumes[time] >= 0.0, start = 1600 + 273.15)

@objective(m, Min, sum(err_V_NG_bound) + sum(err_V_Air_bound) + sum(err_V_Hot_bound)
                + sum(err_w_CH4_bound) + sum(err_w_C2H6_bound) + sum(err_w_C3H8_bound)
                + sum(err_w_CO2_bound) + sum(err_w_H2O_bound) + sum(err_w_N2_bound))

#Constraint on CO2

@NLconstraint(m, [t in time], (measurements.wi_Fumes[1][t]/M_CO2) * (measurements.V_HotFumes[t]/T_HotFumes[t]) * ( (measurements.wi_NaturalGas[1][t]/M_CH4) * (1 + err_V_Hot[t] + err_w_CH4[t] + err_w_CO2[t])
                                                                                                                + (measurements.wi_NaturalGas[2][t]/M_C2H6) * (1 + err_V_Hot[t] + err_w_C2H6[t] + err_w_CO2[t])
                                                                                                                + (measurements.wi_NaturalGas[3][t]/M_C3H8) * (1 + err_V_Hot[t] + err_w_C3H8[t] + err_w_CO2[t])
                                                                                                        )
                                                                                                        ==
                            (measurements.V_NaturalGas[t]/T_NG)   * ( (prop_CO2_CH4 * measurements.wi_NaturalGas[1][t]/M_CH4) * (   (measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                            + (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                            + (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                          ) 
                                                                    + (prop_CO2_C2H6 * measurements.wi_NaturalGas[2][t]/M_C2H6) * ((measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                            +  (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                            +  (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                        ) 
                                                                    + (prop_CO2_C3H8 * measurements.wi_NaturalGas[3][t]/M_C3H8) * ((measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                            +  (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                            +  (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                        )
                                                                    )  
)                                 
         
#Constraint on H2O

@NLconstraint(m, [t in time], (measurements.wi_Fumes[2][t]/M_H2O) * (measurements.V_HotFumes[t]/T_HotFumes[t]) * ( (measurements.wi_NaturalGas[1][t]/M_CH4) * (1 + err_V_Hot[t] + err_w_CH4[t] + err_w_H2O[t])
                                                                                                                + (measurements.wi_NaturalGas[2][t]/M_C2H6) * (1 + err_V_Hot[t] + err_w_C2H6[t] + err_w_H2O[t])
                                                                                                                + (measurements.wi_NaturalGas[3][t]/M_C3H8) * (1 + err_V_Hot[t] + err_w_C3H8[t] + err_w_H2O[t])
                                                                                                                )
                                                                                                        ==
                            (measurements.V_NaturalGas[t]/T_NG)   * ( (prop_H2O_CH4 * measurements.wi_NaturalGas[1][t]/M_CH4) * ( (measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                                + (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                                + (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                          ) 
                                                                    + (prop_H2O_C2H6 * measurements.wi_NaturalGas[2][t]/M_C2H6) * ((measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                            +  (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                            +  (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                        ) 
                                                                    + (prop_H2O_C3H8 * measurements.wi_NaturalGas[3][t]/M_C3H8) * ((measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                            +  (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                            +  (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                        )
                                                            )  
)        

#Constraint on O2
@NLconstraint(m, [t in time], 0.21 * measurements.V_Air[t]/T_Air * ( measurements.wi_NaturalGas[1][t]/M_CH4 * (1 +  err_w_CH4[t] + err_V_Air[t] )
                                                                    + measurements.wi_NaturalGas[2][t]/M_C2H6 * (1 +  err_w_C2H6[t] + err_V_Air[t] )
                                                                    + measurements.wi_NaturalGas[3][t]/M_C3H8 * (1 +  err_w_C3H8[t] + err_V_Air[t] )
                                                                    )
                                                        ==
                                    measurements.V_NaturalGas[t]/T_NG * (  prop_O2_CH4 * measurements.wi_NaturalGas[1][t]/M_CH4 * (1 +  err_w_CH4[t] + err_V_NG[t]) 
                                                                        + prop_O2_C2H6 * measurements.wi_NaturalGas[2][t]/M_C2H6 * (1 +  err_w_C2H6[t] + err_V_NG[t]) 
                                                                        + prop_O2_C3H8 * measurements.wi_NaturalGas[3][t]/M_C3H8 * (1 +  err_w_C3H8[t] + err_V_NG[t]) 
                                                                        )
)

#Registering user defined function
JuMP.register(m, :heat_absorbed, 2, heat_absorbed, autodiff=true)

#Temperature constraint
@NLconstraint(m, [t in time], measurements.V_NaturalGas[t]/T_NG *( measurements.wi_NaturalGas[1][t]/M_CH4 * CH4_heat_comb*((measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                        +  (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                        +  (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_CH4[t] + err_V_NG[t]))
                                                                                                                    ) 
                                                            +  measurements.wi_NaturalGas[2][t]/M_C2H6 * C2H6_heat_comb* ((measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                        + (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                        + (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_C2H6[t] + err_V_NG[t]))
                                                                                                                    ) 
                                                            +  measurements.wi_NaturalGas[3][t]/M_C3H8 * C3H8_heat_comb* ((measurements.wi_Fumes[1][t]/M_CO2 * (1 + err_w_CO2[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                        + (measurements.wi_Fumes[2][t]/M_H2O * (1 + err_w_H2O[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                        + (measurements.wi_Fumes[3][t]/M_N2 * (1 + err_w_N2[t] + err_w_C3H8[t] + err_V_NG[t]))
                                                                                                                    )
                                                    )
                                                    ==
            measurements.V_HotFumes[t]/T_HotFumes[t] * ( measurements.wi_Fumes[1][t]/M_CO2 * heat_absorbed(T_HotFumes[t],1) *((measurements.wi_NaturalGas[1][t]/M_CH4) * (1 + err_V_Hot[t] + err_w_CH4[t] + err_w_CO2[t])
                                                                                                                            + (measurements.wi_NaturalGas[2][t]/M_C2H6) * (1 + err_V_Hot[t] + err_w_C2H6[t] + err_w_CO2[t])
                                                                                                                            + (measurements.wi_NaturalGas[3][t]/M_C3H8) * (1 + err_V_Hot[t] + err_w_C3H8[t] + err_w_CO2[t])
                                                                                                                             )
                                                        + measurements.wi_Fumes[2][t]/M_H2O * heat_absorbed(T_HotFumes[t],2) * (  (measurements.wi_NaturalGas[1][t]/M_CH4) * (1 + err_V_Hot[t] + err_w_CH4[t] + err_w_H2O[t])
                                                                                                                                + (measurements.wi_NaturalGas[2][t]/M_C2H6) * (1 + err_V_Hot[t] + err_w_C2H6[t] + err_w_H2O[t])
                                                                                                                                + (measurements.wi_NaturalGas[3][t]/M_C3H8) * (1 + err_V_Hot[t] + err_w_C3H8[t] + err_w_H2O[t])
                                                                                                                                )
                                                        + measurements.wi_Fumes[3][t]/M_N2 * heat_absorbed(T_HotFumes[t],3) * (   (measurements.wi_NaturalGas[1][t]/M_CH4) * (1 + err_V_Hot[t] + err_w_CH4[t] + err_w_N2[t])
                                                                                                                                + (measurements.wi_NaturalGas[2][t]/M_C2H6) * (1 + err_V_Hot[t] + err_w_C2H6[t] + err_w_N2[t])
                                                                                                                                + (measurements.wi_NaturalGas[3][t]/M_C3H8) * (1 + err_V_Hot[t] + err_w_C3H8[t] + err_w_N2[t])
                                                                                                                                )
                                                        )
)


@constraint(m, [t in time], measurements.wi_NaturalGas[1][t] * (1 + err_w_CH4[t])
             + measurements.wi_NaturalGas[2][t] * (1 + err_w_C2H6[t]) 
             + measurements.wi_NaturalGas[3][t] * (1 + err_w_C3H8[t]) == 1)

@constraint(m,[t in time], measurements.wi_Fumes[1][t] * (1 + err_w_CO2[t])
             + measurements.wi_Fumes[2][t] * (1 + err_w_H2O[t]) 
             + measurements.wi_Fumes[3][t] * (1 + err_w_N2[t]) == 1)


errVNG_constr = @NLconstraint(m,[t in time],  err_V_NG[t]^2 <=  err_V_NG_bound[t])
errVAir_constr = @NLconstraint(m,[t in time], err_V_Air[t]^2<=  err_V_Air_bound[t])
errVHot_constr = @NLconstraint(m,[t in time],  err_V_Hot[t]^2<=  err_V_Hot_bound[t])

errwCH4_constr = @NLconstraint(m,[t in time], err_w_CH4[t]^2 <= err_w_CH4_bound[t])
errwC2J6_constr = @NLconstraint(m,[t in time], err_w_C2H6[t]^2<=  err_w_C2H6_bound[t])
errwC3H8_constr = @NLconstraint(m,[t in time], err_w_C3H8[t]^2<=  err_w_C3H8_bound[t])

errwH2O_constr = @NLconstraint(m,[t in time], err_w_H2O[t]^2<=  err_w_H2O_bound[t])
errwCO2_constr = @NLconstraint(m,[t in time], err_w_CO2[t]^2 <= err_w_CO2_bound[t])
errwN2_constr = @NLconstraint(m, [t in time], err_w_N2[t]^2 <= err_w_N2_bound[t])

# println("The optimization problem to be solved is:")
# print(m)

status = solve(m)
if(status == :Optimal)

    figure()
    plot(time, measurements.V_NaturalGas, linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.V_NaturalGas[t] * (1 +getvalue(err_V_NG[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Natural Gas volume flow")
    legend()

    figure()
    plot(time, measurements.V_Air, linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.V_Air[t] * (1 +getvalue(err_V_Air[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Air Volume flow")
    legend()

    figure()
    plot(time, measurements.V_HotFumes, linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.V_HotFumes[t] * (1 +getvalue(err_V_Hot[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("Hot Fumes volume flow")
    legend()

    figure()
    plot(time, measurements.wi_NaturalGas[1], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_NaturalGas[1][t] * (1 +getvalue(err_w_CH4[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("CH4 mass percent")
    legend()

    figure()
    plot(time, measurements.wi_NaturalGas[2], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_NaturalGas[2][t] * (1 +getvalue(err_w_C2H6[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("C2H6 mass percent")
    legend()

    figure()
    plot(time, measurements.wi_NaturalGas[3], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_NaturalGas[3][t] * (1 +getvalue(err_w_C3H8[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("C3H8 mass percent")
    legend()

    figure()
    plot(time, measurements.wi_Fumes[1], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_Fumes[1][t] * (1 +getvalue(err_w_CO2[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("CO2 mass percent")
    legend()   
    
    
    figure()
    plot(time, measurements.wi_Fumes[2], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_Fumes[2][t] * (1 +getvalue(err_w_H2O[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("H2O mass percent")
    legend()

    figure()
    plot(time, measurements.wi_Fumes[3], linestyle=":",linewidth=2, label="Data")
    plot(time, [ measurements.wi_Fumes[3][t] * (1 +getvalue(err_w_N2[t])) for t = time ], linestyle="-",linewidth=2, label="Clean Data")
    xlabel("Time period")
    ylabel("N2 mass percent")
    legend()
# end

                  
                
               



