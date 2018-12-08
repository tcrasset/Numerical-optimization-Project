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
measurements = loadDataFromFile("q2")

m = Model(solver=ClpSolver())
n_Obs = length(measurements.V_NaturalGas)
time = 1:n_Obs
@variable(m, err_V_Air_bound[time] >= 0.0)
@variable(m, err_V_Hot_bound[time] >= 0.0)
@variable(m, err_V_NG_bound[time] >= 0.0)

@variable(m, err_w_CH4_bound[time] >= 0.0)
@variable(m, err_w_C2H6_bound[time] >= 0.0)
@variable(m, err_w_C3H8_bound[time] >= 0.0)

@variable(m, err_w_CO2_bound[time] >= 0.0)
@variable(m, err_w_H2O_bound[time] >= 0.0)
@variable(m, err_w_O2_bound[time] >= 0.0)


@variable(m,  V_NG[time] >= 0.0)
@variable(m,  V_HotFumes[time] >= 0.0)
@variable(m,  V_Air[time] >= 0.0)

@variable(m, w_CH4[time] >= 0.0)
@variable(m, w_C2H6[time] >= 0.0)
@variable(m, w_C3H8[time] >= 0.0)

@variable(m, w_CO2[time] >= 0.0)
@variable(m, w_H2O[time] >= 0.0)
@variable(m, w_O2[time] >= 0.0)




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

#     wi_NaturalGas::Array{Array{Float64}} # order : CH4, C2H6, C3H8
#     wi_Fumes::Array{Array{Float64}} # order : CO2, H2O, N2

#Constraint de CO2

@constraint(m,  (
                 (prop_CO2_CH4/M_CH4) *  (w_CH4[time] * w_CO2[time] + w_CH4[time] * err_w_CO2_bound[time] + err_w_CH4_bound[time] * w_CO2)
                +(prop_CO2_C2H6/M_C2H6)* (w_C2H6[time] * w_CO2[time] + w_C2H6[time] * err_w_CO2_bound[time] + err_w_C2H6_bound[time] * w_CO2)
                +(prop_CO2_C3H8/M_C3H8)* (w_C3H8[time] * w_CO2[time] + w_C3H8[time] * err_w_CO2_bound[time] + err_w_C3H8_bound[time] * w_CO2)
                )
            *   V_NG/T_NG 
            ==  V_HotFumes/T_HotFumes 
            *   (   
                 (1/(M_CH4 * M_CO2)) * (w_CH4[time] * w_CO2[time] + w_CH4[time] * err_w_CO2_bound[time] + err_w_CH4_bound[time] * w_CO2)
                +(1/(M_C2H6 * M_CO2))* (w_C2H6[time] * w_CO2[time] + w_C2H6[time] * err_w_CO2_bound[time] + err_w_C2H6_bound[time] * w_CO2)
                +(1/(M_C3H8 * M_CO2))* (w_C3H8[time] * w_CO2[time] + w_C3H8[time] * err_w_CO2_bound[time] + err_w_C3H8_bound[time] * w_CO2)
                +(1/(M_CH4 * M_H2O)) * (w_CH4[time] * w_H2O[time] + w_CH4[time] * err_w_H2O_bound[time] + err_w_CH4_bound[time] * w_H2O)
                +(1/(M_C2H6 * M_H2O))* (w_C2H6[time] * w_H2O[time] + w_C2H6[time] * err_w_H2O_bound[time] + err_w_C2H6_bound[time] * w_H2O)
                +(1/(M_C3H8 * M_H2O))* (w_C3H8[time] * w_H2O[time] + w_C3H8[time] * err_w_H2O_bound[time] + err_w_C3H8_bound[time] * w_H2O)
                +(1/(M_CH4 * M_N2))  * (w_CH4[time] * w_N2[time] + w_CH4[time] * err_w_N2_bound[time] + err_w_CH4_bound[time] * w_N2)
                +(1/(M_C2H6 * M_N2)) * (w_C2H6[time] * w_N2[time] + w_C2H6[time] * err_w_N2_bound[time] + err_w_C2H6_bound[time] * w_N2)
                +(1/(M_C3H8 * M_N2)) * (w_C3H8[time] * w_N2[time] + w_C3H8[time] * err_w_N2_bound[time] + err_w_C3H8_bound[time] * w_N2)
                )
            *   1/M_CO2
)

#Constraint de H2O

@constraint(m,  (
                 (prop_H2O_CH4/M_CH4) *  (w_CH4[time] * w_H2O[time] + w_CH4[time] * err_w_H2O_bound[time] + err_w_CH4_bound[time] * w_H2O)
                +(prop_H2O_C2H6/M_C2H6)* (w_C2H6[time] * w_H2O[time] + w_C2H6[time] * err_w_H2O_bound[time] + err_w_C2H6_bound[time] * w_H2O)
                +(prop_H2O_C3H8/M_C3H8)* (w_C3H8[time] * w_H2O[time] + w_C3H8[time] * err_w_H2O_bound[time] + err_w_C3H8_bound[time] * w_H2O)
                )
            *   V_NG/T_NG 
            ==  V_HotFumes/T_HotFumes 
            *   (   
                 (1/(M_CH4 * M_CO2)) * (w_CH4[time] * w_CO2[time] + w_CH4[time] * err_w_CO2_bound[time] + err_w_CH4_bound[time] * w_CO2)
                +(1/(M_C2H6 * M_CO2))* (w_C2H6[time] * w_CO2[time] + w_C2H6[time] * err_w_CO2_bound[time] + err_w_C2H6_bound[time] * w_CO2)
                +(1/(M_C3H8 * M_CO2))* (w_C3H8[time] * w_CO2[time] + w_C3H8[time] * err_w_CO2_bound[time] + err_w_C3H8_bound[time] * w_CO2)
                +(1/(M_CH4 * M_H2O)) * (w_CH4[time] * w_H2O[time] + w_CH4[time] * err_w_H2O_bound[time] + err_w_CH4_bound[time] * w_H2O)
                +(1/(M_C2H6 * M_H2O))* (w_C2H6[time] * w_H2O[time] + w_C2H6[time] * err_w_H2O_bound[time] + err_w_C2H6_bound[time] * w_H2O)
                +(1/(M_C3H8 * M_H2O))* (w_C3H8[time] * w_H2O[time] + w_C3H8[time] * err_w_H2O_bound[time] + err_w_C3H8_bound[time] * w_H2O)
                +(1/(M_CH4 * M_N2))  * (w_CH4[time] * w_N2[time] + w_CH4[time] * err_w_N2_bound[time] + err_w_CH4_bound[time] * w_N2)
                +(1/(M_C2H6 * M_N2)) * (w_C2H6[time] * w_N2[time] + w_C2H6[time] * err_w_N2_bound[time] + err_w_C2H6_bound[time] * w_N2)
                +(1/(M_C3H8 * M_N2)) * (w_C3H8[time] * w_N2[time] + w_C3H8[time] * err_w_N2_bound[time] + err_w_C3H8_bound[time] * w_N2)
                )
            *   1/M_H2O
)


                
                  
                
               
