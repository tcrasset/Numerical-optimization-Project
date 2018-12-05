include("./data_struct.jl")
using Data
using JuMP

Pkg.add("Ipopt")
using Ipopt
# Pkg.add("CPLEX")
# using CPLEX

###################################################

# struct Measurements
#     V_NaturalGas::Array{Float64}
#     V_Air::Array{Float64}
#     V_HotFumes::Array{Float64}
#     wi_NaturalGas::Array{Array{Float64}} # order : CH4, C2H6, C3H8
#     wi_Fumes::Array{Array{Float64}} # order : CO2, H2O, N2
# end


#------------------------ CONSTANTS --------------------------------------------------------

#Masses molaires
M_CH4 = 16.04246
M_02 = 31.9988
M_N2 = 28.0134
M_Air = 28.850334
M_CO2 = 44.0095
M_H20 = 18.01528
#Coefficients réactions
coeff_CH4 = 1
coeff_02 = 2
coeff_N2 = 2 * (79/21)
coeff_Air = coeff_02 + coeff_N2
coeff_CO2 = 1
coeff_H20 = 2


#Nb moles totales
n_reactifs_tots = coeff_CH4 + coeff_Air
n_produits_tots = coeff_CO2 + coeff_H20 + coeff_N2

#Fractions molaires réactifs
x_CH4 = coeff_CH4/n_reactifs_tots
x_Air = coeff_Air/ n_reactifs_tots
x_O2 = coeff_02/ n_reactifs_tots
x_N2 = coeff_N2/n_reactifs_tots

#Fractions molaires produits
x_C02 = coeff_C02/n_produits_tots
x_H20 = coeff_H2O/ n_produits_tots
x_N2 = coeff_N2/n_produits_tots

#Masse molaire moyenne
M_reactifs = x_CH4 * M_CH4 + coeff_Air * M_Air
M_produits = x_C02 * M_CO2 + x_H20 * M_H20 + x_N2 * M_N2




#------------------------ END OF CONSTANTS --------------------------------------------------------


#------------------------FUNCTIONS --------------------------------------------------------------

function V_to_N_gaslaw(P::Float64, T::Float64,V::Float64)
    """Params 
        P - pressure in atm
        T - temperature in Celsius
        V - volume in m3
    """
    return (P * 1000 * V)/(0.08205 * (273.15 + T))
end


#------------------------ MODEL ---------------------------------------------------------------------
measurements = loadDataFromFile("q1_easy")

print(measurements.V_NaturalGas)
print(measurements.V_Air)
print(measurements.V_HotFumes)
print(measurements.wi_NaturalGas)
print(measurements.wi_Fumes)

length(measurements.V_Air)

#Solve for wi_Air




m = Model(solver=IpoptSolver())
n_Obs = length(measurements.V_Air)

for t = 1:n_Obs

    #The measure of V is imprecise
    n_CH4_measured = V_to_N_gaslaw(P_CH4,T_CH4,measurements.V_NaturalGas[t])
    n_Air_measured = V_to_N_gaslaw(P_Air,T_Air,measurements.V_Air[t])

    n_O2_measured = n_Air_measured -  #Comment je calcules n_O2 a partir de N_Air ?
    n_N2_measured #???

    #Quand est ce que je divise par n_reactifs_tots pour avoir la fraction massique ?

    #OU est ce que je mets l'erreur de mesure ? Dans V ? Dans n_measured ? Dans n_consumed ?

    #Assume that there is a little bit more that was consumed than what was measured
    n_CH4_consumed = coeff_CH4 * (n_CH4_measured + err_CH4)  # 1 mol CH4
    n_O2_consumed = n_O2_measured - coeff_02 * (n_CH4_measured + err_CH4) # 2 mol O2
    n_N2_consumed = n_N2_measured - coeff_N2 * (n_CH4_measured + err_CH4) # 2 * 79/21 mol N2


    #j'ai commencé à écrire les choses en fractions massiques puis j'ai changé en haut donc ce qui suit 
    # est faux d'un facteur n_produits_tots
    
    #Verify that what is measured exactly (mass percentage) is actually what is coming out
    @constraint(measurements.wi_Fumes[1][t] == coeff_C02 * (1/M_produits) * x_consumed * M_CH4)     #1 mol CO2
    @constraint(measurements.wi_Fumes[2][t] == coeff_H20 * (1/M_produits) * x_consumed * M_H20) # 2 mol H20
    @constraint(measurements.wi_Fumes[3][t] == coeff_N2 * (1/M_produits) * x_consumed * M_N2) # 2 * 79/21 mol N2


    
    solve(m)
end