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
M_CH4 = 16.043
M_O = 15.9994
M_N = 14.0067
M_Air = 28.9647
M_CO2 = 44.01
M_H20 = 18.01528

#Nb moles totales
X = 1 + 2 + 2*79/21

#Nb moles reactifs
n_CH4 = 1
n_O2 = 2
n_N2 = 2 * 79/21

#Fractions molaires réactifs
x_CH4 = n_CH4/ X
x_Air = (n_O2 + n_N2)/X 
x_O2 = n_O2/X 
x_N2 = n_N2/X

#Test == 1 :
x_CH4 + x_Air == 1.0
x_CH4 + x_O2 + x_N2 == 1.0

#Masse molaire moyenne
M = x_CH4* M_CH4 + x_O2*2*M_O + x_N2*2*M_N 

#Fractions massiques
w_O2 = x_O2 * 2*M_O / M
w_CH4 = x_CH4 * M_CH4 / M
w_N2 = x_N2 * 2*M_N / M

#Test == 1 :
w_O2 + w_CH4 + w_N2 == 1.0

#Loi des gaz parfaits
pV = nRT
V = nRT #p == 1 atm
# V = n * 8.2057*10.0^(-5) * (celsius + 273.15) 

V_O2 =  x_O2 * 8.2057*10.0^(-5) * (25 + 273.15) 
V_CH4 =  x_CH4 * 8.2057*10.0^(-5) * (25 + 273.15) 
V_N2 =  x_N2 * 8.2057*10.0^(-5) * (25 + 273.15) 

V_O2+ V_CH4 + V_N2 == V #Not true bc of floating point arithmetic



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
w_i_sum_Fumes = [sum(measurements.wi_Fumes[1],1),sum(measurements.wi_Fumes[2],1),sum(measurements.wi_Fumes[3],1)]

print(sum_w_i_Fumes[1])

print(w_i_sum_Fumes[1][1])
#Solve for wi_Air




m = Model(solver=IpoptSolver())
# m = Model(solver=CplexSolver())

@variable(m, wi_Air >= 0)
for t = 1:length(measurements.V_Air[1])
    
    #Doesnt add up 
    # @constraint(m, measurements.V_NaturalGas + measurements.V_Air - measurements.V_HotFumes == 0)
    
    #Doesn"t work for iteration 3
    #@constraint(m, measurements.wi_NaturalGas[1][t] - wi_Air - w_i_sum_Fumes[t][1] == 0)
    
    solve(m)
    println(getvalue(wi_Air))
end

#t = 1
m = Model(solver=IpoptSolver())
@variable(m, wi_Air >= 0)
@constraint(m, measurements.wi_NaturalGas[1][1] - wi_Air - w_i_sum_Fumes[1][1] == 0)
solve(m)
println(getvalue(wi_Air))
#t = 2
m = Model(solver=IpoptSolver())
@variable(m, wi_Air >= 0)
@constraint(m, measurements.wi_NaturalGas[1][2] - wi_Air - w_i_sum_Fumes[2][1] == 0)
solve(m)
println(getvalue(wi_Air))
#t = 3
m = Model(solver=IpoptSolver())
@variable(m, wi_Air >= 0)
@constraint(m, measurements.wi_NaturalGas[1][3] - wi_Air - w_i_sum_Fumes[3][1] == 0)
solve(m)
println(getvalue(wi_Air))
#t = 4
m = Model(solver=IpoptSolver())
@variable(m, wi_Air >= 0)
@constraint(m, measurements.wi_NaturalGas[1][4] - wi_Air - w_i_sum_Fumes[4][1] == 0)
solve(m)
println(getvalue(wi_Air))

