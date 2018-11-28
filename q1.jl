include("./data_struct.jl")
Pkg.add("Ipopt")
using Data
using Ipopt
using JuMP

measurements = loadDataFromFile("q1_easy")

print(measurements.V_Air)

print(measurements.wi_NaturalGas)
length(measurements.V_Air)
print(measurements.V_Air.len)
w_i_sum_Fumes = [sum(measurements.wi_Fumes[1],1),sum(measurements.wi_Fumes[2],1),sum(measurements.wi_Fumes[3],1)]


print(sum_w_i_Fumes[1])
phi_i = x_i # Ideal gas law

M_CH4 = 16.043
M_O = 15.9994
M_N = 14.0067
M_Air = 28.9647
M_CO2 = 44.01
M_H20 = 18.01528

M = 

for i = 1:numCtr
    M += 
end

w_Air = x_i * M_i / (sum x_i M_i)
w_CH4 = x_i * M_i / (sum x_i M_i)
w_CH4 = x_i * M_i / (sum x_i M_i)
w_CH4 = x_i * M_i / (sum x_i M_i)



print(w_i_sum_Fumes[1][1])
#Solve for wi_Air



for t = 1:length(measurements.V_Air[1])

	m = Model(solver=IpoptSolver())

	@variable(m, wi_Air >= 0)
	
    @constraint(m, measurements.wi_NaturalGas[1][t] - wi_Air - w_i_sum_Fumes[t][1] == 0)
        
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

    # struct Measurements
    #     V_NaturalGas::Array{Float64}
    #     V_Air::Array{Float64}
    #     V_HotFumes::Array{Float64}
    #     wi_NaturalGas::Array{Array{Float64}} # order : CH4, C2H6, C3H8
    #     wi_Fumes::Array{Array{Float64}} # order : CO2, H2O, N2
    # end