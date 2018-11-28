using JLD

###################################################

module Data

struct Measurements
	V_NaturalGas::Array{Float64}
	V_Air::Array{Float64}
	V_HotFumes::Array{Float64}
	wi_NaturalGas::Array{Array{Float64}} # order : CH4, C2H6, C3H8
	wi_Fumes::Array{Array{Float64}} # order : CO2, H2O, N2
end

end

###################################################

function loadDataFromFile(instanceName::String)
	fileName = string(instanceName, ".jld")
	return load(fileName, "instance")

	# fileName = string(instanceName, ".bin")
	# in = open(fileName,"r")
	# data = deserialize(in)
	# close(in)
	# return data
end

function saveDataToFile(data::Data.Measurements, instanceName::String)
	fileName = string(instanceName, ".jld")
	save(fileName, "instance", data)

	# fileName = string(instanceName, ".bin")
	# out = open(fileName,"w")
	# serialize(out, data)
	# close(out)

	println("Data saved to $fileName")
end
;
