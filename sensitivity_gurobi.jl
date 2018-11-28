using JuMP
using Gurobi

const LinConstrRef = ConstraintRef{Model,LinearConstraint}

function getBasisBoundsGurobi(c::LinConstrRef)
	if length(c.m.linconstrDuals) != MathProgBase.numlinconstr(c.m)
		error("Sensitivity information not available.")
	end
	solverModel = MathProgBase.getrawsolver(c.m)
	beginIdx = c.idx
	numElem = 1
	lower = Gurobi.get_dblattrarray(solverModel, "SARHSLow", beginIdx, numElem)
	upper = Gurobi.get_dblattrarray(solverModel, "SARHSUp", beginIdx, numElem)
	return lower[1], upper[1]
end

function getBasisBoundsGurobi{DIM}(c::Array{JuMP.ConstraintRef, DIM})
	if length(c[1].m.linconstrDuals) != MathProgBase.numlinconstr(c[1].m)
		error("Sensitivity information not available.")
	end
	numCtr = length(c)
	lowers = Array{Float64}(numCtr)
	uppers = Array{Float64}(numCtr)
	for i = 1:numCtr
		# Gotta call the C API for each ctr and not for all of them at once
		# because the index may not be in order or contiguous
		lowers[i], uppers[i] = getBasisBoundsGurobi(c[i])
	end
	return lowers, uppers
end

function getBasisBoundsGurobi{DIM, ANY}(c::JuMP.JuMPArray{JuMP.ConstraintRef, DIM, ANY})
	return getBasisBoundsGurobi(c.innerArray)
end
