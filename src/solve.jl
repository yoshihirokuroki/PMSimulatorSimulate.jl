using DifferentialEquations

# OOP Solves
function solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; evs::PMSimulatorBase.PMEvents = PMSimulatorBase.PMEvents(PMSimulatorBase.PMInstance[]), kwargs...)
    sol_out = PMParameterizedSolve.PMSolution[]
    for instance in evs.instances
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl)
        sol_i = PMParameterizedSolve.solve(mdl, alg; callback = cbs, kwargs)
        push!(sol_out, sol_i)
    end
    if length(sol_out) == 1
        return sol_out[1]
    else
        return sol_out
    end
end
    

function solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; data::UnionPMSimulatorBase.DataFrames.AbstractDataFrame = DataFrame(), kwargs...)
    sol_out = PMParameterizedSolve.PMSolution[]
    dfevs = PMSimulatorBase.df2evs(data)
    for instance in dfevs.instances
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl)
        sol_i = PMParameterizedSolve.solve(mdl, alg; callback = cbs, kwargs)
        push!(sol_out, sol_i)
    end
    if length(sol_out) == 1
        return sol_out[1]
    else
        return sol_out
    end
end
    
# IIP Solves
function solve!(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; evs::PMSimulatorBase.PMEvents = PMSimulatorBase.PMEvents(PMSimulatorBase.PMInstance[]) , kwargs...)
    mdls_out = PMParameterizedBase.PMModel[]
    for instance in evs.instances
        mdl_i = deepcopy(mdl)
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl)
        PMParameterizedSolve.solve!(mdl_i, alg; callback = cbs, kwargs)
        push!(mdls_out, mdl_i)
    end
    if length(mdls_out) == 1
        return mdls_out[1]
    else
        return mdls_out
    end
end
    

function solve!(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; data::UnionPMSimulatorBase.DataFrames.AbstractDataFrame = DataFrame(), kwargs...)
    mdls_out = PMParameterizedBase.PMModel[]
    dfevs = PMSimulatorBase.df2evs(data)
    for instance in dfevs.instances
        mdl_i = deepcopy(mdl)
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl)
        PMParameterizedSolve.solve!(mdl_i, alg; callback = cbs, kwargs)
        push!(mdls_out, mdl_i)
    end
    if length(mdls_out) == 1
        return mdls_out[1]
    else
        return mdls_out
    end
end