using DifferentialEquations

# Function to help check for duplicate IDs
hasduplicates(xs) = !allunique(xs)


# OOP Solves
function solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; evs::PMSimulatorBase.PMEvents = PMSimulatorBase.PMEvents(PMSimulatorBase.PMInstance[]), kwargs...)
    sol_out = Dict{Union{Symbol, Int64}, PMParameterized.PMSolution}
    IDs = [instance.id for instance in evs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in events") : nothing
    for instance in evs.instances
        mdl_i = deepcopy(mdl) # Create copy to prevent modification of OG model
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl_i)
        sol_i = PMParameterizedSolve.solve(mdl_i, alg; callback = cbs, kwargs)
        sol_out[instance.ID] = sol_i
    end
    if length(sol_out) == 1
        return Base.values(sol_out)[1]
    else
        return sol_out
    end
end
    

function solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; evs::Vector{PMEvent} = PMSimulatorBase.PMEvent[], kwargs...)
    mdl_i = deepcopy(mdl)
    cbs = collect_evs(evs, mdl_i)
    sol = PMParameterizedSolve.solve(mdl_i, alg; callback = cbs, kwargs)
    return sol
end

function solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; data::UnionPMSimulatorBase.DataFrames.AbstractDataFrame = DataFrame(), kwargs...)
    sol_out = Dict{Union{Symbol, Int64}, PMParameterized.PMSolution}
    dfevs = PMSimulatorBase.df2evs(data)
    IDs = [instance.id for instance in dfevs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in dataframe") : nothing
    for instance in dfevs.instances
        mdl_i = deepcopy(mdl) # Create copy to prevent modification of OG model
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl_i)
        sol_i = PMParameterizedSolve.solve(mdl_i, alg; callback = cbs, kwargs)
        # push!(sol_out, sol_i)
        sol_out[instance.ID] = sol_i
    end
    if length(sol_out) == 1
        return Base.values(sol_out)[1]
    else
        return sol_out
    end
end
    
# IIP Solves
function solve!(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; evs::PMSimulatorBase.PMEvents = PMSimulatorBase.PMEvents(PMSimulatorBase.PMInstance[]) , kwargs...)
    mdls_out = Dict{Symbol, PMParameterizedBase.PMModel}
    IDs = [instance.id for instance in evs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in events") : nothing
    for instance in evs.instances
        mdl_i = deepcopy(mdl)
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl_i)
        PMParameterizedSolve.solve!(mdl_i, alg; callback = cbs, kwargs)
        # push!(mdls_out, mdl_i)
        mdls_out[instance.ID, mdl_i]
    end
    if length(mdls_out) == 1
        return Base.values(mdls_out)[1]
    else
        return mdls_out
    end
end

function solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; data::UnionPMSimulatorBase.DataFrames.AbstractDataFrame = DataFrame(), kwargs...)
    mdls_out = Dict{Union{Symbol, Int64}, PMParameterized.PMModel}
    dfevs = PMSimulatorBase.df2evs(data)
    IDs = [instance.id for instance in dfevs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in dataframe") : nothing
    for instance in dfevs.instances
        mdl_i = deepcopy(mdl) # Create copy to prevent modification of OG model
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl_i)
        PMParameterizedSolve.solve!(mdl_i, alg; callback = cbs, kwargs)
        mdls_out[instance.ID] = mdl_i
    end
    if length(mdls_out) == 1
        return Base.values(mdls_out)[1]
    else
        return mdls_out
    end
end
    


function solve!(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; evs::Vector{PMEvent} = PMSimulatorBase.PMEvent[], kwargs...)
    mdl_out = deepcopy(mdl)
    cbs = collect_evs(evs, mdl_out)
    PMParameterizedSolve.solve!(mdl_out, alg; callback = cbs, kwargs)
    return mdl_out
end

