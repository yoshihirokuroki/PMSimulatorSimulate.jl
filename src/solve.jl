using DifferentialEquations

# Function to help check for duplicate IDs
hasduplicates(xs) = !allunique(xs)



# OOP Solves
function _solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing}, evs::PMSimulatorBase.PMEvents; kwargs...)
    sol_out = Dict{Union{Symbol, Int64}, PMParameterizedSolve.PMSolution}()
    IDs = [instance.ID for instance in evs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in events") : nothing
    for instance in evs.instances
        mdl_i = deepcopy(mdl) # Create copy to prevent modification of OG model
        evi = vcat(instance.inputs, instance.updates)
        cbs = collect_evs(evi, mdl_i)
        sol_i = PMParameterizedSolve.solve(mdl_i, alg; callback = cbs, kwargs...)
        println(instance.ID)
        println(typeof(sol_i))
        sol_out[instance.ID] = sol_i
    end
    if length(sol_out) == 1
        return Base.values(sol_out)[1]
    else
        return sol_out
    end
end
    

function _solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing}, evs::Vector{PMEvent}; kwargs...)
    mdl_i = deepcopy(mdl)
    cbs = collect_evs(evs, mdl_i)
    sol = PMParameterizedSolve.solve(mdl_i, alg; callback = cbs, kwargs...)
    return sol
end

function _solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing}, data::PMSimulatorBase.DataFrames.AbstractDataFrame; kwargs...)
    sol_out = Dict{Union{Symbol, Int64}, PMParameterizedSolve.PMSolution}()
    dfevs = PMSimulatorBase.df2evs(data)
    IDs = [instance.ID for instance in dfevs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in dataframe") : nothing
    for instance in dfevs.instances
        mdl_i = deepcopy(mdl) # Create copy to prevent modification of OG model
        evi = vcat(instance.inputs, instance.updates)
        cbs = collect_evs(evi, mdl_i)
        sol_i = PMParameterizedSolve.solve(mdl_i, alg; callback = cbs, kwargs...)
        # push!(sol_out, sol_i)
        sol_out[instance.ID] = sol_i
    end
    if length(sol_out) == 1
        return Base.values(sol_out)[1]
    else
        return sol_out
    end
end

function DifferentialEquations.solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; kwargs...)
    if :evs ∈ keys(kwargs) && :data ∈ keys(kwargs)
        error("Cannot define both evs and data kwargs")
    elseif :evs ∈ keys(kwargs)
        kwin = Dict(key => kwargs[key] for key in keys(kwargs) if key != :evs)
        out = _solve(mdl, alg, kwargs[:evs]; kwin...)
    elseif :data ∈ keys(kwargs)
        kwin = Dict(key => kwargs[key] for key in keys(kwargs) if key != :data)
        out = _solve(mdl, alg, kwargs[:data]; kwin...)
    else
        out = PMParameterizedSolve.solve(mdl, alg; kwargs...)
    end
    return out
end




    
# IIP Solves
function _solve!(mdl::PMModel, alg::Union{DEAlgorithm,Nothing}, evs::PMSimulatorBase.PMEvents; kwargs...)
    mdls_out = Dict{Symbol, PMParameterizedBase.PMModel}()
    IDs = [instance.ID for instance in evs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in events") : nothing
    for instance in evs.instances
        mdl_i = deepcopy(mdl)
        evi = vcat(instance.inputs, instance.updates)
        cbs = collect_evs(evi, mdl_i)
        PMParameterizedSolve.solve!(mdl_i, alg; callback = cbs, kwargs...)
        # push!(mdls_out, mdl_i)
        mdls_out[instance.ID, mdl_i]
    end
    if length(mdls_out) == 1
        return Base.values(mdls_out)[1]
    else
        return mdls_out
    end
end

function _solve!(mdl::PMModel, alg::Union{DEAlgorithm,Nothing}, evs::Vector{PMEvent}; kwargs...)
    mdl_out = deepcopy(mdl)
    cbs = collect_evs(evs, mdl_out)
    PMParameterizedSolve.solve!(mdl_out, alg; callback = cbs, kwargs...)
    return mdl_out
end

function _solve!(mdl::PMModel, alg::Union{DEAlgorithm,Nothing}, data::PMSimulatorBase.DataFrames.AbstractDataFrame; kwargs...)
    mdls_out = Dict{Union{Symbol, Int64}, PMParameterizedBase.PMModel}()
    dfevs = PMSimulatorBase.df2evs(data)
    IDs = [instance.id for instance in dfevs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in dataframe") : nothing
    for instance in dfevs.instances
        mdl_i = deepcopy(mdl) # Create copy to prevent modification of OG model
        evi = vcat(instance.inputs, instance.update)
        cbs = collect_evs(evi, mdl_i)
        PMParameterizedSolve.solve!(mdl_i, alg; callback = cbs, kwargs...)
        mdls_out[instance.ID] = mdl_i
    end
    if length(mdls_out) == 1
        return Base.values(mdls_out)[1]
    else
        return mdls_out
    end
end
    

function DifferentialEquations.solve!(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; kwargs...)
    if :evs ∈ keys(kwargs) && :data ∈ keys(kwargs)
        error("Cannot define both evs and data kwargs")
    elseif :evs ∈ keys(kwargs)
        kwin = Dict(key => kwargs[key] for key in keys(kwargs) if key != :evs)
        out = _solve!(mdl, alg, kwargs[:evs]; kwin...)
        return out
    elseif :data ∈ keys(kwargs)
        kwin = Dict(key => kwargs[key] for key in keys(kwargs) if key != :data)
        out = _solve!(mdl, alg, kwargs[:data]; kwin...)
        return out
    else
        PMParameterizedSolve.solve!(mdl, alg; kwargs...)
        return mdl
    end
end


