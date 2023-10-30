using DifferentialEquations

# Function to help check for duplicate IDs
hasduplicates(xs) = !allunique(xs)


# OOP Solves
function PMParameterizedSolve.solve(mdl::PMModel, evs::PMSimulatorBase.PMEvents, alg::Union{DEAlgorithm,Nothing} = nothing; kwargs...)
    sol_out = Dict{Union{Symbol, Int64}, PMParameterizedSolve.PMSolution}()
    IDs = [instance.ID for instance in evs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in events") : nothing
    # Save parameters and inputs prior to solution. This will let us restore them after. # This could (PROBABLY WILL) have implications/cause problems for parallel solution...
    initP = deepcopy(mdl.parameters.values)
    initU = deepcopy(mdl.states.values)
    initIn = deepcopy(mdl._inputs.values)
    
    for instance in evs.instances
        evi = vcat(instance.inputs, instance.updates)
        cbs = collect_evs(evi, mdl)
        sol_i = PMParameterizedSolve.solve(mdl, alg; callback = cbs, kwargs...)
        # Restore values.
        mdl.parameters.values[:] = initP
        mdl.states.values[:] = initU
        mdl.states.parameters.values[:] = initP
        mdl._inputs.values[:] = initIn
        sol_out[instance.ID] = sol_i
    end

    if length(sol_out) == 1
        return Base.values(sol_out)[1]
    else
        return sol_out
    end
end


    

function PMParameterizedSolve.solve(mdl::PMModel, evs::Vector{PMEvent}, alg::Union{DEAlgorithm,Nothing} = nothing; kwargs...)
    # Save parameters and inputs prior to solution. This will let us restore them after. # This could (PROBABLY WILL) have implications/cause problems for parallel solution...
    initP = deepcopy(mdl.parameters.values)
    initU = deepcopy(mdl.states.values)
    initIn = deepcopy(mdl._inputs.values)
    cbs = collect_evs(evs, mdl)
    sol = solve(mdl, alg; callback = cbs, kwargs...)
    # Restore values.
    mdl.parameters.values[:] = initP
    mdl.states.values[:] = initU
    mdl.states.parameters.values[:] = initP
    mdl._inputs.values[:] = initIn
    return sol
end

function PMParameterizedSolve.solve(mdl::PMModel, data::PMSimulatorBase.DataFrames.AbstractDataFrame, alg::Union{DEAlgorithm,Nothing} = nothing; kwargs...)
    sol_out = Dict{Union{Symbol, Int64}, PMParameterizedSolve.PMSolution}()
    dfevs = PMSimulatorBase.df2evs(data)
    IDs = [instance.ID for instance in dfevs.instances]
    hasduplicates(IDs) ? error("Duplicated IDs detected in dataframe") : nothing
    # Save parameters and inputs prior to solution. This will let us restore them after. # This could (PROBABLY WILL) have implications/cause problems for parallel solution...
    initP = deepcopy(mdl.parameters.values)
    initU = deepcopy(mdl.states.values)
    initIn = deepcopy(mdl._inputs.values)
    for instance in dfevs.instances
        evi = vcat(instance.inputs, instance.updates)
        cbs = collect_evs(evi, mdl)
        sol_i = PMParameterizedSolve.solve(mdl, alg; callback = cbs, kwargs...)
        # push!(sol_out, sol_i)
        sol_out[instance.ID] = sol_i
        # Restore values.
        mdl.parameters.values[:] = initP
        mdl.states.values[:] = initU
        mdl.states.parameters.values[:] = initP
        mdl._inputs.values[:] = initIn
    end
    if length(sol_out) == 1
        return Base.values(sol_out)[1]
    else
        return sol_out
    end
end