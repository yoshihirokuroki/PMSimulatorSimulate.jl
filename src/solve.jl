using DifferentialEquations
function solve(mdl::PMModel, alg::Union{DEAlgorithm,Nothing} = nothing; evs::PMSimulatorBase.PMEvents, kwargs...)
    for instance in evs.instances
        evi = vcat(instance.input, instance.update)
        cbs = collect_evs(evi, mdl)
        sol_i = PMParameterizedSolve.solve(mdl, alg; callback = cbs, kwargs)
    end
end
    