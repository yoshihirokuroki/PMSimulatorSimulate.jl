
InputOrUpdate = Union{PMInput, PMUpdate}



indexof(sym::Symbol, syms) = findfirst(isequal(sym), syms)

function getMTKindex(mdl::PMModel, sym::Symbol)
    psyms = Symbol.(parameters(mdl.model))
    ssyms = [x.metadata[ModelingToolkit.VariableSource][2] for x in states(mdl.model)]
    pindex = indexof(sym, psyms)
    sindex = indexof(sym, ssyms)
    if !isnothing(pindex) && !isnothing(sindex)
        error("Found $sym in parameters and states, cannot get index")
    elseif isnothing(pindex) && isnothing(sindex)
        error("Cannot locate $sym in parameters or states")
    end
    indices = [pindex, sindex]
    index = indices[indices.!=nothing][1] # Get only the non-nothing index
    return index
end




function updateParameterOrState!(mdl::PMModel, sym::Symbol, val::Float64)
    if sym in mdl.parameters.names
        mdl.parameters[sym] = val
    elseif sym in mdl.states.names
        mdl.states[sym] = val
    else
        error("Cannot find $sym in parameters or states")
    end
end

function updateInput!(mdl::PMModel, sym::Symbol, val::Float64)
    idx = mdl._inputs.sym_to_val[sym]
    pair = getindex(getfield(mdl._inputs,:values),idx)
    pairnew = Pair(pair.first, val)
    setindex!(getfield(getfield(mdl,:_inputs),:values),pairnew, idx)
end


function generateInputCB(mdl::PMModel, tstart::Float64, tinf::Union{Float64,Nothing}, amt::Float64, inputP::Symbol, input::Symbol)
    cbset = DiscreteCallback[]
    indexP = getMTKindex(mdl, inputP)
    indexS = getMTKindex(mdl, input)
    if tstart == 0.0
        if tinf > 0.0
            updateInput!(mdl, input, amt)
            indexP = getMTKindex(mdl, inputP)
            cbend = PresetTimeCallback(tstart+tinf, (integrator) -> integrator.p[indexP] = integrator.p[indexP] - amt/tinf)
            push!(cbset, cbend)
        else
            # WE NEED TO FIX THINGS HERE!
            updateParameterOrState!(mdl, input, amt)
        end
    else
        if tinf > 0.0
            cbstartinf = PresetTimeCallback(tstart, (integrator) -> integrator.p[indexP] =  integrator.p[indexP] + amt/tinf)
            push!(cbset, cbstartinf)
        elseif tinf < 0.0
            error("Cannot have negative infusion time")
        elseif tinf == 0.0
            cbstartbolus = PresetTimeCallback(tstart, (integrator) -> integrator.u[indexS] += amt)
            push!(cbset, cbstartbolus)
        else
            error("Some other error")
        end
    end
    # Infusion end
    if tinf != 0.0
        indexP = getMTKindex(mdl, inputP)
        cbend = PresetTimeCallback(tstart+tinf, (integrator) -> integrator.p[indexP] = integrator.p[indexP] - amt/tinf)
        push!(cbset, cbend)
    end
    return cbset
end

function generateUpdateCB(mdl::PMModel, update::PMUpdate)
    time = update.time
    quantity = update.quantity
    value = update.value
    if time == 0.0
        updateParameterOrState!(mdl, quantity, value)
        return nothing
    else
        index = getMTKindex(mdl, quantity)
        cbtmp = PresetTimeCallback(time, (integrator) -> integrator.p[index] = value)
        return cbtmp
    end
end





        

    