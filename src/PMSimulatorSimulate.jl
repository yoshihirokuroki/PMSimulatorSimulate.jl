module PMSimulatorSimulate
using PMParameterizedBase
using DiffEqCallbacks
using PMParameterizedSolve
using PMSimulatorBase
using DifferentialEquations
parameters = PMSimulatorBase.ModelingToolkit.parameters
states = PMSimulatorBase.ModelingToolkit.states
PMModel = PMParameterizedBase.PMModel
PMSolution = PMParameterizedSolve.PMSolution
partialSol = PMParameterizedSolve.partialSol
PMEvent = PMSimulatorBase.PMEvent
collect_evs = PMSimulatorBase.collect_evs

include("assemble.jl")
include("solve.jl")

export solve
export collect_evs
export solve!
end
