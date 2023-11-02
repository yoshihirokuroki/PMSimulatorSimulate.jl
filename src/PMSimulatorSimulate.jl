module PMSimulatorSimulate
using PMParameterizedBase
using PMParameterizedSolve
import PMParameterizedSolve: solve
import PMParameterizedSolve: solve!
using PMSimulatorBase
using DifferentialEquations
parameters = PMSimulatorBase.ModelingToolkit.parameters
states = PMSimulatorBase.ModelingToolkit.states
PMModel = PMParameterizedBase.PMModel
PMSolution = PMParameterizedSolve.PMSolution
partialSol = PMParameterizedSolve.partialSol
PMEvent = PMParameterizedBase.PMEvent
collect_evs = PMSimulatorBase.collect_evs
include("solve.jl")

export solve
export collect_evs
export solve!
end
