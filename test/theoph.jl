using DataFrames
using CSV
using TidierData
using PMSimulatorSimulate
include("house.jl");
theoph = DataFrame(CSV.File("exTheoph.csv"));

theoph = @chain theoph begin
    @rename(input = cmt)
    @mutate(input = ifelse(input == 1, :GUT, :nothing))
    @mutate(tinf = 2.0)
end;

theoph_ev = df2evs(theoph);
# house2 = deepcopy(house);
house.states.GUT = 0.0


### CREATE COPY OF MODEL FIRST?
#THIS UPDATES THE MODEL WHICH SEEMS BAD!

cbs = PMSimulator.collect_evs([theoph_ev.instances[1].inputs..., theoph_ev.instances[1].updates...], house);
sol = PMParameterized.solve(house)#, callback = cbs);
# sol2 = PMParameterized.solve(house2);
using Plots
plot(sol.t, sol.GUT)
