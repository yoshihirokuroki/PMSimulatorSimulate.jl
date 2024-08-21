using Revise
using .PMParameterizedBase
using .PMSimulatorBase
using .PMSimulatorSimulate
using .PMParameterizedSensitivity
using Unitful
using DifferentialEquations
using Plots
using CSV
using DataFrames

singh = @model Singh begin
    @IVs t [unit = u"hr", description = "Independent variable (time in hours)", tspan = (0.0, 24.0)]
    Dt = Differential(t)

    @constants begin
        mg_to_g = 1E-3, [unit = u"mg/g", description = "convert mg to g"]
        mol_to_nmol = 1E9, [unit = u"mol/nmol", description = "convert mol to nmol"]
        ug_to_g = 1E-6, [unit = u"μg/g", description = "convert ug to g"]
        mL_to_L = 1E-3, [unit = u"mL/L", description = "convert mL to L"]
        um_to_cm = 1E-4, [unit = u"μm/cm", description = "convert um to cm"]
        day_to_h = 24, [unit = u"d/hr", description = "convert day to h"]
        cm3_to_mm3 = 1E3, [unit = u"cm^3/mm^3", description = "convert cm^3 to mm^3"]
        mm3_to_L = 1E-6, [unit = u"L*mm^-3", description = "convert mm^3 to L"]
        MW = 148781, [unit = "g/mol", description = "T-DM1 molecular weight (Scheuher et al 2022)"]
    end


     @parameters begin 
        BW = 70, [unit = u"kg", description = "Body Weight"]
        k_on_ADC = 0.37, [unit = u"(nmol*L^-1)^-1*hr^-1", description = "2nd-order association rate constant between T-DM1 and HER2 antigen"]
        k_off_ADC = 0.097, [unit = u"hr^-1", description = "dissociation rate constant between T-DM1 and HER2 antigen"]
        k_int_ADC = 0.09, [unit=u"hr^-1", description ="internalization rate of the HER2-T-DM1 complex inside the cell"]
        k_deg_ADC = 0.03, [unit=u"hr^-1", description="proteasomal degradation rate of T-DM1 in endosomal/lysosomal space"]
        k_on_Tub = 0.03, [unit=u"(nmol*L^-1)^-1*hr^-1", description ="2nd-order association rate constant between DM1 catabolites and intracellular tubulin protein"]
        k_off_Tub = 0.0285, [unit=u"hr^-1", description="dissociation rate constant between tubulin catabolites and tubulin protein"]
        Tub_total = 65, [unit=u"nmol*L^-1", description ="total concentration of intracellular tubulin protein"]
        k_dec_ADC = 0.0223, [unit=u"hr^-1" description="non-specific deconjugation rate of T-DM1 from the extracellular space"]
        k_diff_Drug = 0.092, [unit=u"hr^-1", description="bidirectional diffusion rate constant for DM1 catabolites in the intracellular and extracellular space"]
        k_out_Drug = 0, [unit=u"hr^-1", description="active efflux rate constant of DM1 catabolites from the intracellular space to extracellular space"]
        Ag_total = 1660, [unit=u"nmol*L^-1", description="(Ag_total_HER2 3+) total antigen expression levels derived using the number of HER2 receptors per cell and the number of cells packed in a liter volume assuming the volume of a cell as 1 pl"]
        # Parameters associated with tumor disposition of T-DM1
        R_Cap = 8.0*um_to_cm, [unit=u"cm", description="radius of the tumor blood capillary"]
        R_Krogh = 75.0*um_to_cm, [unit=u"cm", description="an average distance between two capillaries"]
        P_ADC = 334*um_to_cm/day_to_h, [unit=u"cm*hr^-1", description="rate of permeability of T-DM1 across the blood vessels"]
        P_Drug = 21000*um_to_cm/day_to_h, [unit=u"cm*hr^-1", description="rate of permeability of DM1 catabolites across the blood vessels"]
        D_ADC = 0.022/day_to_h, [unit=u"mm^2*hr^-1", description="rate of diffusion of T-DM1 across the blood vessels"]
        D_Drug = 0.25/day_to_h, [unit=u"mm^2*hr^-1", description="rate of diffusion of DM1 catabolites across the blood vessels"]
        varepsilon_ADC = 0.24,[description="tumor void volume for T-DM1"]
        varepsilon_Drug = 0.44,[description="tumor void volume for DM1 catabolites"]
        # Parameters associated with systemic pharmacokinetics in Human
        CL_ADC = 0.0043/day_to_h , [unit=u"L*hr^-1*kg^-1", description="central clearance of ADC"]
#        CL_ADC = 6*0.0043/day_to_h , [unit=u"L*hr^-1*kg^-1", description="central clearance of ADC"]
        CLD_ADC = 0.014/day_to_h, [unit=u"L*hr^-1*kg^-1", description="distributional clearance of ADC"]
        V1_ADC = 0.034, [unit=u"L*kg^-1", description="central volume of distribution of ADC"]
        V2_ADC = 0.04, [unit=u"L*kg^-1", description="peripheral volume of distribution of ADC"]
        CL_Drug = 2.23/day_to_h, [unit=u"L*hr^-1*kg^-1", description="central clearance of DM1 catabolites"]
        CLD_Drug = 1.0/day_to_h, [unit=u"L*hr^-1*kg^-1", description="distributional clearance of DM1 catabolites"]
        V1_Drug = 0.034, [unit=u"L*kg^-1", description="central volume of distribution of DM1 catabolites"]
        V2_Drug = 5.0, [unit=u"L*kg^-1", description="peripheral volume of distribution of DM1 catabolites"]
        k_dec_ADC_plasma = 0.241/day_to_h, [unit=u"hr^-1", description="non-specific deconjugation rate constant for T-DM1 in the systemic circulation"]
        # Clinically reported breast cancer-related parameters used to build the translated PK-PD model
        k_growth_exponential = log(2)/(25/day_to_h), [unit=u"hr^-1", description="exponential growth phase of the tumor"]
        k_growth_linear = 621/day_to_h, [unit=u"mm^3*hr^-1", description="linear growth phase of the tumor"]
        Psi = 20.0#,[description="switch between exponential growth and linear growth phases"]
        V_max = 523.8*cm3_to_mm3, [unit=u"mm^3", description="maximum achievable tumor volume"]
        k_kill_max = 0.139/day_to_h, [unit=u"hr^-1" description="linear killing constant (Scheuher et al 2022)"]
        KC_50 = 23.8, [unit = u"nM", description = "concentration of drug (payload) corresponding to a killing rate constant of half maximum value (Scheuher et al 2022)"]
     end

     @variables begin
        (X1_ADC_nmol(t)=0.0), [unit = u"nmol", description = "amount of T-DM1 in the plasma central compartment,[nmol]"]
        (X2_ADC_nmol(t)=0.0), [unit = u"nmol", description = "amount of T-DM1 in the peripheral compartment"]
        (C_ADC_f_ex_nM(t)=0.0), [unit = u"nmol*L^-1", description = "concentration of T-DM1 in the tumor extracellular space"]
        (C_ADC_b_ex_nM(t)=0.0), [unit = u"nmol*L^-1", description = "concentration of T-DM1 bound to tumor cell surface,[nM]"]
        (C_ADC_endolyso_cell_nM(t)=0.0), [unit = u"nmol*L^-1", description = "concentration of internalized ADC within the lysosomal/endosomal space,[nM]"]
        (C_Drug_f_cell_nM(t)=0.0), [unit = u"nmol*L^-1", description = "concentration of unconjugated intracellular free DM1,[nM]"]
        (C_Drug_b_cell_nM(t)=0.0), [unit = u"nmol*L^-1", description = "concentration of unconjugated intracellular tubulin-bound DM1,[nM]"]
        (C_Drug_f_ex_nM(t)=0.0), [unit = u"nmol*L^-1", description = "concentration of free DM1 in the extracellular space of tumor,[nM]"]
        (C1_Drug_nM(t)=0.0), [unit = u"nmol*L^-1", description = "concentration of DM1 in the plasma central compartment,[nM]"]
        (C2_Drug_nM(t)=0.0), [unit = u"nmol*L^-1", description = "concentration of DM1 in the peripheral compartment,[nM]"]
        (DAR(t)=0.0), [description = "average number of DM1 molecules conjugated to Trastuzumab (unitless)"]
        (TV_mm3(t)=1.0), [unit = u"mm^3", description = "tumor volume [mm^3]"]
     end
          
    @observed R_tumor = (3 * TV_mm3/(4*pi))^(1/3)


    @eq Dt(X1_ADC_nmol) ~ (
        -(CL_ADC/V1_ADC)*X1_ADC_nmol - (CLD_ADC/V1_ADC)*X1_ADC_nmol + (CLD_ADC/V2_ADC)*X2_ADC_nmol - k_dec_ADC_plasma*X1_ADC_nmol 
        - ((2*P_ADC*R_Cap)/(R_Krogh^2))*((X1_ADC_nmol/(V1_ADC*BW))*varepsilon_ADC - C_ADC_f_ex_nM)*(TV_mm3*mm3_to_L)
        - ((6*D_ADC)/((R_tumor)^2))*((X1_ADC_nmol/(V1_ADC*BW))*varepsilon_ADC - C_ADC_f_ex_nM)*(TV_mm3*mm3_to_L)
    )

    @eq Dt(X2_ADC_nmol) ~ (CLD_ADC/V1_ADC)*X1_ADC_nmol - (CLD_ADC/V2_ADC)*X2_ADC_nmol

    @eq Dt(C_ADC_f_ex_nM) ~ (
        ((2*P_ADC*R_Cap)/(R_Krogh^2))*((X1_ADC_nmol/(V1_ADC*BW))*varepsilon_ADC - C_ADC_f_ex_nM)
        + ((6*D_ADC)/((R_tumor)^2))*((X1_ADC_nmol/(V1_ADC*BW))*varepsilon_ADC - C_ADC_f_ex_nM) 
        - k_on_ADC*C_ADC_f_ex_nM*((Ag_total - C_ADC_b_ex_nM)/varepsilon_ADC) + k_off_ADC*C_ADC_b_ex_nM - k_dec_ADC*C_ADC_f_ex_nM
    )

    @eq Dt(C_ADC_b_ex_nM) ~ (
        k_on_ADC*C_ADC_f_ex_nM*((Ag_total - C_ADC_b_ex_nM)/varepsilon_ADC)
        - (k_off_ADC + k_int_ADC + k_dec_ADC)*C_ADC_b_ex_nM
    )

    @eq Dt(C_ADC_endolyso_cell_nM) ~ k_int_ADC*C_ADC_b_ex_nM - k_deg_ADC*C_ADC_endolyso_cell_nM

    @eq Dt(C_Drug_f_cell_nM) ~ (
        k_deg_ADC*DAR*C_ADC_endolyso_cell_nM - k_on_Tub*C_Drug_f_cell_nM*(Tub_total - C_Drug_b_cell_nM)
        + k_off_Tub*C_Drug_b_cell_nM - k_out_Drug*C_Drug_f_cell_nM + k_diff_Drug*(C_Drug_f_ex_nM - C_Drug_f_cell_nM)
    )

    @eq Dt(C_Drug_b_cell_nM) ~ k_on_Tub*C_Drug_f_cell_nM*(Tub_total - C_Drug_b_cell_nM) - k_off_Tub*C_Drug_b_cell_nM

    @eq Dt(C_Drug_f_ex_nM) ~ (
        ((2*P_Drug*R_Cap)/(R_Krogh^2))*(C1_Drug_nM*varepsilon_Drug - C_Drug_f_ex_nM)
        + ((6*D_Drug)/((R_tumor)^2))*(C1_Drug_nM*varepsilon_Drug - C_Drug_f_ex_nM)
        + k_out_Drug*C_Drug_f_cell_nM + k_dec_ADC*DAR*(C_ADC_f_ex_nM + C_ADC_b_ex_nM) - k_diff_Drug*(C_Drug_f_ex_nM - C_Drug_f_cell_nM)
    )

    @eq Dt(C1_Drug_nM) ~ (
        -(CL_Drug/V1_Drug)*C1_Drug_nM - (CLD_Drug/V1_Drug)*C1_Drug_nM
        + (CLD_Drug/V1_Drug)*C2_Drug_nM + (X1_ADC_nmol*DAR*k_dec_ADC_plasma)/(V1_Drug*BW)
        + (CL_ADC*DAR*(X1_ADC_nmol/V1_ADC))/(V1_Drug*BW)
        - ((2*P_Drug*R_Cap)/(R_Krogh^2))*(C1_Drug_nM*varepsilon_Drug - C_Drug_f_ex_nM)
        - ((6*D_Drug)/((R_tumor)^2))*(C1_Drug_nM*varepsilon_Drug - C_Drug_f_ex_nM)
    )
    @eq Dt(C2_Drug_nM) ~ (CLD_Drug/V2_Drug)*C1_Drug_nM - (CLD_Drug/V2_Drug)*C2_Drug_nM

    @eq Dt(DAR) ~ -k_dec_ADC_plasma*DAR

    @eq Dt(TV_mm3) ~ (
        ((k_growth_exponential*(1-(TV_mm3/V_max)))/((1+(k_growth_exponential*TV_mm3/k_growth_linear)^Psi^(1/Psi)))
        - (k_kill_max/(KC_50 + (C_Drug_f_cell_nM + C_Drug_b_cell_nM)))*(C_Drug_f_cell_nM + C_Drug_b_cell_nM))*TV_mm3
    )

    # convert model predicted amount [nmol] of T-DM1 in plasma to concentration [nM] and include in solution output
    @observed C1_ADC_nM = X1_ADC_nmol/(V1_ADC*BW);

    # convert model predicted amount [nmol] of T-DM1 in peripheral to concentration [nM]
    @observed C2_ADC_nM = X2_ADC_nmol/(V2_ADC*BW); # nM

    # convert model predicted tumor volume [mm3] to tumor diameter [cm]
    @observed TumorDiameter = 2*((3*(TV_mm3/cm3_to_mm3)/4*pi).^(1/3)); # cm
end;


println(getUnit(singh.parameters,:D_ADC))
println(getDescription(singh.parameters,:D_ADC))
println(getUnit(singh.states,:TV_mm3))


singh.tspan = (0.0, 65.0*singh.constants.day_to_h)
# define dosing events
# use dose = 3.6 mg/kg Q3W in HER2-positive metastic breast cancer patients
dose_in_mgkg = 3.6; # [mg/kg] 
# convert dose to nmol
dose = (dose_in_mgkg*singh.constants.mg_to_g*singh.parameters.BW)./((singh.constants.MW+0.0)/singh.constants.mol_to_nmol); # [nmol] 

DAR_0 = 3.5 # T-DM1 DAR is 3.5
L_tumor = 19; # [mm] initial tumor lesion length
B_tumor = 16; # [mm] = initial tumor lesion breadth
TV_0 = (L_tumor*B_tumor^2)/2; # [mm^3], initial tumor volume
singh.states.DAR = DAR_0
singh.states.TV_mm3 = TV_0;



evs = PMEvent[] # Create an array to hold our dosing events
dosetimes = Float64.([0, 3*7*24, 6*7*24, 9*7*24]); # once every 3 weeks

for dtime in dosetimes # Create dosing events in a loop
    push!(evs, PMInput(time = dtime, amt = dose, input = :X1_ADC_nmol))
    push!(evs, PMInput(time = dtime, amt = DAR_0, input = :DAR))
end

# push!(evs, PMUpdate(time = 3*7*24, quantity = :BW, value = 0.5))
cbssingh = PMSimulatorSimulate.collect_evs(evs, singh);

num_sol = PMSimulatorSimulate.solve(singh, evs, alg=AutoTsit5(Rosenbrock23()));
# convert model time points in h to day
tspan_day = num_sol.t/singh.constants.day_to_h; # convert time from h to day



model_plasma_ADC_Drug_plot = plot(title = "Plasma Human PK (3.6 mg/kg Q3W)", legend=:bottomright, yaxis=:log10,);
plot!(tspan_day, num_sol.C1_ADC_nM, label = "ADC (T-DM1)", color=:blue, lw=3);
plot!(tspan_day, num_sol.C1_Drug_nM, label = "Drug (Payload, DM-1)", color=:green, lw=3);
xlabel!("time (day)"); ylabel!("Concentration (nM)");
ylims!(0.01, 1000);
yticks!([0.01, 0.1, 1, 10, 100, 1000], ["0.01", "0.1", "1", "10", "100", "1000"]);
xticks!([0, 21, 42, 63], ["0", "21", "42", "63"]);
display(model_plasma_ADC_Drug_plot)

model_TotalDrugInTumor_plot = plot(title = "Tumor Payload (3.6 mg/kg Q3W)", xlabel = "time (days)", ylabel="Concentration (nM)");
plot!(tspan_day,  num_sol.C_Drug_f_cell_nM + num_sol.C_Drug_b_cell_nM, label = "Intracellular payload", lw=3);
plot!(tspan_day,  num_sol.C_Drug_f_ex_nM, label = "Extracellular free payload", lw=3, linestyle=:dashdot);
hline!([singh.parameters.KC_50+0.0],lw=1,color=:gray,linestyle=:dash, label="KC 50");
xticks!([0, 21, 42, 63], ["0", "21", "42", "63"]);
plot!(legend = :outerbottom, legendcolumns=3);
display(model_TotalDrugInTumor_plot)

model_TumorDiameter_plot = plot(tspan_day, num_sol.TumorDiameter, xlabel = "time (days)", ylabel="Tumor Diameter (cm)", title = "Tumor diameter (cm)", label = "3.6 mg/kg Q3W", lw=3,legend=:bottomright);
# hline!([R_tumor(singhP.k_growth_linear / singhP.k_growth_exponential)/ cm_to_mm * 2],lw=1,color=:gray,linestyle=:dash, label="exponential to linear growth"); 
# hline!([R_tumor(singhP.V_max/2)/ cm_to_mm * 2],lw=1,color=:gray,linestyle=:dash, label="half maximum growth speed"); 
xticks!([0, 21, 42, 63], ["0", "21", "42", "63"]);
display(model_TumorDiameter_plot)

prob_sens = PMParameterizedSensitivity.ODEForwardSensitivityProblem(singh, evs);

sol_result = DataFrame(
    Time_day = tspan_day, 
    ADC_plasma = num_sol.C1_ADC_nM, 
    Drug_plasma = num_sol.C1_Drug_nM, 
    Drug_free_int_cell = num_sol.C_Drug_f_cell_nM, 
    Drug_bound_int_cell = num_sol.C_Drug_b_cell_nM, 
    Drug_free_ex_cell = num_sol.C_Drug_f_ex_nM, 
    Tumor_diameter = num_sol.TumorDiameter
    )

CSV.write("PMSimulatorSimulate.jl/test/ADC_solve.csv", sol_result)