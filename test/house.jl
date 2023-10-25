using PMParameterizedBase

house = @model house begin
    @IVs t [description = "Independent variable time (hrs)", tspan = (0.0, 120.0)]
    D = Differential(t)
    @parameters begin
        CL = 1.0, [description = "Clearance (1/hr)"]
        VC = 20.0, [description = "Volume of distribution (L)"]
        KA = 1.2, [description = "absorption rate (1/hr)"]
        F1 = 1.0, [description = "Bioavailability fraction"]
        D1 = 2.0, [description = "Infusion duration (hr)"]
        WTCL = 0.75, [description = "Exponent WT on CL"]
        WTVC = 1.0, [description = "Exponent WT on VC"]
        SEXCL = 0.7, [description = "Prop cov effect on CL"]
        SEXVC = 0.85, [description = "Prop cov effect on VC"]
        KIN = 100.0, [description = "Resp prod rate constant (1/hr)"]
        KOUT = 2.0, [description = "Resp elim rate constant (1/hr)"]
        IC50 = 10.0, [description = "Conc giving 50% max resp (ng/mL)"]
        WT = 70.0, [description = "Weight (kg)"]
        SEX = 0, [description = "Covariate female sex"]
    end

    KOUTi = exp(log(KOUT));
    CLi   = exp(log(CL)   + WTCL*log(WT/70) + log(SEXCL)*SEX);
    VCi   = exp(log(VC)   + WTVC*log(WT/70) + log(SEXVC)*SEX);
    KAi   = exp(log(KA));
    RESP_0 = KIN/KOUTi;

    @variables begin
        (GUT(t) = 0.0), [description = "Dosing compartment"]
        (CENT(t) = 0.0), [description = "Central compartment"]
        (RESP(t) = RESP_0), [description = "response"]
    end

    @observed CP = CENT/VCi
    INH =  (CP/(IC50+CP))


    @eq D(GUT) ~ -KAi*GUT;
    @eq D(CENT) ~ KAi*GUT - (CLi/VCi)*CENT;
    @eq D(RESP) ~ KIN*(1-INH) - KOUTi*RESP
end;