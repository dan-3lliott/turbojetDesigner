clear all
clc

%=====CONSTANTS AND INPUTS=====

%top-level givens
givens.M = 0.80;
givens.mdot_air = 5.02; %kg/s
givens.mu_air_ref = 1.825e-5; %kg/m*s
givens.mu_air_tref = 273.15+20; %K
givens.mu_comb_ref = 5.6e-5; %kg/m*s
givens.mu_comb_tref = 1550; %K
properties.Ta = 249; %K
properties.Pa = 45.1e3; %Pa
properties.Rbar = 8314; %J/gmolK
properties.gamma_air = 1.4; %K
properties.MW_air = 28.97; %g/mol
properties.R_air = (properties.Rbar/properties.MW_air);
properties.Cp_air = (properties.R_air)*(properties.gamma_air/(properties.gamma_air-1));
properties.gamma_comb = 1.28;
properties.MW_comb = 28.7; %g/mol
properties.R_comb = (properties.Rbar/properties.MW_comb);
properties.Cp_comb = (properties.R_comb)*(properties.gamma_comb/(properties.gamma_comb-1));
properties.To4 = 1680; %K

%assumptions - diffuser
diffuser.Rd = 0.985; %pressure recovery coefficient
diffuser.M2max = 0.5;

%assumptions - compressor
compressor.Dfmax = 0.6;
compressor.Rebmin = 2e5; %LEAVING THIS LOWER FOR NOW
compressor.Mrelmax = 0.79;
compressor.bladeARmin = 1.1;
compressor.bladeARmax = 1.6;
compressor.wbarfactor = 2.3;
compressor.taperRatio = 0.8;
compressor.thicknessMax = 0.1; %proportion of chord b
compressor.ai = 0; %deg - coming in from the diffuser before the IGV
compressor.af = 0; %deg - leaving after the OGV
compressor.a3StageMax = 45; %deg
compressor.a3a1offset = 2.5; %deg
compressor.Pr_igv = 0.9921; %stagnation pressure ratio for inlet guide vane
compressor.Pr_ogv = 0.9917; %same thing for ogv
compressor.Prreq = 10.9;
compressor.Nstages = 7; %PICK THIS!!!!!
compressor.Pr = 1; %initial value only
compressor.W_total = 0; %initial value only

%assumptions - combustor
combustor.etaComb = 0.99;
combustor.dhR = 43.8e6; %J/kg
combustor.tr = 2.1e-3; %sec
combustor.msn = 0.65; %snout air mass flow ratio
combustor.dPoLoverPo3 = 0.041;
combustor.dgr = 0.8; %dump gap ratio - minimizing length
combustor.rdDLr = 0.25; %ratio of dump radius to liner depth - minimizing length
combustor.Ar_prediff = 1.3; %prediffuser area ratio, manually select
combustor.theta = 10; %deg
combustor.rmr = 0.9; %mean radius ratio from compressor -> combustor, manually select
combustor.ALr = 0.375; %liner area ratio, manually select to target stiffness

%assumptions - turbine
turbine.taperRatio = 0.5;
turbine.thicknessMax = 0.2; %proportion of chord b
turbine.zweif = 0.8;
turbine.Retmin = 1e5;
turbine.Retmax = 5e5;
turbine.etaShaft = 0.995;
turbine.aFinal = 0; %deg
turbine.Nstages = 1; %PICK THIS!!!!!

%assumptions - nozzle
nozzle.gamma = 1.35;
nozzle.eta_n = 0.95;

%calculate stagnation conditions at inlet
properties.To1 = properties.Ta*(1+((properties.gamma_air-1)/2)*givens.M^2);
properties.To2 = properties.To1; %assuming adiabatic diffuser
properties.Po1 = properties.Pa*(1+((properties.gamma_air-1)/2)*givens.M^2)^(properties.gamma_air/(properties.gamma_air-1));
properties.Po2 = properties.Po1*diffuser.Rd;

%=====OPTIMIZATION=====
[out, thrust, finalEngine] = optimizeEngine(givens, properties, diffuser, compressor, combustor, turbine, nozzle);

function [out, thrust, finalEngine] = optimizeEngine(givens, properties, diffuser, compressor, combustor, turbine, nozzle)
    %initial guess
    Pr1_comp_guess = 1.53; %pressure ratio for first compressor stage
    Wr_turb_guess = 1.0; %work ratio for first turbine stage
    x0 = [Pr1_comp_guess, Wr_turb_guess];
    
    %upper and lower bounds
    lb = [(compressor.Prreq^(1/compressor.Nstages)) 0];
    ub = [1.7                                       2.0];
    
    %cache variable
    x_last = [];
    diffuser_last = [];
    compressor_last = [];
    combustor_last = [];
    turbine_last = [];
    nozzle_last = [];
    
    %optimization functions
    options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "ConstraintTolerance",1e-3,...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    "MaxFunctionEvaluations",1e6,...
    "MaxIterations",1e6);
    objectiveFunction = @(x) engineThrust(x);
    nonlcon = @(x) constraints(x);
    [out, thrust, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], lb, ub, nonlcon, options);
    thrust = -thrust; %flip to positive since we're done tricking fmincon

    %save engine and output final graphs/figures
    finalEngine = saveEngine(givens, properties, diffuser_last, compressor_last, combustor_last, turbine_last, nozzle_last);
    
    %objective function
    function thrust = engineThrust(x)
        %ensure current cache is up-to-date
        updateCache(x);
        thrust = -nozzle_last.thrust;
    end
    
    %nonlinear constraint function
    function [c, ceq] = constraints(x)
        %ensure current cache is up-to-date
        updateCache(x);
        %unpack design variables
        [Pr1_comp, Wr_turb] = unpack(x);

        %=====COMPRESSOR CONSTRAINTS=====
        %constrain total compressor pressure ratio
        ceq(1) = compressor.Prreq - compressor_last.Pr;
        %constrain number of compressor stages to be an integer
        %ceq(1) = mod(Nstages_comp,1);
        %=====TURBINE CONSTRAINTS=====
        %constrain number of turbine stages to be an integer
        %ceq(3) = mod(Nstages_turb,1);
        c = [];
    end
    
    %cache updating function
    function updateCache(x)
        %reuse cached result only if x matches
        if isempty(x_last) || any(x~=x_last)
            %unpack optimization variables
            [Pr1_comp, Wr_turb] = unpack(x);
            %cache component performances

            %=====COMPRESSOR=====
            compressor_last = comp(Pr1_comp, diffuser, givens, properties, compressor);
            properties.Po3 = compressor_last.Po3;
            properties.To3 = compressor_last.To3;

            %=====COMBUSTOR=====O
            combustor_last = comb(givens, properties, combustor, compressor_last);
            properties.Po4 = combustor_last.Po4;

            %=====TURBINE=====
            %set up turbine requirements from compressor and combustor
            turbine.W_total = compressor_last.W_total;
            turbine.Wr1 = Wr_turb;
            
            %cz now chosen inside first_turb_stage (as DV), not from compressor
            turbine_last = turb(0, turbine.Nstages, compressor_last.N, givens, properties, turbine);
            turbine_last.ogv = ogv_turb(givens, properties, turbine, turbine_last.allStages(end));

            %=====NOZZLE=====
            nozzle.To5   = turbine_last.ogv.To_out;
            nozzle.Po5   = turbine_last.ogv.Po_out;
            nozzle.cp_n  = (nozzle.gamma*properties.R_comb)/(nozzle.gamma - 1);
            nozzle_last = nozz(givens, properties, nozzle, combustor_last);

            %=====DIFFUSER=====
            diffuser_last = diffu(diffuser, givens, properties, compressor_last);

            %design variables
            x_last = x;
        end
    end

    %unpacking function
    function [Pr1_comp, Wr_turb] = unpack(x)
        Pr1_comp = x(1);
        Wr_turb = x(2);
    end
end
