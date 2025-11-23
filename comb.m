function combustor_out = comb(givens, properties, combustor, compressor)
    %starting points
    M3e_guess = 0.15;
    x0 = [M3e_guess];
    
    %lower and upper bounds
    lb = [0.05];
    ub = [0.2];

    %perform optimization
    options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "ConstraintTolerance",1e-5,...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    "MaxFunctionEvaluations",1e6,...
    "MaxIterations",1e6);
    objectiveFunction = @(x) optimizeComb(x, givens, properties, combustor, compressor);
    nonlcon = @(x) constraints(x);
    [out, objective, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], lb, ub, nonlcon, options);

    %post process and return combustor structure
    combustor_out = analyzeComb(out, givens, properties, combustor, compressor);

    %objective function
    function objective = optimizeComb(x, givens, properties, combustor, compressor)
        combustor_result = analyzeComb(x, givens, properties, combustor, compressor);
        objective = -(combustor_result.Po4/properties.Po3);
    end

    %nonlinear constraint function
    function [c, ceq] = constraints(x)
        c = [];
        ceq = [];
    end

    %primary analysis function
    function combustor_result = analyzeComb(x, givens, properties, combustor, compressor)
        %unpack design variables and inputs
        [M3e] = unpack(x);
        %calculate fuel ratio
        combustor.f = ((properties.To4 - properties.To3)/((combustor.dhR*combustor.etaComb/properties.Cp_comb) - properties.To4)); %double check etaComb application
        %calculate diffuser inlet dims - assuming == comp OGV exit
        combustor.A3 = pi*(compressor.rt^2 - compressor.stages(end).rhs^2);
        combustor.rmi = (compressor.rt + compressor.stages(end).rhs)/2;
        combustor.hi = compressor.rt - compressor.stages(end).rhs;

        %calculate prediffuser exit dims using choice of Ar_prediff
        combustor.he = combustor.hi*(combustor.rmr)*(combustor.Ar_prediff);
        combustor.L_prediff = abs(combustor.he - combustor.hi)/(2*tand(combustor.theta));
        combustor.rme = combustor.rmi / combustor.rmr; %prediffuser and main annular region mean radius, assuming they're equal

        %calculate mach number at prediffuser exit
        combustor.M3pi = compressor.stages(end).M1rel; %not actually M1rel - made this M2 for the ogv only
        M3pe = sym('M3pe');
        implicitSol = (combustor.Ar_prediff) == (combustor.M3pi/M3pe)*((1+((properties.gamma_air-1)/2)*M3pe^2)/(1+((properties.gamma_air-1)/2)*combustor.M3pi^2))^((properties.gamma_air+1)/(2*(properties.gamma_air-1)));
        combustor.M3pe = double(vpasolve(implicitSol, M3pe, [0 1]));

        %calculate diffuser exit dims using M3i since this should model
        %both dump and prediffuser loss
        A3e = sym('A3e');
        term2 = (combustor.M3pi/M3e);
        num = (1+((properties.gamma_air-1)/2)*M3e^2);
        den = (1+((properties.gamma_air-1)/2)*combustor.M3pi^2);
        term3 = (num/den)^((properties.gamma_air+1)/(2*(properties.gamma_air-1)));
        implicitSol = (A3e/combustor.A3) == exp(((properties.gamma_air*combustor.M3pi^2)/2)*((1 - (combustor.A3/A3e))^2 + (1 - (combustor.A3/A3e))^6))*term2*term3;
        combustor.A3e = double(vpasolve(implicitSol, A3e, [0 combustor.A3*100]));

        %calculate diffuser po3 loss and assume adiabatic diffuser
        combustor.Po3e = properties.Po3 * exp(((-properties.gamma_air*combustor.M3pi^2)/2)*((1 - (combustor.A3/combustor.A3e))^2 + (1 - (combustor.A3/combustor.A3e))^6));
        combustor.To3e = properties.To3;
        combustor.dPo_cold = (combustor.Po3e - properties.Po3) - (combustor.dPoLoverPo3*properties.Po3);

        %calculate flame tube length
        combustor.P3e = combustor.Po3e*(1+((properties.gamma_air-1)/2)*M3e^2)^(-properties.gamma_air/(properties.gamma_air-1));
        combustor.T3e = combustor.To3e*(1+((properties.gamma_air-1)/2)*M3e^2)^-1;
        combustor.rho3e = combustor.P3e/(combustor.T3e*properties.R_air);
        combustor.L_ft = (3*givens.mdot_air*combustor.tr)/(combustor.rho3e*combustor.A3e);

        %calculate hot losses
        combustor.plf_hot = (properties.To4/properties.To3)*log(properties.To4/properties.To3);
        term2 = (combustor.A3/combustor.A3e)^2;
        num = (properties.gamma_air/2)*(combustor.M3pi^2);
        den = (1+((properties.gamma_air-1)/2)*combustor.M3pi^2)^(properties.gamma_air/(properties.gamma_air-1));
        term3 = num/den;
        combustor.dPo_hot = -(combustor.plf_hot * term2 * term3) * properties.Po3;

        %calculate Po4
        combustor.Po4 = properties.Po3 + combustor.dPo_cold + combustor.dPo_hot;

        %calculate hole area in liner from prescribed pressure drop
        combustor.plf_liner = combustor.dPoLoverPo3 / (term2*term3);
        combustor.ALh_eff = sqrt(1/combustor.plf_liner)*combustor.A3e;

        %calculate liner area
        combustor.AL = combustor.A3e * combustor.ALr;

        %calculate and check liner static pressure delta using areas
        czL = sym('czL');
        TL = properties.To4 - ((czL^2)/(2*properties.Cp_comb));
        PL = combustor.Po4 * (TL/properties.To4)^(properties.gamma_comb/(properties.gamma_comb-1));
        implicitSol = (PL/(properties.R_comb*TL)) == (combustor.msn*(1+combustor.f)*givens.mdot_air)/(combustor.AL*czL);
        combustor.czL = double(vpasolve(implicitSol, czL, [0 inf]));
        combustor.TL = properties.To4 - ((combustor.czL^2)/(2*properties.Cp_comb));
        combustor.PL = combustor.Po4 * (combustor.TL/properties.To4)^(properties.gamma_comb/(properties.gamma_comb-1)); %static pressure inside of liner

        czA = sym('czA');
        TA = combustor.To3e - ((czA^2)/(2*properties.Cp_air));
        PA = combustor.Po3e * (TA/combustor.To3e)^(properties.gamma_air/(properties.gamma_air-1));
        implicitSol = (PA/(properties.R_air*TA)) == ((1-combustor.msn)*givens.mdot_air)/((combustor.A3e-combustor.AL)*czA);
        combustor.czA = double(vpasolve(implicitSol, czA, [0 inf]));
        combustor.TA = combustor.To3e - ((combustor.czA^2)/(2*properties.Cp_air));
        combustor.PA = combustor.Po3e *(combustor.TA/combustor.To3e)^(properties.gamma_air/(properties.gamma_air-1));

        combustor.linerStiffness = (combustor.PA - combustor.PL)/combustor.PL;

        %calculate liner depth
        combustor.DL = combustor.AL/(2*pi*combustor.rme);

        %calculate total diffuser length
        combustor.L_dump = (combustor.dgr*combustor.hi) + (combustor.rdDLr * combustor.DL);
        combustor.L_diff = combustor.L_prediff + combustor.L_dump;

        %calculate total combustor length
        combustor.L = combustor.L_diff + combustor.L_ft;

        %calculate final geometric quantities
        ri3e = sym('ri3e');
        implicitSol = (combustor.A3e) == pi*((2*combustor.rme - ri3e)^2 - ri3e^2);
        combustor.ri3e = double(vpasolve(implicitSol, ri3e, [0 compressor.rt]));
        combustor.ro3e = 2*combustor.rme - combustor.ri3e;

        %pass analyzed combustor back out as return
        combustor.M3e = M3e;
        combustor_result = combustor;
    end

    %unpacking helper function
    function [M3e] = unpack(x)
        M3e = x(1); %mach after diffuser
    end
end