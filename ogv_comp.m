function guide_vane_stage = ogv_comp(Prreq, givens, properties, compressor, a2_design, cz, const_cz_bool, To1, Po1, a1)

%physical constants
R = properties.R_air; %J/kg*K
gamma = properties.gamma_air;
Cp = R*(gamma/(gamma-1));
wbarfactor = compressor.wbarfactor;

%givens
mdot = givens.mdot_air; %kg/s
Rebmin = compressor.Rebmin;

%other qties
i = 2.2; %deg NEED TO CHECK THIS
m = 0.25; %carter's rule coefficient
rt = compressor.rt;

%constraints and requirements
Dfmax = compressor.Dfmax;
Cpressmax = 100; %was 0.60
M1relmax = compressor.Mrelmax;
a2a1offset = a2_design - a1;

%lower and upper bounds
lb = [0.1];
ub = [2.0];

%starting points
sigmaSguess = 1; %unitless
x0 = [sigmaSguess];

%perform optimization
options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "ConstraintTolerance",1e-5,...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    "MaxFunctionEvaluations",1e6,...
    "MaxIterations",1e6);
objectiveFunction = @(x) optimizeStage(x, To1, Po1, gamma, R, Cp, wbarfactor, a2a1offset, cz);
nonlcon = @(x) constraints(x, mdot, M1relmax, Po1, To1, gamma, R, Cp, Dfmax, Cpressmax, wbarfactor, givens, Rebmin, a2a1offset, cz);
[out, totalBladeCount, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], lb, ub, nonlcon, options);

%output final datapoints
guide_vane_stage = postProcess(out, mdot, Cp, To1, Po1, gamma, R, wbarfactor, compressor, cz);

%optimization function
function totalBladeCount = optimizeStage(x, To1, Po1, gamma, R, Cp, wbarfactor, a2a1offset, cz)
    %unpack optimization variables
    [sigmaS] = unpack(x);

    %perform velocity triangle analysis
    [c1, a2, c2, cz2, T1, M1, To2, M2] = velocityTriangle(cz, a1, a2a1offset, To1, gamma, R, Cp, const_cz_bool);

    %calculate pressure ratios
    [Dfs, wbars, Prs] = dfAndPr(a1, a2, sigmaS, wbarfactor, gamma, M1);

    %calculate stator radii
    P1 = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    [rhs, rms] = statorRadii(rt, mdot, P1, R, T1, cz);

    %other blade dimensions and reynold's numbers
    [ss, bs, Rebs, totalBladeCount] = bladeDims(c1, givens, rt, rhs, rms, R, T1, P1, sigmaS, compressor);
end

%constraint function
function [c, ceq] = constraints(x, mdot, M1relmax, Po1, To1, gamma, R, Cp, Dfmax, Cpressmax, wbarfactor, givens, Rebmin, a2a1offset, cz)
    %unpack optimization variables
    [sigmaS] = unpack(x);

    %perform velocity triangle analysis
    [c1, a2, c2, cz2, T1, M1, To2, M2] = velocityTriangle(cz, a1, a2a1offset, To1, gamma, R, Cp, const_cz_bool);

    %calculate stator radii
    P1 = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    [rhs, rms] = statorRadii(rt, mdot, P1, R, T1, cz);

    %M1rel max constraints - rel is reg frame since not moving
    c(1) = M1 - M1relmax;

    %Df constraints for diffusion factor
    [Dfs, wbars, Prs] = dfAndPr(a1, a2, sigmaS, wbarfactor, gamma, M1);
    c(end+1) = Dfs - Dfmax;

    %Cp constraint for pressure rise coefficients
    P1 = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    Po2 = Prs * Po1;
    P2 = Po2 * (1+((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
    Cpresss = 1 - (Po2-P2)/(Po1-P1) - wbars;
    %c(end+1) = Cpresss - Cpressmax;

    %make sure that vane is hitting pressure ratio
    ceq(1) = Prreq - Prs;


    %constrain blade reynolds number
    [ss, bs, Rebs, bladeCountS] = bladeDims(c1, givens, rt, rhs, rms, R, T1, P1, sigmaS, compressor);
    c(end+1) = Rebmin - Rebs;

    %DEBUGGING - remove re constraint
    %c(end-1:end) = [];
    %enforce min rh
    c(end+1) = 0.001 - rhs;

    %verify that both matrices are column vectors
    c = c(:);
    ceq = ceq(:);
    end

%radius calculation function for stator
function [rhs, rms] = statorRadii(rt, mdot, P1, R, T1, cz)
    %calculate hub and pitchline radii
    rho1 = (P1/(R*T1));
    rhs = sqrt(rt^2 - mdot/(rho1*cz*pi));
    rms = sqrt((rhs^2 + rt^2)/2);
end

%blade dimension calculations for reporting and constraining
function [ss, bs, Rebs, bladeCountS] = bladeDims(c1, givens, rt, rhs, rms, R, T1, P1, sigmaS, compressor)
    %stator
    mus = givens.mu_air_ref * (T1/givens.mu_air_tref)^0.7;
    rho1 = (P1/(R*T1));
    hs = (rt - rhs);
    bs = hs/compressor.bladeARmin;
    Rebs = (rho1*c1*bs/mus);
    ss = bs/sigmaS;
    bladeCountS = 2*pi*rms/ss;
end

%velocity triangle analysis function assuming constant axial velocity
%returns flow angles and mach numbers for pressure ratio calcs
function [c1, a2, c2, cz2, T1, M1, To2, M2] = velocityTriangle(cz, a1, a2a1offset, To1, gamma, R, Cp, const_cz_bool)
    %calculate angles
    a2 = a1 + a2a1offset;

    %calculate velocities
    c1 = cz/cosd(a1);
    if const_cz_bool
        cz2 = cz; %assuming cz const
        c2 = cz2/cosd(a2);
    else
        c2 = c1;
        cz2 = c2*cosd(a2); %assuming cz changes
    end

    %calculate mach numbers
    T1 = To1 - ((c1)^2)/(2*Cp);
    M1 = (abs(c1)/sqrt(gamma*R*T1));

    To2 = To1;
    T2 = To2 - ((c2)^2)/(2*Cp);
    M2 = (abs(c2)/sqrt(gamma*R*T2));
end

%diffusion factor and pressure ratio calculation function
function [Dfs, wbars, Prs] = dfAndPr(a1, a2, sigmaS, wbarfactor, gamma, M1)
    %diffusion factors
    Dfs = (1 - (cosd(a1)/cosd(a2)) + (1/sigmaS)*(cosd(a1)/2)*(tand(a1)-tand(a2)));

    %pressure ratios
    wbars = wbarfactor * ((cosd(a1)/cosd(a2))^2)*(sigmaS/cosd(a2))*(0.012 + 0.0004*exp(7.5*Dfs));
    Prs = (1 - wbars*(1 - (1 + ((gamma-1)/2)*(M1^2))^(gamma/(1-gamma))));
end

%unpacking function
function [sigmaS] = unpack(x)
    sigmaS = x(1);
end

%carter's rule function
function [chi1, chi2] = metalAngles(i, out, m, To1, gamma, R, Cp, cz)
    %unpack optimized variables
    [sigmaS] = unpack(out);

    %rerun velocity triangle analysis
    [c1, a2, c2, cz2, T1, M1, To2, M2] = velocityTriangle(cz, a1, a2a1offset, To1, gamma, R, Cp, const_cz_bool);
    
    %calculate chi1, chi2 using optimum incidence and carter's rule for
    %deviation
    chi1 = a1 - i;
    chi2 = (m*chi1*sqrt(sigmaS) - a2)/(m*sqrt(sigmaS) - 1);
end

%output and post-processing function
function compressor_stage = postProcess(out, mdot, Cp, To1, Po1, gamma, R, wbarfactor, compressor, cz)
    %unpack optimized variables
    [sigmaS] = unpack(out);

    %rerun velocity triangle analysis
    [c1, a2, c2, cz2, T1, M1, To2, M2] = velocityTriangle(cz, a1, a2a1offset, To1, gamma, R, Cp, const_cz_bool);

    %rerun diffusion factor and pressure ratio calcs
    [Dfs, wbars, Prs] = dfAndPr(a1, a2, sigmaS, wbarfactor, gamma, M1);

    %rerun cp calcs 
    P1 = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    Po2 = Prs * Po1;
    P2 = Po2 * (1+((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
    Cpresss = 1 - (Po2-P2)/(Po1-P1) - wbars;

    %reobtain stator radii
    [rhs, rms] = statorRadii(rt, mdot, P1, R, T1, cz);

    %other blade dimensions and reynold's numbers
    [ss, bs, Rebs, bladeCountS] = bladeDims(c1, givens, rt, rhs, rms, R, T1, P1, sigmaS, compressor);

    %blade metal angles
    [chi1, chi2] = metalAngles(i, out, m, To1, gamma, R, Cp, cz);

    %write file
    compressor_stage.Pr = Prs;
    compressor_stage.Po3 = Po2;
    compressor_stage.To3 = To2;
    compressor_stage.rhr = rhs;
    compressor_stage.rhs = rhs;
    compressor_stage.rmr = 0;
    compressor_stage.rms = rms;
    compressor_stage.psi = 0;
    compressor_stage.phi = 0;
    compressor_stage.eta = 0;
    compressor_stage.M1rel = M2; %USING THIS FOR COMBUSTOR CODE ONLY
    compressor_stage.a2 = a1;
    compressor_stage.a3 = a2;
    compressor_stage.chi2r = chi1;
    compressor_stage.chi3 = chi2;
    compressor_stage.Dfr = 0;
    compressor_stage.Dfs = Dfs;
    compressor_stage.sigmaR = 0;
    compressor_stage.sigmaS = (bs/(2*pi*rms/roundBladeCount(bladeCountS)));
    compressor_stage.bladeCountR = 0;
    compressor_stage.bladeCountS = roundBladeCount(bladeCountS);
    compressor_stage.Rebr = 0;
    compressor_stage.Rebs = Rebs;
    compressor_stage.bladeARr = 0;
    compressor_stage.bladeARs = compressor.bladeARmin;
    compressor_stage.ANsq = 0;
    compressor_stage.sigmaCent = 0;
    compressor_stage.sigmaBend = 0;
    compressor_stage.Po2 = 0;
    compressor_stage.Wdot = 0;
    compressor_stage.a1 = a1;
    compressor_stage.degR = 0;
end

%blade count rounding function
function bladeCount = roundBladeCount(bladeCount)
    if (mod(floor(bladeCount),2) == 0)
        bladeCount = ceil(bladeCount);
    else
        bladeCount = floor(bladeCount);
    end
end

end