function ogv = ogv_turb(givens, properties, turbine, last_stage)
% Models the turbine outlet guide vane (OGV) as a diffuser


    % Physical constants 
    R     = properties.R_comb;                 % J/kg/K
    gamma = properties.gamma_comb;
    Cp    = R*(gamma/(gamma-1));

    % Givens
    mdot   = givens.mdot_air;                 % kg/s 
    Rebmin = turbine.Retmin;                  % use turbine Re threshold as a min Re

    % Inlet from last turbine stage (rotor exit / OGV inlet)
    To1  = last_stage.To3;                    % K (total at rotor exit)
    Po1  = last_stage.Po3;                    % Pa (total at rotor exit)
    a1   = last_stage.a3;                     % deg (absolute flow angle into OGV)
    cz   = last_stage.cz;                     % m/s (axial velocity at OGV inlet)

    M1   = last_stage.M2;                  
    rh   = last_stage.rh;                     % m (hub radius)

    % Design-variable bounds
    % x = [sigmaN, M2_exit]
    lb = [0.1,  0.10];                        % min solidity, min exit Mach
    ub = [2.0,  1.00];                        % max solidity, max subsonic Mach

    % initial guess
    sigma_guess = 1.0;
    M2_guess    = 0.40;
    x0          = [sigma_guess, M2_guess];

    % Optimization setup
    options = optimoptions("fmincon",...
        "Algorithm","interior-point",...
        "ConstraintTolerance",1e-5,...
        "EnableFeasibilityMode",true,...
        "SubproblemAlgorithm","cg",...
        "MaxFunctionEvaluations",1e6,...
        "MaxIterations",1e6);

    objectiveFunction = @(x) optimizeStage(x, To1, Po1, gamma, R, Cp, ...
                                           mdot, cz, a1, rh, M1);
    nonlcon = @(x) constraints(x, To1, Po1, gamma, R, Cp, ...
                               mdot, cz, a1, rh, Rebmin, givens, M1);

    [out, ~, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], ...
                                      lb, ub, nonlcon, options);

    % Post-process optimized design
    ogv = postProcess(out, To1, Po1, gamma, R, Cp, ...
                      mdot, cz, a1, rh, givens, M1);
end


% objective: minimize total-pressure loss (maximize Po2/Po1)
function J = optimizeStage(x, To1, Po1, gamma, R, Cp, mdot, cz, a1, rh, M1)
    [sigmaN, M2] = unpack(x);

    % velocity triangle with fixed a2 = 0 deg and fixed M1
    [c1, a2, c2, cz2, T1, To2] = ...
        velocityTriangle_ogv(cz, a1, To1, gamma, R, Cp, M1, M2);

    % diffusion factor & pressure ratio across OGV
    wbarfactor = 1.0; % dimensionless scaling for loss correlation
    [Dfs, wbars, Prs] = dfAndPr(a1, a2, sigmaN, wbarfactor, gamma, M1);

    % objective = total-pressure loss = 1 - Po2/Po1
    Po2 = Prs * Po1;
    Pr  = Po2/Po1;
    J   = 1 - Pr;           % minimize loss to minimize (1 - Pr)
    if ~isfinite(J)
        J = 1e3;
    end
end

% constraints: Re, geometry sanity (M1 is fixed)
function [c, ceq] = constraints(x, To1, Po1, gamma, R, Cp, ...
                                mdot, cz, a1, rh, Rebmin, givens, M1)
    [sigmaN, M2] = unpack(x);

    % velocity triangle with fixed M1
    [c1, a2, c2, cz2, T1, To2] = ...
        velocityTriangle_ogv(cz, a1, To1, gamma, R, Cp, M1, M2);

    % static P1 at OGV inlet using fixed M1
    P1 = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));

    % stator radii with fixed hub, solve for mean & tip
    [rt, rms] = statorRadii(rh, mdot, P1, R, T1, cz);

    % blade dimensions and Re at stator
    [ss, bs, Rebs, bladeCountS] = bladeDims(c1, givens, rh, rt, rms, R, T1, P1, sigmaN); %#ok<ASGLU>

    c   = [];
    ceq = [];

    % Re constraint: Re_b >= Rebmin
    c(end+1) = Rebmin - Rebs;

    % enforce a minimum hub radius
    c(end+1) = 0.001 - rh;

    c   = c(:);
    ceq = ceq(:);
end

% velocity triangle for turbine OGV (diffuser), a2 = 0, fixed M1
function [c1, a2, c2, cz2, T1, To2] = ...
    velocityTriangle_ogv(cz, a1, To1, gamma, R, Cp, M1, M2)

    % inlet: use M1 from last turbine stage
    T1  = To1 / (1 + ((gamma-1)/2)*(M1^2));
    c1  = cz/cosd(a1);  % use cz from last stage

    % exit: enforce a2 = 0 deg (zero swirl)
    a2  = 0;

    % assume stator (no work): To2 ~ To1
    To2 = To1;

    % given M2 (design variable), obtain T2 and c2
    T2  = To2 / (1 + ((gamma-1)/2)*(M2^2));
    c2  = M2 * sqrt(gamma*R*T2);
    cz2 = c2*cosd(a2);      % = c2, since a2 = 0
end

% diffusion factor and pressure ratio (same style as ogv_comp)
function [Dfs, wbars, Prs] = dfAndPr(a1, a2, sigmaS, wbarfactor, gamma, M1)
    % diffusion factor
    Dfs = (1 - (cosd(a1)/cosd(a2)) + ...
          (1/sigmaS)*(cosd(a1)/2)*(tand(a1)-tand(a2)));

    % pressure loss coefficient
    wbars = wbarfactor * ((cosd(a1)/cosd(a2))^2) * ...
             (sigmaS/cosd(a2)) * (0.012 + 0.0004*exp(7.5*Dfs));

    % total-pressure ratio across stator
    Prs = (1 - wbars*(1 - (1 + ((gamma-1)/2)*(M1^2))^(gamma/(1-gamma))));
end

% stator radii with fixed hub (solve for rt, rms)
function [rt, rms] = statorRadii(rh, mdot, P1, R, T1, cz)
    % Use continuity with fixed hub radius to get mean radius & tip radius
    rho1 = (P1/(R*T1));

    A    = mdot/(2*pi*rho1*cz);
    rms  = sqrt(rh^2 + A);
    rt   = sqrt(2*rms^2 - rh^2);
end

% blade dimensions and Reynolds number for stator
function [ss, bs, Rebs, bladeCountS] = bladeDims(c1, givens, rh, rt, rms, R, T1, P1, sigmaS)
    mus  = givens.mu_comb_ref * (T1/givens.mu_comb_tref)^0.7;
    rho1 = (P1/(R*T1));
    hs   = (rt - rh);

    % choose a modest aspect ratio for OGV blades
    bladeAR_target = 2.0;          % h / chord
    bs   = hs/bladeAR_target;      % chord
    Rebs = (rho1*c1*bs/mus);
    ss   = bs/sigmaS;              % pitch
    bladeCountS = 2*pi*rms/ss;
end

% round blade count to nearest odd integer
function bladeCount = roundBladeCount(bladeCount)
    if (mod(floor(bladeCount),2) == 0)
        bladeCount = ceil(bladeCount);
    else
        bladeCount = floor(bladeCount);
    end
end

% unpack design variables
function [sigmaN, M2] = unpack(x)
    sigmaN = x(1);
    M2     = x(2);
end

% post-processing: build OGV struct
function ogv = postProcess(out, To1, Po1, gamma, R, Cp, ...
                           mdot, cz, a1, rh, givens, M1)

    [sigmaN, M2] = unpack(out);

    % recompute velocity triangle with fixed M1
    [c1, a2, c2, cz2, T1, To2] = ...
        velocityTriangle_ogv(cz, a1, To1, gamma, R, Cp, M1, M2);

    % recompute Df, wbar, Prs
    wbarfactor = 1.0;
    [Dfs, wbars, Prs] = dfAndPr(a1, a2, sigmaN, wbarfactor, gamma, M1);

    % static pressures and Cp
    P1  = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    Po2 = Prs * Po1;
    P2  = Po2 * (1+((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));

    Cpresss = 1 - (Po2-P2)/(Po1-P1) - wbars; 

    % radii with fixed hub
    [rt, rms] = statorRadii(rh, mdot, P1, R, T1, cz);

    % blade dims
    [ss, bs, Rebs, bladeCountS] = bladeDims(c1, givens, rh, rt, rms, R, T1, P1, sigmaN);
    bladeCountS_rounded = roundBladeCount(bladeCountS);

    % derived aspect ratio actually used
    hs        = (rt - rh);
    bladeARs  = hs/bs;

    % ogv struct 
    ogv.Po_in         = Po1;
    ogv.To_in         = To1;
    ogv.a_in          = a1;

    ogv.Po_out        = Po2;
    ogv.To_out        = To2;      
    ogv.a_out         = a2;    
    ogv.chi_out       = 0;        

    ogv.rhn           = rh;
    ogv.rtn           = rt;
    ogv.rmn           = rms;
    ogv.bladeARn      = bladeARs;
    ogv.Nn            = bladeCountS_rounded;
    ogv.Reon          = Rebs;
    ogv.sigmaN        = sigmaN;

    ogv.Po_out_Po_in  = Po2/Po1;
    ogv.Mexit         = M2;

   
end
