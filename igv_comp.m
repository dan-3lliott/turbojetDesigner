function guide_vane_stage = igv_comp(Prreq, givens, properties, compressor, a2_design, cz, To1, Po1, a1)

%givens
mdot = givens.mdot_air; %kg/s
gamma = properties.gamma_air;
mu = givens.mu_air_ref; %Ns/m^2
Tref_mu = givens.mu_air_tref; %K
R = properties.R_air; %J/kgK
Cp = (gamma*R)/(gamma-1); %J/kgK
rt = compressor.rt;

%constraints and requirements
zweif = 0.8;
Remin = 1e5;
Remax = 5e5;

%lower and upper bounds
lb = [Remin * (1e-6)];
ub = [Remax * (1e-6)];

%starting points
Reonguess = 0.1; %reynolds number at nozzle throat
x0 = [Reonguess];

%perform optimization
options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "ConstraintTolerance",1e-4,...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    "MaxFunctionEvaluations",1e6,...
    "MaxIterations",1e6);
objectiveFunction = @(x) optimizeStage(x, R, To1, Po1, cz, a1, Cp, gamma, mdot, zweif, mu, Tref_mu, a2_design, rt);
nonlcon = @(x) constraints(x, R, To1, Po1, cz, a1, Cp, gamma, Prreq, mdot, zweif, mu, Tref_mu, a2_design, rt);
[out, bladeCount, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], lb, ub, nonlcon, options);

%post process and output the necessary data for the optimizer
guide_vane_stage = postProcess(out, R, To1, Po1, cz, a1, Cp, gamma, mdot, zweif, mu, Tref_mu, a2_design, rt);

%stage efficiency optimization function
function Nn = optimizeStage(x, R, To1, Po1, cz, a1, Cp, gamma, mdot, zweif, mu, Tref_mu, a2, rt)
    %unpack design variables
    [Reon] = unpack(x);

    %perform velocity triangle analysis with chosen parameters
    [c1, ctheta1, M1, c2, ctheta2, M2] = velocityTriangleAnalysis(R, cz, a1, a2, Po1, To1, Cp, gamma, mdot);

    %determine nozzle radii
    P1 = Po1*(1+((gamma-1)/2)*M1^2)^(-gamma/(gamma-1));
    T1 = To1*(1+((gamma-1)/2)*M1^2)^-1;
    [rhn, rmn] = nozzleRadii(rt, mdot, P1, R, T1, cz);
    
    %calculate optimum solidities using zweifel criteria
    [sigmanz, sigman] = sigmaZweif(cz, ctheta1, ctheta2, a2, M2, zweif, gamma);

    %calculate blade aspect ratio and pressure ratio for nozzle
    h1 = rt - rhn;
    [bladeARn, Po2, Nn, zetan] = nozzle_blade_aspect_ratio(rmn, Reon, mu, R, Cp, To1, Po1, gamma, M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1, Tref_mu);
end

%nonlinear constraint function
function [c, ceq] = constraints(x, R, To1, Po1, cz, a1, Cp, gamma, Prreq, mdot, zweif, mu, Tref_mu, a2, rt)
    %unpack design variables
    [Reon] = unpack(x);

    %perform velocity triangle analysis with chosen parameters
    [c1, ctheta1, M1, c2, ctheta2, M2] = velocityTriangleAnalysis(R, cz, a1, a2, Po1, To1, Cp, gamma, mdot);

    %determine nozzle radii
    P1 = Po1*(1+((gamma-1)/2)*M1^2)^(-gamma/(gamma-1));
    T1 = To1*(1+((gamma-1)/2)*M1^2)^-1;
    [rhn, rmn] = nozzleRadii(rt, mdot, P1, R, T1, cz);

    %calculate optimum solidities using zweifel criteria
    [sigmanz, sigman] = sigmaZweif(cz, ctheta1, ctheta2, a2, M2, zweif, gamma);

    %calculate blade aspect ratio and pressure ratio for nozzle
    h1 = rt - rhn;
    [bladeARn, Po2, Nn, zetan] = nozzle_blade_aspect_ratio(rmn, Reon, mu, R, Cp, To1, Po1, gamma, M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1, Tref_mu);

    %constrain nozzle pressure ratio
    ceq = (Po2/Po1) - Prreq;
    c = [];
end

%unpacking function
function [Reon] = unpack(x)
    Reon = x(1);
end

%solidity calculation function
function [sigmanz, sigman] = sigmaZweif(cz, ctheta1, ctheta2, a2, M2, zweif, gamma)
    %solve for nozzle solidity
    anmean = atand(((ctheta1+ctheta2)/2)/cz); % stagger angle
    sigmanz = ((ctheta1/ctheta2 - 1)*sind(2*a2)*(((gamma/2)*(M2^2))/((1 + ((gamma-1)/2)*(M2^2))^(gamma/(gamma-1)) - 1)))/(-zweif); % axial solidity
    sigman = sigmanz/cosd(anmean); % regular solidity
end

%radius calculation function for nozzle
function [rhn, rmn] = nozzleRadii(rt, mdot, P1, R, T1, cz)
    %calculate hub and pitchline radii
    rho1 = (P1/(R*T1));
    rhn = sqrt(rt^2 - mdot/(rho1*cz*pi));
    rmn = sqrt((rhn^2 + rt^2)/2);
end

%velocity triangle analysis function
function [c1, ctheta1, M1, c2, ctheta2, M2] = velocityTriangleAnalysis(R, cz, a1, a2, Po1, To1, Cp, gamma, mdot)
    %nozzle inlet
    c1 = cz/cosd(a1);
    ctheta1 = cz*tand(a1);
    T1 = To1 - (c1^2)/(2*Cp);
    M1 = c1/sqrt(gamma*R*T1);

    %nozzle outlet
    To2 = To1; %assuming adiabatic across nozzle
    ctheta2 = cz*tand(a2);
    c2 = cz/cosd(a2);
    T2 = To2 - (c2^2)/(2*Cp);
    M2 = c2/sqrt(gamma*R*T2);
end

function [bladeARn, Po2, Nn, zetan] = nozzle_blade_aspect_ratio(rm, Reon, mu, R, Cp, To1, Po1, gamma, M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1, Tref_mu)
    iterations = 2000; % iterations for both loops
    %blade aspect ratio - nozzle
    PoIterable = Po1;
    for i = 1:iterations
        T2 = To1 * (1 + ((gamma-1)/2)*(M2^2))^-1;
        P2 = PoIterable(i) * (1 + ((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
        rho2 = P2/(R*T2);
        
        mun = mu*(T2/Tref_mu)^0.7;
        nun = mun/rho2;
        
        c2 = sqrt(ctheta2^2 + cz^2);
        on = Reon * (1e6) * nun / c2;
        sn = on / cosd(a2);
        Nn = getIntegerBladeCount(2*pi*rm/sn); %num of blades for nozzle
        sn = (2*pi*rm)/Nn;
        bn = sigman*sn;
        bladeARn = h1/bn;
        
        %stagnation pressure ratio - nozzle
        bzn = sigmanz * sn;
        Cnozzle = 0.993 + 0.021*(bzn/h1);
        
        zetastarn = 1.04 + 0.06*((a1 + a2)/100)^2;
        Dhn = ((2*sn*h1*cosd(a2))/(sn*cosd(a2)+h1));
        Ren = (rho2 * c2 * Dhn) / mun;
        zetan = (zetastarn*Cnozzle - 1)*((10^5)/Ren)^0.25;
        
        num = (1 - ((c2^2)/(2*Cp*To1)) * (1/(1-zetan)));
        den = (1 - ((c2^2)/(2*Cp*To1)));
        Po2 = Po1 - ((1 - (num/den)^(gamma/(gamma-1)))*Po1);
        PoIterable(i+1) = Po2;
    end
end
 
%blade count rounding function
function Nrounded = getIntegerBladeCount(Nunrounded)
    Nrounded = floor(Nunrounded);
    if (mod(Nrounded,2) == 0)
        Nrounded = Nrounded + 1;
    end
end

%post-processing and output function
function compressor_stage = postProcess(x, R, To1, Po1, cz, a1, Cp, gamma, mdot, zweif, mu, Tref_mu, a2, rt)
    %unpack design variables
    [Reon] = unpack(x);

    %perform velocity triangle analysis with chosen parameters
    [c1, ctheta1, M1, c2, ctheta2, M2] = velocityTriangleAnalysis(R, cz, a1, a2, Po1, To1, Cp, gamma, mdot);

    %calculate nozzle radii
    P1 = Po1*(1+((gamma-1)/2)*M1^2)^(-gamma/(gamma-1));
    T1 = To1*(1+((gamma-1)/2)*M1^2)^-1;
    [rhn, rmn] = nozzleRadii(rt, mdot, P1, R, T1, cz);

    %calculate optimum solidities using zweifel criteria
    [sigmanz, sigman] = sigmaZweif(cz, ctheta1, ctheta2, a2, M2, zweif, gamma);

    %calculate blade aspect ratio and pressure ratio for nozzle
    h1 = rt - rhn;
    [bladeARn, Po2, Nn, zetan] = nozzle_blade_aspect_ratio(rmn, Reon, mu, R, Cp, To1, Po1, gamma, M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1, Tref_mu);

    %calculate blade aspect ratio and pressure ratio for rotor
    To2 = To1;

    %for nozzle outputs
    chi1 = a1; %ASSUME NO INCIDENCE
    chi2 = a2 + cartersDeviation((a2-a1), sigman);
    Po2_Po1 = Po2/Po1;

    %create and pass structure
    compressor_stage.Pr = Po2_Po1;
    compressor_stage.Po3 = Po2;
    compressor_stage.To3 = To2;
    compressor_stage.rhr = rhn;
    compressor_stage.rhs = rhn;
    compressor_stage.rmr = 0;
    compressor_stage.rms = rmn;
    compressor_stage.psi = 0;
    compressor_stage.phi = 0;
    compressor_stage.eta = 0;
    compressor_stage.M1rel = M1;
    compressor_stage.a2 = a1;
    compressor_stage.a3 = a2;
    compressor_stage.chi2r = chi1;
    compressor_stage.chi3 = chi2;
    compressor_stage.Dfr = 0;
    compressor_stage.Dfs = 0;
    compressor_stage.sigmaR = 0;
    compressor_stage.sigmaS = sigman;
    compressor_stage.bladeCountR = 0;
    compressor_stage.bladeCountS = Nn;
    compressor_stage.Rebr = 0;
    compressor_stage.Rebs = Reon * (1e6);
    compressor_stage.bladeARr = 0;
    compressor_stage.bladeARs = bladeARn;
    compressor_stage.ANsq = 0;
    compressor_stage.sigmaCent = 0;
    compressor_stage.sigmaBend = 0;
    compressor_stage.Po2 = 0;
    compressor_stage.Wdot = 0;
end

function d = cartersDeviation(flowAngleDelta, sigma)
    m = 1/8;
    d = m*(abs(flowAngleDelta)/sigma);
end

end