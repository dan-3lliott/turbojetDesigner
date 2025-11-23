function [compressor_stage, rt_comp, N_rpm, a1_comp, cz_comp, M2_diff] = first_comp_stage(Prreq, diffuser, givens, properties, compressor, previous_stage, amax)

%physical constants
R = properties.R_air; %J/kg*K
gamma = properties.gamma_air;
Cp = R*(gamma/(gamma-1));
wbarfactor = compressor.wbarfactor;

%givens
To1 = previous_stage.To3; %K
Po1 = previous_stage.Po3; %Pa
%determine cz given choice of M2 and a1 - relates to diffuser and IGV
mdot = givens.mdot_air; %kg/s
Rebmin = compressor.Rebmin;

%other qties
i = 2.2; %deg NEED TO CHECK THIS
m = 0.25; %carter's rule coefficient

%constraints and requirements
Dfmax = compressor.Dfmax;
Cpressmax = 100; %was 0.60
M1relmax = compressor.Mrelmax;
a3a1offset = compressor.a3a1offset; %deg

%lower and upper bounds
lb = [0.1  0.001  0.1 0.1 0.1 0.1 0];
ub = [amax 0.5    0.9 0.9 2.0 2.0 diffuser.M2max];

%starting points
a1guess = 16; %deg
rtcompGuess = 0.15; %m
phiGuess = 0.6; %unitless
psiGuess = 0.5; %unitless
sigmaRguess = 1.1; %unitless
sigmaSguess = 1; %unitless
M2_diff_guess = diffuser.M2max; %unitless
x0 = [a1guess, rtcompGuess, phiGuess, psiGuess, sigmaRguess, sigmaSguess, M2_diff_guess];

%perform optimization
options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "ConstraintTolerance",1e-5,...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    "MaxFunctionEvaluations",1e6,...
    "MaxIterations",1e6);
objectiveFunction = @(x) optimizeStage(x, To1, Po1, gamma, R, Cp, wbarfactor, a3a1offset);
nonlcon = @(x) constraints(x, amax, mdot, M1relmax, Po1, To1, gamma, R, Cp, Dfmax, Cpressmax, wbarfactor, givens, Rebmin, a3a1offset);
[out, etaST, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], lb, ub, nonlcon, options);
etaST = -etaST; %flip to positive since we're done tricking fmincon

%apply carter's rule to get metal angles
[chi1, chi2r, chi2s, chi3] = metalAngles(i, out, m, To1, gamma, R, Cp);

%output final datapoints
[N_rpm, compressor_stage] = postProcess(out, etaST, mdot, Cp, To1, Po1, gamma, R, chi1, chi2r, chi2s, chi3, wbarfactor, compressor);
a1_comp = out(1);
rt_comp = out(2);
cz_comp = calculateCz(out(7), properties);
M2_diff = out(end);

%optimization function
function etaST = optimizeStage(x, To1, Po1, gamma, R, Cp, wbarfactor, a3a1offset)
    %unpack optimization variables
    [a1, rt, phi, psi, sigmaR, sigmaS, M2_diff] = unpack(x);
    cz = calculateCz(M2_diff, properties);

    %obtain rotor radii
    [rhr, rmr, N] = rotorRadii(rt, a1, mdot, cz, To1, Po1, Cp, gamma, R, phi);

    %perform velocity triangle analysis
    [a2, a3, b1, b2, M1, M1rel, M2, M2rel, M3, T1, T2, To2, T3] = velocityTriangle(cz, a1, a3a1offset, psi, phi, To1, gamma, R, Cp);

    %calculate pressure ratios
    [~, ~, wbarrr, wbars, Prr, Prs] = dfAndPr(a2, a3, b1, b2, sigmaR, sigmaS, wbarfactor, gamma, M1, M1rel, M2, M2rel);

    %calculate stage efficiency
    etaST = -((Prr * Prs)^((gamma-1)/gamma)-1)/((((cz/phi)^2)*psi)/(Cp*To1));
end

%constraint function
    function [c, ceq] = constraints(x, a1max, mdot, M1relmax, Po1, To1, gamma, R, Cp, Dfmax, Cpressmax, wbarfactor, givens, Rebmin, a3a1offset)
    %unpack optimization variables
    [a1, rt, phi, psi, sigmaR, sigmaS, M2_diff] = unpack(x);
    cz = calculateCz(M2_diff, properties);

    %obtain rotor radii
    [rhr, rmr, N] = rotorRadii(rt, a1, mdot, cz, To1, Po1, Cp, gamma, R, phi);

    %perform velocity triangle analysis
    [a2, a3, b1, b2, M1, M1rel, M2, M2rel, M3, T1, T2, To2, T3] = velocityTriangle(cz, a1, a3a1offset, psi, phi, To1, gamma, R, Cp);

    %M1rel max constraint
    c(1) = M1rel - M1relmax;

    %a3 constraint to make sure next stage will be okay
    c(end+1) = a3 - a1max;

    %Df constraints for diffusion factor
    [Dfr, Dfs, wbarr, wbars, Prr, Prs] = dfAndPr(a2, a3, b1, b2, sigmaR, sigmaS, wbarfactor, gamma, M1, M1rel, M2, M2rel);
    c(end+1) = Dfr - Dfmax;
    c(end+1) = Dfs - Dfmax;

    %Cp constraint for pressure rise coefficients
    P1 = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    Po2 = Prr * Po1;
    P2 = Po2 * (1+((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
    Cpressr = 1 - (Po2-P2)/(Po1-P1) - wbarr;
    %c(end+1) = Cpressr - Cpressmax;
    Po3 = Prs * Po2;
    P3 = Po3 * (1+((gamma-1)/2)*(M3^2))^(-gamma/(gamma-1));
    Cpresss = 1 - (Po3-P3)/(Po2-P2) - wbars;
    %c(end+1) = Cpresss - Cpressmax;

    %make sure that compressor is hitting pressure ratio
    ceq(1) = Prreq - (Prr * Prs);

    %calculate stator radii
    [rhs, rms] = statorRadii(rt, a2, mdot, cz, To2, Po2, Cp, gamma, R);

    %constrain blade reynolds number
    [sr, br, Rebr, bladeCountR, ss, bs, Rebs, bladeCountS] = bladeDims(cz, b1, T1, P1, givens, rt, rhr, rhs, rmr, rms, R, a2, T2, P2, sigmaR, sigmaS, compressor);
    c(end+1) = Rebmin - Rebr;
    c(end+1) = Rebmin - Rebs;

    %DEBUGGING - remove re constraint
    %c(end-1:end) = [];
    %enforce min rh
    c(end+1) = 0.001 - rhr;

    %verify that both matrices are column vectors
    c = c(:);
    ceq = ceq(:);
end

%radius calculation function for rotor
function [rhr, rmr, N] = rotorRadii(rt, a1, mdot, cz, To1, Po1, Cp, gamma, R, phi)
    %calculate hub and pitchline radii
    ctheta1 = cz*tand(a1);
    c1 = ctheta1/sind(a1);
    T1 = To1 - (c1^2)/(2*Cp);
    M1 = c1/sqrt(gamma*R*T1);
    P1 = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    rho1 = (P1/(R*T1));
    rhr = sqrt(rt^2 - mdot/(rho1*cz*pi));
    rmr = sqrt((rhr^2 + rt^2)/2);
    %calculate N from flow coefficient
    U = cz/phi;
    N = 30*U/(pi*rmr);
end

%radius calculation function for stator
function [rhs, rms] = statorRadii(rt, a2, mdot, cz, To2, Po2, Cp, gamma, R)
    %calculate hub and pitchline radii
    ctheta2 = cz*tand(a2);
    c2 = ctheta2/sind(a2);
    T2 = To2 - (c2^2)/(2*Cp);
    M2 = c2/sqrt(gamma*R*T2);
    P2 = Po2 * (1+((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
    rho2 = (P2/(R*T2));
    rhs = sqrt(rt^2 - mdot/(rho2*cz*pi));
    rms = sqrt((rhs^2 + rt^2)/2);
end

%blade dimension calculations for reporting and constraining
    function [sr, br, Rebr, bladeCountR, ss, bs, Rebs, bladeCountS] = bladeDims(cz, b1, T1, P1, givens, rt, rhr, rhs, rmr, rms, R, a2, T2, P2, sigmaR, sigmaS, compressor)
    %rotor
    mur = givens.mu_air_ref * (T1/givens.mu_air_tref)^0.7;
    rho1 = (P1/(R*T1));
    w1 = cz/cosd(b1);
    hr = (rt - rhr);
    br = hr/compressor.bladeARmin;
    Rebr = (rho1*w1*br/mur);
    sr = br/sigmaR;
    bladeCountR = 2*pi*rmr/sr;

    %stator
    mus = givens.mu_air_ref * (T2/givens.mu_air_tref)^0.7;
    rho2 = (P2/(R*T2));
    c2 = cz/cosd(a2);
    hs = (rt - rhs);
    bs = hs/compressor.bladeARmin;
    Rebs = (rho2*c2*bs/mus);
    ss = bs/sigmaS;
    bladeCountS = 2*pi*rms/ss;
end

%velocity triangle analysis function assuming constant axial velocity
%returns flow angles and mach numbers for pressure ratio calcs
function [a2, a3, b1, b2, M1, M1rel, M2, M2rel, M3, T1, T2, To2, T3] = velocityTriangle(cz, a1, a3a1offset, psi, phi, To1, gamma, R, Cp)
    %rotor inlet - casing reference
    ctheta1 = cz*tand(a1);
    c1 = ctheta1/sind(a1);
    
    %rotor inlet - rotor reference
    b1 = atand(tand(a1)-(1/phi));
    wtheta1 = ctheta1 - (cz/phi);
    w1 = wtheta1/sind(b1);
    
    %stator inlet - casing reference
    b2 = atand(tand(b1)+(psi/phi));
    ctheta2 = (cz/phi) + cz*tand(b2);

    a2 = atand(tand(b2)+(1/phi));
    c2 = ctheta2/sind(a2);
    
    %stator inlet - rotor reference
    wtheta2 = cz*tand(b2);
    w2 = wtheta2/sind(b2);
    
    %stage exit - casing reference
    a3 = a1 + a3a1offset;
    ctheta3 = cz*tand(a3);
    c3 = ctheta3/sind(a3);
    
    %calculate mach numbers
    T1 = To1 - ((c1)^2)/(2*Cp);
    M1 = (abs(c1)/sqrt(gamma*R*T1));
    M1rel = (abs(w1)/sqrt(gamma*R*T1));
    
    To2 = (cz/phi)*(ctheta2-ctheta1)/Cp + To1;
    T2 = To2 - ((c2)^2)/(2*Cp);
    M2 = (abs(c2)/sqrt(gamma*R*T2));
    M2rel = (abs(w2)/sqrt(gamma*R*T2));

    To3 = To2;
    T3 = To3 - ((c3)^2)/(2*Cp);
    M3 = (abs(c3)/sqrt(gamma*R*T3));
end

%diffusion factor and pressure ratio calculation function
function [Dfr, Dfs, wbarr, wbars, Prr, Prs] = dfAndPr(a2, a3, b1, b2, sigmaR, sigmaS, wbarfactor, gamma, M1, M1rel, M2, M2rel)
    %diffusion factors
    Dfr = (1 - (cosd(b1)/cosd(b2)) + (1/sigmaR)*(cosd(b1)/2)*(tand(b2)-tand(b1)));
    Dfs = (1 - (cosd(a2)/cosd(a3)) + (1/sigmaS)*(cosd(a2)/2)*(tand(a2)-tand(a3)));

    %pressure ratios
    wbarr = wbarfactor * ((cosd(b1)/cosd(b2))^2)*(sigmaR/cosd(b2))*(0.012 + 0.0004*exp(7.5*Dfr));
    wbars = wbarfactor * ((cosd(a2)/cosd(a3))^2)*(sigmaS/cosd(a3))*(0.012 + 0.0004*exp(7.5*Dfs));
    Prr = (1 - wbarr*(1 - (1 + ((gamma-1)/2)*(M1rel^2))^(gamma/(1-gamma))))*(((1+(gamma-1)*(M2^2)/2)/(1+(gamma-1)*(M2rel^2)/2))^(gamma/(gamma-1))) * (((1+(gamma-1)*(M1rel^2)/2)/(1+(gamma-1)*(M1^2)/2))^(gamma/(gamma-1)));
    Prs = (1 - wbars*(1 - (1 + ((gamma-1)/2)*(M2^2))^(gamma/(1-gamma))));
end

%unpacking function
function [a1, rt, phi, psi, sigmaR, sigmaS, M2_diff] = unpack(x)
    a1 = x(1);
    rt = x(2);
    phi = x(3);
    psi = x(4);
    sigmaR = x(5);
    sigmaS = x(6);
    M2_diff = x(7);
end

%carter's rule function
function [chi1, chi2r, chi2s, chi3] = metalAngles(i, out, m, To1, gamma, R, Cp)
    %unpack optimized variables
    [a1, rt, phi, psi, sigmaR, sigmaS, M2_diff] = unpack(out);
    cz = calculateCz(M2_diff, properties);

    %obtain rotor radii
    [rhr, rmr, N] = rotorRadii(rt, a1, mdot, cz, To1, Po1, Cp, gamma, R, phi);

    %rerun velocity triangle analysis
    [a2, a3, b1, b2, M1, M1rel, M2, M2rel, M3, T1, T2, To2, T3] = velocityTriangle(cz, a1, a3a1offset, psi, phi, To1, gamma, R, Cp);
    
    %calculate chi1, chi2s using optimum incidence
    chi1 = b1 + i;
    chi2s = a2 - i;

    %calculate chi2r, chi3 using carter's rule for deviation
    chi2r = (m*chi1*sqrt(sigmaR) + b2)/(1 + m*sqrt(sigmaR));
    chi3 = (m*chi2s*sqrt(sigmaS) - a3)/(m*sqrt(sigmaS) - 1);
end

%cz calculation functiono
function cz = calculateCz(M2_diff, properties)
    %calculate T_in
    T_in = properties.To2*(1+((properties.gamma_air-1)/2)*M2_diff^2)^-1;
    cz = M2_diff*sqrt(properties.gamma_air*properties.R_air*T_in);
end

%output and post-processing function
function [N, compressor_stage] = postProcess(out, etaST, mdot, Cp, To1, Po1, gamma, R, chi1, chi2r, chi2s, chi3, wbarfactor, compressor)
    %unpack optimized variables
    [a1, rt, phi, psi, sigmaR, sigmaS, M2_diff] = unpack(out);
    cz = calculateCz(M2_diff, properties);
    
    %reobtain rotor radii
    [rhr, rmr, N] = rotorRadii(rt, a1, mdot, cz, To1, Po1, Cp, gamma, R, phi);

    %rerun velocity triangle analysis
    [a2, a3, b1, b2, M1, M1rel, M2, M2rel, M3, T1, T2, To2, T3] = velocityTriangle(cz, a1, a3a1offset, psi, phi, To1, gamma, R, Cp);

    %rerun diffusion factor and pressure ratio calcs
    [Dfr, Dfs, wbarr, wbars, Prr, Prs] = dfAndPr(a2, a3, b1, b2, sigmaR, sigmaS, wbarfactor, gamma, M1, M1rel, M2, M2rel);

    %rerun cp calcs
    P1 = Po1 * (1+((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    Po2 = Prr * Po1;
    P2 = Po2 * (1+((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
    Cpressr = 1 - (Po2-P2)/(Po1-P1) - wbarr;
    Po3 = Prs * Po2;
    P3 = Po3 * (1+((gamma-1)/2)*(M3^2))^(-gamma/(gamma-1));
    Cpresss = 1 - (Po3-P3)/(Po2-P2) - wbars;

    %reobtain stator radii
    To2 = T2*(1 + ((gamma-1)/2)*(M2^2));
    To3 = To2; %assuming adiabatic stator does no work
    [rhs, rms] = statorRadii(rt, a2, mdot, cz, To2, Po2, Cp, gamma, R);

    %other blade dimensions and reynold's numbers
    [sr, br, Rebr, bladeCountR, ss, bs, Rebs, bladeCountS] = bladeDims(cz, b1, T1, P1, givens, rt, rhr, rhs, rmr, rms, R, a2, T2, P2, sigmaR, sigmaS, compressor);

    %blade stresses
    [ANsq, sigmaCent, sigmaBend] = bladeStresses(rt, rhr, compressor, N, cz, sigmaR, br, Cp, To1, To3, P1, P2);

    %general outputs
    row1 = {'rm,1 (m)', 'rh,1 (m)', 'rt,1 (m)', 'rm,3 (m)', 'rh,3 (m)', 'N (rpm)', 'cz (m/s)', 'phi', 'psi', 'a1 (deg)', 'etaST (%)', 'Wdot (MW)'};
    rm1 = sqrt(((rt^2) + (rhr^2))/2); %mean radius, rotor
    rm3 = sqrt(((rt^2) + (rhs^2))/2); %mean radius, stator
    N = (30*cz)/(phi*pi*rm1); %rpm, using rotor mean radius
    etaST = etaST; %stage efficiency
    Wdot = mdot*Cp*(To2 - To1); %work required
    row2 = {rm1, rhr, rt, rm3, rhs, N, cz, phi, psi, a1, etaST, Wdot};

    %rotor outputs
    row3 = {'camberR (deg)', 'chi1 (deg)', 'chi2r (deg)', 'b2 (deg)', 'sigmaR', 'wbarr', 'Dfr', 'Cpressr', 'M1rel', 'Po2/Po1', '', ''};
    camberR = chi1 - chi2r;
    row4 = {camberR, chi1, chi2r, b2, sigmaR, wbarr, Dfr, Cpressr, M1rel, Prr, '', ''};

    %stator outputs
    row5 = {'camberS (deg)', 'chi2s (deg)', 'chi3 (deg)', 'a2 (deg)', 'sigmaS', 'wbars', 'Dfs', 'Cpresss', 'M2', 'Po3/Po2', '', ''};
    camberS = chi2s - chi3;
    row6 = {camberS, chi2s, chi3, a2, sigmaS, wbars, Dfs, Cpresss, M2, Prs, '', ''};
    
    %write file
    compressor_stage.Pr = Prr * Prs;
    compressor_stage.Po3 = Po3;
    compressor_stage.To3 = To3;
    compressor_stage.rhr = rhr;
    compressor_stage.rhs = rhs;
    compressor_stage.rmr = rmr;
    compressor_stage.rms = rms;
    compressor_stage.psi = psi;
    compressor_stage.phi = phi;
    compressor_stage.eta = etaST;
    compressor_stage.M1rel = M1rel;
    compressor_stage.a2 = a2;
    compressor_stage.a3 = a3;
    compressor_stage.chi2r = chi2r;
    compressor_stage.chi3 = chi3;
    compressor_stage.Dfr = Dfr;
    compressor_stage.Dfs = Dfs;
    compressor_stage.sigmaR = (br/(2*pi*rmr/roundBladeCount(bladeCountR)));
    compressor_stage.sigmaS = (bs/(2*pi*rms/roundBladeCount(bladeCountS)));
    compressor_stage.bladeCountR = roundBladeCount(bladeCountR);
    compressor_stage.bladeCountS = roundBladeCount(bladeCountS);
    compressor_stage.Rebr = Rebr;
    compressor_stage.Rebs = Rebs;
    compressor_stage.bladeARr = compressor.bladeARmin;
    compressor_stage.bladeARs = compressor.bladeARmin;
    compressor_stage.ANsq = ANsq;
    compressor_stage.sigmaCent = sigmaCent;
    compressor_stage.sigmaBend = sigmaBend;
    compressor_stage.Po2 = Po2;
    compressor_stage.Wdot = Wdot;
    compressor_stage.a1 = a1;
    compressor_stage.degR = (T2 - T1) / (T3 - T1);
end

%blade count rounding function
function bladeCount = roundBladeCount(bladeCount)
    if (mod(floor(bladeCount),2) == 0)
        bladeCount = ceil(bladeCount);
    else
        bladeCount = floor(bladeCount);
    end
end

%blade stress calculation function
function [ANsq, sigmaCent, sigmaBend] = bladeStresses(rt, rhr, compressor, N, cz, sigmaR, br, Cp, To1, To3, P1, P2)
    %AN^2
    Az = 2*pi*(rt^2 - rhr^2);
    ANsq = Az*N^2;
    %centrifugal stress
    term1 = (rt + 2*rhr)/6;
    term2 = (2*rt + rhr)/6;
    sigmaCent = ((N*2*pi/60)^2)*(rt - rhr)*(term1 + compressor.taperRatio*term2);
    %bending stress
    term1 = (cz/(rt*N*pi/30));
    term2 = (Cp*abs(To3-To1))/(Cp*To1);
    term3 = 1/(2*sigmaR);
    term4 = (rt/(br*compressor.thicknessMax));
    p_avg = (P1+P2)/2;
    sigmaBend = p_avg * term1 * term2 * term3 * term4;
end

end