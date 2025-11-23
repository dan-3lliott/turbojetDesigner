function turbine_stage = first_turb_stage(N, givens, properties, turbine, previous_stage)

% ---- Givens / properties ----
mdot   = givens.mdot_air;            % kg/s
To1    = previous_stage.To3;         % K
Po1    = previous_stage.Po3;         % Pa
gamma  = properties.gamma_comb;     
R      = properties.R_comb;          % J/kg/K
Cp     = (gamma*R)/(gamma-1);        % J/kg/K

etaTR  = turbine.etaShaft;          
zweif  = turbine.zweif;             
Remin  = turbine.Retmin;           
Remax  = turbine.Retmax;            

% constants
a1   = previous_stage.a3;            % deg
Wout = turbine.WperStage;            % [W] portion of shaft work for this stage

%M2 = 1 for the first stage
M2   = 1.0;

% ---------- optimizer bounds and initial guess ----------
a3min  = -45;  
a3max  =  45;

rh_min = 0.05;
rh_max = 0.15;

cz_min   = 150;    % m/s 
cz_max   = 200;   % m/s

lb = [a3min  rh_min   Remin   Remin  cz_min];
ub = [a3max  rh_max   Remax   Remax  cz_max];

% Initial guesses
a3guess   = -6;        % deg
rh_guess  = 0.05;      % m
Reonguess = 1e5;    % Re at nozzle throat
Reorguess = 1e5;    % Re at rotor "throat"
cz_guess  = 180;       % m/s axial velocity guess

% Clamp initial guesses to bounds
a3guess   = min(max(a3guess,   a3min), a3max);
rh_guess  = min(max(rh_guess,  rh_min), rh_max);
Reonguess = min(max(Reonguess, Remin),  Remax);
Reorguess = min(max(Reorguess, Remin),  Remax);
cz_guess  = min(max(cz_guess,  cz_min), cz_max);

x0 = [a3guess rh_guess Reonguess Reorguess cz_guess];

% Optional fixed exit angle for last stage
if isfield(turbine, 'force_a3')
    a3_fixed = turbine.force_a3;
    lb(1)    = a3_fixed;
    ub(1)    = a3_fixed;
    x0(1)    = a3_fixed;
end

% -------------------------------------------------------
% perform optimization (M2 is fixed inside, not a DV)
objectiveFunction = @(x) optimizeStage(x, R, To1, Po1, a1, Cp, ...
                                       gamma, N, mdot, zweif, givens, ...
                                       etaTR, Wout, M2);

nonlcon = @(x) constraints(x, R, To1, Po1, a1, Cp, ...
                           gamma, N, mdot, zweif, givens, ...
                           etaTR, Wout, M2);

    options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "ConstraintTolerance",1e-3,...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    "MaxFunctionEvaluations",1e6,...
    "MaxIterations",1e6);

[opt_x, ~, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], ...
                                    lb, ub, nonlcon, options);

% return stage struct using postProcess
turbine_stage = postProcess(opt_x, R, To1, Po1, a1, Cp, ...
                            gamma, N, mdot, zweif, givens, ...
                            etaTR, Wout, M2, turbine);

end  % main function

% =================== nested functions ===================

% ---- objective: maximize etaST (via minimizing -etaST) ----
function etaST = optimizeStage(x, R, To1, Po1, a1, Cp, ...
                               gamma, N, mdot, zweif, givens, ...
                               etaTR, Wout, M2)
    etaST = 1e6; % penalty default
    try
        % unpack design variables
        [a3, rh, Reon, Reor, cz] = unpack(x);

        % velocity triangles + radii based on rh (M2 fixed to 1)
        [c1, ctheta1, c2, ctheta2, a2, ...
         w2, wtheta2, b2, w3, wtheta3, b3, ...
         c3, ctheta3, M1, rtn, rmn, U, rho1, P1] = ...
            velocityTriangleAnalysis(cz, a1, a3, To1, Po1, M2, Cp, ...
                                     gamma, R, mdot, rh, N);

        if ~all(isfinite([rtn, rmn, U, rho1, P1]))
            return;
        end

        % psi, phi
        [psi, phi, ~] = psiPhiR(cz, U, ctheta3, ctheta2, b2, b3);
        if ~all(isfinite([psi, phi]))
            return;
        end

        % Zweifel solidities
        [sigmanz, sigmarz, sigman, sigmar] = ...
            sigmaZweif(w3, To1, U, Cp, R, cz, ctheta1, ctheta2, ...
                       wtheta2, wtheta3, a2, b3, M2, psi, c3, zweif, gamma);
        if ~all(isfinite([sigmanz, sigmarz, sigman, sigmar]))
            return;
        end

        % Nozzle loss model
        h1 = rtn - rh;
        [bladeARn, Po2, Nn, zetan] = ...
            nozzle_blade_aspect_ratio(rmn, Reon, givens, R, Cp, To1, Po1, ...
                                      gamma, M2, sigman, sigmanz, a1, a2, ...
                                      ctheta2, cz, h1);
        if ~all(isfinite([bladeARn, Po2, Nn, zetan])) || Po2 <= 0
            return;
        end

        % Rotor geometry & loss model
        To2 = To1;
        [rtr, rmr] = rotorRadii(cz, rh, mdot, To2, Po2, M2, gamma, R);
        h2 = rtr - rh;
        if ~all(isfinite([rtr, rmr, h2])) || h2 <= 0
            return;
        end

        [bladeARr, Prtot, To2, To3, rho3, Nr, M2rel, zetar, br] = ...
            rotor_blade_aspect_ratio(U, sigmar, sigmarz, rmr, Reor, ...
                                     givens, R, Cp, To1, Po1, Po2, ...
                                     gamma, M2, c2, w2, w3, h2, c3, ...
                                     psi, b2, b3);
        if ~(isfinite(Prtot)) || Prtot <= 0 || ...
           ~all(isfinite([To2, To3, rho3, Nr, M2rel, zetar, bladeARr, br]))
            return;
        end

        % Stage efficiency (Euler / polytropic relation)
        num = psi*(gamma-1)*((U^2)/(gamma*R*To1));
        den = Prtot^((1-gamma)/gamma) - 1;
        if ~(isfinite(num)) || ~(isfinite(den)) || den == 0
            return;
        end

        eta_tmp = -num/den;  % negative so fmincon maximizes
        if isfinite(eta_tmp)
            etaST = eta_tmp;
        end
    catch
        etaST = 1e6;
    end
end

% ---- nonlinear constraints: enforce work requirement *without* guards ----
function [c, ceq] = constraints(x, R, To1, Po1, a1, Cp, ...
                                gamma, N, mdot, ~, givens, ...
                                etaTR, Wout, M2)

    % no inequality constraints
    c = [];

    % unpack design variables
    [a3, rh, ~, ~, cz] = unpack(x);

    % get U, ctheta2, ctheta3 from the velocity triangles
    [c1, ctheta1, c2, ctheta2, a2, ...
     w2, wtheta2, b2, w3, wtheta3, b3, ...
     c3, ctheta3, M1, rtn, rmn, U, rho1, P1] = ...
        velocityTriangleAnalysis(cz, a1, a3, To1, Po1, M2, Cp, ...
                                 gamma, R, mdot, rh, N);

    % compute actual stage work
    Wdot = calculateWorkDone(etaTR, mdot, U, ctheta3, ctheta2);

    % equality constraint: require Wdot = Wout
    ceq = Wdot - Wout;
end


% ---- work calculation ----
function Wdot = calculateWorkDone(etaTR, mdot, U, ctheta3, ctheta2)
    % positive Wdot for extracted shaft work
    Wdot = -(etaTR * mdot * U * (ctheta3 - ctheta2));
end

% ---- unpacking function for design variables ----
function [a3, rh, Reon, Reor, cz] = unpack(x)
    a3   = x(1);
    rh   = x(2);
    Reon = x(3);
    Reor = x(4);
    cz   = x(5);
end

% ---- solidity calculation (Zweifel) ----
function [sigmanz, sigmarz, sigman, sigmar] = ...
    sigmaZweif(w3, To1, U, Cp, R, cz, ctheta1, ctheta2, ...
               wtheta2, wtheta3, a2, b3, M2, psi, c3, zweif, gamma)

    % nozzle solidity
    anmean = atand(((ctheta1+ctheta2)/2)/cz);
    sigmanz = ((ctheta1/ctheta2 - 1)*sind(2*a2)*(((gamma/2)*(M2^2))/ ...
               ((1 + ((gamma-1)/2)*(M2^2))^(gamma/(gamma-1)) - 1))) / (-zweif);
    sigman  = sigmanz/cosd(anmean);

    % rotor solidity (uses M3rel)
    To2   = To1;
    To3   = To2 + (psi * U^2)/Cp;
    T3    = To3 - (c3^2)/(2*Cp);
    M3    = c3/sqrt(gamma*R*T3);
    M3rel = (M3/c3)*w3;

    brmean = atand(((wtheta2+wtheta3)/2)/cz);
    sigmarz = ((wtheta2/wtheta3 - 1)*sind(2*b3)*(((gamma/2)*(M3rel^2))/ ...
               ((1 + ((gamma-1)/2)*(M3rel^2))^(gamma/(gamma-1)) - 1))) / (zweif);
    sigmar  = sigmarz/cosd(brmean);
end

% ---- stator (nozzle) radii and U with rh as DV, using Po1 ----
function [rtn, rmn, U, rho1, P1] = ...
    nozzleRadiiAndU(cz, a1, To1, Cp, gamma, R, mdot, rh, N, Po1)

    % inlet state
    c1  = cz/cosd(a1);
    T1  = To1 - c1^2/(2*Cp);
    M1  = c1/sqrt(gamma*R*T1);

    % Use Po1 to get P1 (no more 1 bar hack)
    P1  = Po1 * (1 + ((gamma-1)/2)*M1^2)^(-gamma/(gamma-1));
    T1  = To1 * (1 + ((gamma-1)/2)*M1^2)^(-1);
    rho1 = P1/(R*T1);

    % annulus area from continuity (given rh)
    A   = mdot/(2*pi*rho1*cz);
    rmn = sqrt(rh^2 + A);
    rtn = sqrt(2*rmn^2 - rh^2);
    U   = rmn*N*pi/30;
end

% ---- rotor radii with constant hub (for given Po2, M2) ----
function [rtr, rmr] = rotorRadii(cz, rh, mdot, To2, Po2, M2, gamma, R)
    P2   = Po2 * (1 + ((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
    T2   = To2 * (1 + ((gamma-1)/2)*(M2^2))^(-1);
    rho2 = P2/(R*T2);
    rmr = sqrt(rh^2 + (mdot/(2*rho2*cz*pi)));
    rtr = sqrt(2*rmr^2 - rh^2);
end

% ---- velocity triangles for first stage (M2 fixed, rh passed in) ----
function [c1, ctheta1, c2, ctheta2, a2, ...
          w2, wtheta2, b2, w3, wtheta3, b3, ...
          c3, ctheta3, M1, rtn, rmn, U, rho1, P1] = ...
    velocityTriangleAnalysis(cz, a1, a3, To1, Po1, M2, Cp, ...
                             gamma, R, mdot, rh, N)

    % nozzle inlet
    c1       = cz/cosd(a1);
    ctheta1  = cz*tand(a1);
    T1       = To1 - (c1^2)/(2*Cp);
    M1       = c1/sqrt(gamma*R*T1);

    % nozzle exit (M2 fixed)
    To2      = To1; 
    T2       = To2*(1 + ((gamma-1)/2)*M2^2)^-1;
    c2       = M2 * sqrt(gamma*R*T2);
    a2       = acosd(cz/c2);
    ctheta2  = cz*tand(a2);

    % rotor outlet (stationary frame)
    ctheta3  = cz*tand(a3);
    c3       = cz/cosd(a3);

    % radii + U from rh and continuity
    [rtn, rmn, U, rho1, P1] = nozzleRadiiAndU(cz, a1, To1, Cp, ...
                                              gamma, R, mdot, rh, N, Po1);

    % rotor frame velocities
    wtheta3  = ctheta3 - U;  
    w3       = sqrt(wtheta3^2 + cz^2);  
    b3       = atand(wtheta3/cz);

    wtheta2  = ctheta2 - U;  
    w2       = sqrt(wtheta2^2 + cz^2);  
    b2       = atand(wtheta2/cz);
end

% ---- psi, phi, reaction ----
function [psi, phi, Rdeg] = psiPhiR(cz, U, ctheta3, ctheta2, b2, b3)
    phi  = cz/U; 
    psi  = (ctheta3-ctheta2)/U; 
    Rdeg = (-phi/2)*(tand(b2) + tand(b3));
end

% ---- rotor loss / aspect ratio model ----
function [bladeARr, Prtot, To2, To3, rho3, Nr, M2rel, zetar, br] = ...
    rotor_blade_aspect_ratio(U, sigmar, sigmarz, rm, Reor, givens, R, Cp, ...
                             To1, Po1, Po2, gamma, M2, c2, w2, w3, h2, c3, psi, b2, b3)

    iterations = 2000;
    To2 = To1;
    T2  = To2 - (c2^2)/(2*Cp);
    To3 = To2 + (psi * U^2)/Cp;

    % transform to rotor frame
    M2rel   = (M2/c2)*w2;
    Po2rel  = Po2*((1+((gamma-1)/2)*M2rel^2)^(gamma/(gamma-1))) / ...
                  ((1+((gamma-1)/2)*M2^2)^(gamma/(gamma-1)));
    To2rel  = T2*(1+((gamma-1)/2)*M2rel^2);
    PoIterable = Po2rel;

    T3    = To3 - (c3^2)/(2*Cp);
    M3    = c3/sqrt(gamma*R*T3);
    M3rel = (M3/c3)*w3;

    for i = 1:iterations
        T3    = To3 * (1 + ((gamma-1)/2)*(M3^2))^-1;
        P3    = PoIterable(i) * (1 + ((gamma-1)/2)*(M3rel^2))^(-gamma/(gamma-1));
        rho3  = P3/(R*T3);
        mur   = givens.mu_comb_ref * (T3/givens.mu_comb_tref)^0.7;
        nur   = mur/rho3;

        or = Reor * nur / w3;
        sr = or / max(cosd(b3), 1e-6);
        Nr = getIntegerBladeCount(2*pi*rm/sr);
        sr = (2*pi*rm)/Nr;
        br = sigmar*sr;                 % chord
        bladeARr = h2/br;

        % rotor loss correlation
        bzr    = sigmarz * sr;
        Crotor = 0.975 + 0.075*(bzr/h2);

        zetastarr = 1.04 + 0.06*((b2 + b3)/100)^2;
        Dhr       = ((2*sr*h2*cosd(b3))/(sr*cosd(b3)+h2));
        Rer       = (rho3 * w3 * Dhr) / mur;
        zetar     = (zetastarr*Crotor - 1)*((10^5)/Rer)^0.25;

        num    = (1 - ((w3^2)/(2*Cp*To2rel)) * (1/(1-zetar)));
        den    = (1 - ((w3^2)/(2*Cp*To2rel)));
        Po3rel = Po2rel - ((1 - (num/den)^(gamma/(gamma-1)))*Po2rel);

        PoIterable(i+1) = Po3rel;
    end

    % back to absolute frame
    Po3   = Po3rel*((1+((gamma-1)/2)*M3^2)^(gamma/(gamma-1))) / ...
                   ((1+((gamma-1)/2)*M3rel^2)^(gamma/(gamma-1)));
    Prtot = Po1/Po3;
end

% ---- nozzle loss / aspect ratio model ----
function [bladeARn, Po2, Nn, zetan] = ...
    nozzle_blade_aspect_ratio(rm, Reon, givens, R, Cp, To1, Po1, gamma, ...
                              M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1)

    iterations = 2000;
    PoIterable = Po1;

    for i = 1:iterations
        T2   = To1 * (1 + ((gamma-1)/2)*(M2^2))^-1;
        P2   = PoIterable(i) * (1 + ((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
        rho2 = P2/(R*T2);

        mun = givens.mu_comb_ref * (T2/givens.mu_comb_tref)^0.7;
        nun = mun/rho2;

        c2  = sqrt(ctheta2^2 + cz^2);
        on  = Reon * nun / c2;
        sn  = on / max(cosd(a2), 1e-6);
        Nn  = getIntegerBladeCount(2*pi*rm/sn);
        sn  = (2*pi*rm)/Nn;
        bn  = sigman*sn;
        bladeARn = h1/bn;

        bzn     = sigmanz * sn;
        Cnozzle = 0.993 + 0.021*(bzn/h1);

        zetastarn = 1.04 + 0.06*((a1 + a2)/100)^2;
        Dhn       = ((2*sn*h1*cosd(a2))/(sn*cosd(a2)+h1));
        Ren       = (rho2 * c2 * Dhn) / mun;
        zetan     = (zetastarn*Cnozzle - 1)*((10^5)/Ren)^0.25;

        num = (1 - ((c2^2)/(2*Cp*To1)) * (1/(1-zetan)));
        den = (1 - ((c2^2)/(2*Cp*To1)));
        Po2 = Po1 - ((1 - (num/den)^(gamma/(gamma-1)))*Po1);

        PoIterable(i+1) = Po2;
    end
end

% ---- blade count rounding ----
function Nrounded = getIntegerBladeCount(Nunrounded)
    if ~isreal(Nunrounded) || ~isfinite(Nunrounded) || Nunrounded < 3
        Nrounded = 3;
        return;
    end
    Nrounded = floor(Nunrounded);
    if mod(Nrounded,2) == 0
        Nrounded = Nrounded + 1;
    end
end

% ---- Rotor blade stress calculation ----
function [ANsq, sigmaCent, sigmaBend] = ...
    bladeStresses(rtr, rh, turbine, N, cz, sigmar, br, Cp, To1, To3, P1, P2)
    % AN^2
    Az   = 2*pi*(rtr^2 - rh^2);
    ANsq = Az * N^2;

    % centrifugal stress 
    term1     = (rtr + 2*rh)/6;
    term2     = (2*rtr + rh)/6;
    sigmaCent = ((N*pi/30)^2) * (rtr - rh) * (term1 + turbine.taperRatio * term2);

    % bending stress 
    term1b     = (cz / (rtr * N * pi/30));
    term2b     = (Cp * abs(To3 - To1)) / (Cp * To1);
    term3b     = 1/(2*sigmar);
    term4b     = (rtr / (br * turbine.thicknessMax));
    p_avg      = (P1 + P2)/2;
    sigmaBend  = p_avg * term1b * term2b * term3b * term4b;
end

% ---- post-processing and output for first stage ----
function turbine_stage = postProcess(x, R, To1, Po1, a1, Cp, ...
                                     gamma, N, mdot, zweif, givens, ...
                                     etaTR, Wout, M2, turbine)
    [a3, rh, Reon, Reor, cz] = unpack(x);

    [c1, ctheta1, c2, ctheta2, a2, ...
     w2, wtheta2, b2, w3, wtheta3, b3, ...
     c3, ctheta3, M1, rtn, rmn, U, rho1, P1] = ...
        velocityTriangleAnalysis(cz, a1, a3, To1, Po1, M2, Cp, ...
                                 gamma, R, mdot, rh, N);

    [psi, phi, Rdeg] = psiPhiR(cz, U, ctheta3, ctheta2, b2, b3);

    [sigmanz, sigmarz, sigman, sigmar] = ...
        sigmaZweif(w3, To1, U, Cp, R, cz, ctheta1, ctheta2, ...
                   wtheta2, wtheta3, a2, b3, M2, psi, c3, zweif, gamma);

    h1 = rtn - rh;
    [bladeARn, Po2, Nn, zetan] = ...
        nozzle_blade_aspect_ratio(rmn, Reon, givens, R, Cp, To1, Po1, ...
                                  gamma, M2, sigman, sigmanz, a1, a2, ...
                                  ctheta2, cz, h1);

    To2 = To1;
    [rtr, rmr] = rotorRadii(cz, rh, mdot, To2, Po2, M2, gamma, R);

    % Static P2 from Po2 and M2
    P2 = Po2 * (1 + ((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));

    h2 = rtr - rh;
    [bladeARr, Prtot, To2, To3, rho3, Nr, M2rel, zetar, br] = ...
        rotor_blade_aspect_ratio(U, sigmar, sigmarz, rmr, Reor, givens, ...
                                 R, Cp, To1, Po1, Po2, gamma, M2, ...
                                 c2, w2, w3, h2, c3, psi, b2, b3);

    % Efficiency
    num   = psi*(gamma-1)*((U^2)/(gamma*R*To1));
    den   = Prtot^((1-gamma)/gamma) - 1;
    etaST = num/den;

    Po3 = Po1/Prtot;
    num_poly = log(To1/To3);
    den_poly = log(1 + (1/etaST)*(To1/To3 - 1));
    etaPoly  = num_poly/den_poly;
    rho3_rho1 = rho3/rho1;

    % Blade stresses
    [ANsq, sigmaCent, sigmaBend] = bladeStresses(rtr, rh, turbine, N, ...
                                                 cz, sigmar, br, Cp, ...
                                                 To1, To3, P1, P2);

    % Ratios
    Po2_Po1 = Po2/Po1;
    Po3_Po2 = (1/Prtot) * (Po2_Po1)^(-1);
    AN2 = pi*(rtr^2 - rh^2) * (N^2);
    chi2 = a2; 
    chi3 = b3 - cartersDeviation((b3-b2), sigmar); 

    % Actual work this stage
    Wdot = calculateWorkDone(etaTR, mdot, U, ctheta3, ctheta2);

    % Excel output
    output = {
        '',   '°R', 'U (m/s)', 'rm (m)', 'phi', 'psi', 'To3 (K)', 'Po3 (kPa)', 'etaST (%)', 'etaPoly (%)', 'rho3/rho1', '', '';
        '',    Rdeg,   U,        rmn,      phi,  psi,  To3,       (Po3/1000), (etaST*100), (etaPoly*100), rho3_rho1,  '', '';
        '',   'chi2 (deg)', 'a2 (deg)', 'sigmaN', 'Nn', 'Reon', 'rhn', 'rtn', 'bladeARn', 'zetan', 'Po2/Po1', 'M1', '';
        '',    chi2,         a2,         sigman,   Nn,  Reon,   rh,   rtn,   bladeARn,   zetan,   Po2_Po1,   M1,  '';
        'b3 (deg)', 'chi3 (deg)', 'a3 (deg)', 'sigmaR', 'Nr', 'Reor', 'rhr', 'rtr', 'bladeARr', 'zetar', 'Po3/Po2', 'M2rel', 'AN^2 (m^2 RPM^2)';
        b3,          chi3,        a3,         sigmar,    Nr,   Reor,  rh,    rtr,   bladeARr,   zetar,   Po3_Po2,   M2rel,   AN2
        };
    %writecell(output, 'turbineData_firstStage.xlsx');

    % write struct
    turbine_stage.Pr        = Po3/Po1;             
    turbine_stage.Po3       = Po3;                 
    turbine_stage.To3       = To3;                
    turbine_stage.Po2       = Po2;                 
    turbine_stage.To2       = To2;                 
    turbine_stage.rh        = rh;                
    turbine_stage.rtr       = rtr;                
    turbine_stage.rtn       = rtn;                
    turbine_stage.rmn       = rmn;                
    turbine_stage.rmr       = rmr;                
    turbine_stage.psi       = psi;
    turbine_stage.phi       = phi;
    turbine_stage.eta       = etaST;               
    turbine_stage.eta_poly  = etaPoly;             
    turbine_stage.a3        = a3;                  
    turbine_stage.b3        = b3;                  
    turbine_stage.sigmaR    = sigmar;              
    turbine_stage.sigmaN    = sigman;
    turbine_stage.Nr        = Nr;
    turbine_stage.Nn        = Nn;
    turbine_stage.Reor      = Reor;
    turbine_stage.Reon      = Reon;
    turbine_stage.AN2       = AN2;
    turbine_stage.Po2_Po1   = Po2_Po1;
    turbine_stage.Po3_Po2   = Po3_Po2;
    turbine_stage.a2        = a2;
    turbine_stage.chi2      = chi2;
    turbine_stage.cz        = cz;

    turbine_stage.M2        = M2;          % store fixed M2 = 1 for first stage

    turbine_stage.ANsq      = ANsq;
    turbine_stage.sigmaCent = sigmaCent;
    turbine_stage.sigmaBend = sigmaBend;

    % work info
    turbine_stage.Wdot      = Wdot;
    turbine_stage.Wdot_req  = Wout;

    turbine_stage.bladeARr = bladeARr;
    turbine_stage.bladeARn = bladeARn;
    turbine_stage.chi3 = chi3;
    turbine_stage.a1 = a1;
end

function d = cartersDeviation(flowAngleDelta, sigma)
    m = 1/8;
    d = m*(abs(flowAngleDelta)/sigma);
end
