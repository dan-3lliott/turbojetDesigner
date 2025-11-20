clear; clc;

% ====== givens/properties ======
givens.M = 0.80;
givens.mdot_air     = 5.02;     % [kg/s] air mass flow
givens.mu_comb_ref  = 5.6e-5;   % [kg/m·s] (for turbine/OGV loss calcs)
givens.mu_comb_tref = 1550;     % [K]     (for turbine/OGV loss calcs)

combustor.f = 0.03;
compressor.cz = 166;

properties.gamma_comb = 1.28;   % [-] hot-gas gamma
properties.MW_comb    = 28.7;   % [g/mol]
properties.Rbar       = 8314;   % [J/(kmol·K)]
properties.R_comb     = properties.Rbar / properties.MW_comb;              % [J/(kg·K)]
properties.Cp_comb    = (properties.gamma_comb*properties.R_comb) ...
                        /(properties.gamma_comb - 1);                       % [J/(kg·K)]

% Ambient / flight (for nozzle/thrust)
properties.Ta   = 249;          % [K]
properties.Minf = 0.80;         % [-]
properties.Pa   = 45.1e3;       % [Pa] (set as needed for thrust; replace with actual ambient)

% Combustor data you may need elsewhere
properties.eta_c = 0.99;        % combustor efficiency [-]
properties.hv    = 43.8e6;      % fuel lower heating value [J/kg]

% --- Turbine inlet totals
properties.To4 = 1680;                               % [K] total T at turbine inlet (station 4)
properties.To3 = 589.4;                                % [K] compressor exit (for f calc elsewhere if needed)
properties.Po4 = 6.5835e5;               % [Pa] example Po4

% ====== Turbine operating guesses ======
N_rpm   = 2.7992e4;       % [rpm] shaft speed
rt_turb = 0.25;       % [m]   tip radius at turbine inlet (kept for API)
a1_turb = 0;         % [deg] absolute inlet flow angle to stage 1

% ====== Turbine setup ======
turbine.Nstages      = 2;         % number of stages
turbine.etaShaft     = 0.995;     % shaft mechanical efficiency
turbine.zweif        = 0.80;      % Zweifel target
turbine.Retmin       = 1e5;       % Reynolds lower bound
turbine.Retmax       = 5e5;       % Reynolds upper bound
%turbine.W_total      = 1.2e6;     % [W] total shaft work to distribute
turbine.W_total      = 1557440.56464342;     % [W] total shaft work to distribute

turbine.taperRatio   = 0.5;       % used in blade stress model
turbine.thicknessMax = 0.2;       % used in blade stress model (lowercase "max")

% linear "work ratio" parameter for stage 1, Wr1_turb.
% turb.m will build a sequence Wr(k) = Wr1 - d*(k-1)
% with product_k Wr(k) = 1, then normalize to Wfrac.
turbine.Wr1          = 1.05;       % dimensionless "work ratio" for stage 1

% Run turbine (turb.m will build linear Wfrac internally) ======
% cz is now chosen as a design variable in first_turb_stage, not passed in here
T = turb(a1_turb, turbine.Nstages, N_rpm, givens, properties, turbine);

nStages = numel(T.allStages)

% ====== Outlet Guide Vane (OGV) — zero swirl, stator loss, constant hub ====
lastStage = T.allStages(end);           % turbine exit (rotor outlet)
OGV       = ogv_turb(givens, properties, turbine, lastStage);

% ====== Nozzle inlet totals (station 5) from OGV ======
nozzle.To5   = OGV.To_out;
nozzle.Po5   = OGV.Po_out;
nozzle.gamma = 1.25;                                       % set per your model
nozzle.eta_n = 0.95;                                       % nozzle efficiency
nozzle.cp_n  = (nozzle.gamma*properties.R_comb)/(nozzle.gamma - 1);

% ====== Run nozzle model (nozz.m) ======
conv_nozzle = nozz(givens, properties, nozzle, combustor);

% ====== Report results ======
fprintf('\n=== Turbine Results ===\n');
fprintf('Stages                 : %d\n', T.Nstages);
fprintf('Total shaft work (W)   : %.3e\n', T.W_total);
fprintf('Overall polytropic efficiency   : %.4f\n', T.polyEfficiency);
fprintf('Overall Po_out/Po_in   : %.4f\n', T.Pr);

% Work fractions from turb.m (linear internal distribution)
if isfield(T, 'Wfrac')
    Wfrac = T.Wfrac(:).';
    src   = 'linear (internal in turb.m)';
else
    Wfrac = nan(1, turbine.Nstages);
    src   = 'not available';
end
fprintf('\nPer-stage work fractions (%s, sum=%.6f):\n', src, sum(Wfrac));
for k = 1:numel(Wfrac)
    fprintf('  Stage %d: %.4f (%.1f%%)\n', k, Wfrac(k), 100*Wfrac(k));
end

% Show turbine-exit (pre-OGV) and OGV-exit (post-OGV) flow angles
if isfield(lastStage,'a3') && isfinite(lastStage.a3)
    fprintf('\nTurbine exit a3 (deg)  : %.3f\n', lastStage.a3);
else
    fprintf('\nTurbine exit a3 (deg)  : (not available)\n');
end
fprintf('OGV exit a_out (deg)   : %.3f\n', OGV.a_out);
fprintf('OGV Po_out/Po_in       : %.4f\n', OGV.Po_out_Po_in);

fprintf('\nPer-stage efficiencies:\n');
for k = 1:(T.Nstages)   % skip the seed (index 1)
    eta_k = T.allStages(k).eta;
    fprintf('  Stage %d: %.4f\n', k, eta_k);
end

fprintf('\nPer-stage M2:\n');
for k = 1:(T.Nstages)   % skip the seed (index 1)
    M2_k = T.allStages(k).M2;
    fprintf('  Stage %d: %.4f\n', k, M2_k);
end

% New: per-stage phi and cz (now cz is DV in stage 1, propagated)
fprintf('\nPer-stage phi and cz:\n');
for k = 1:(T.Nstages)
    st = T.allStages(k);
    fprintf('  Stage %d: phi = %.4f,  cz = %.3f m/s\n', ...
            k, st.phi, st.cz);
end

% ====== Overall total-pressure ratio across turbine ======
fprintf('\n=== Turbine Total-Pressure Ratio ===\n');
fprintf('Overall Po_out/Po_in (T.Pr) : %.5f\n', T.Pr);

% Compute the product of all individual stage Pr values for comparison
Pr_product = 1.0;
for k = 1:(T.Nstages)
    Pr_k = T.allStages(k).Pr;
    fprintf('  Stage %d Pr = %.5f\n', k, Pr_k);
    Pr_product = Pr_product * Pr_k;
end
fprintf('Product of stage Pr values  : %.5f\n', Pr_product);

% Show % loss
fprintf('Total pressure *loss* (1 - T.Pr) : %.3f (%.2f%%)\n', ...
        1 - T.Pr, 100*(1 - T.Pr));

% ====== Shaft Mechanical Efficiency ======
if isfield(turbine, 'etaShaft')
    fprintf('\n=== Shaft Mechanical Efficiency ===\n');
    fprintf('Specified η_shaft: %.4f  (%.2f%%)\n', turbine.etaShaft, 100*turbine.etaShaft);
else
    fprintf('\n(No etaShaft field found in turbine struct)\n');
end

% ====== Turbine geometric / loading requirements & actuals ======
fprintf('\n=== Turbine Geometric / Loading Requirements ===\n');
fprintf('Taper ratio (target)                : %.3f\n', turbine.taperRatio);
fprintf('Max thickness-to-chord (t/c target) : %.3f\n', turbine.thicknessMax);
fprintf('Zweifel loading |psi_z| target      : %.3f\n', turbine.zweif);
fprintf('Throat Reynolds number target range : [%.2e, %.2e]\n', ...
        turbine.Retmin, turbine.Retmax);

% Per-stage actual throat Reynolds numbers (from turb_stage / first_turb_stage)
fprintf('\nPer-stage throat Reynolds numbers (actual):\n');
for k = 1:(T.Nstages)   % skip the seed (index 1)
    st = T.allStages(k);

    Re_nozzle = NaN;
    Re_rotor  = NaN;
    if isfield(st, 'Reon'), Re_nozzle = st.Reon; end
    if isfield(st, 'Reor'), Re_rotor  = st.Reor; end

    fprintf('  Stage %d: Re_nozzle = %.3e,  Re_rotor = %.3e\n', ...
            k, Re_nozzle, Re_rotor);
end
fprintf('\n');

% Per-stage actual throat Reynolds numbers (from turb_stage / first_turb_stage)
fprintf('\nPer-stage tip radii (actual):\n');
for k = 1:(T.Nstages)   % skip the seed (index 1)
    st = T.allStages(k);

    rt_nozz = NaN;
    rt_rotor  = NaN;
    if isfield(st, 'rtn'), rt_nozz = st.rtn; end
    if isfield(st, 'rtr'), rt_rotor  = st.rtr; end

    fprintf('  Stage %d: rtn = %.3e,  rtr = %.3e\n', ...
            k, rt_nozz, rt_rotor);
end
fprintf('\n');

fprintf('\nPer-stage mean radii (actual):\n');
for k = 1:(T.Nstages)   % skip the seed (index 1)
    st = T.allStages(k);

    rm_nozz = NaN;
    rm_rotor  = NaN;
    if isfield(st, 'rmn'), rm_nozz = st.rmn; end
    if isfield(st, 'rmr'), rm_rotor  = st.rmr; end

    fprintf('  Stage %d: rmn = %.3e,  rmr = %.3e\n', ...
            k, rm_nozz, rm_rotor);
end
fprintf('\n');

sumW = 0;
for k = 1:(T.Nstages)
    st = T.allStages(k);
    if isfield(st,'Wdot')
        sumW = sumW + st.Wdot;
        fprintf('Stage %d: Wreq = %.3e W, Wdot = %.3e W\n', ...
                k, st.Wdot_req, st.Wdot);
    end
end
fprintf('Sum of stage work: %.3e W (target %.3e W)\n', sumW, T.W_total);


% Quick nozzle inlet confirmation
fprintf('\n=== Nozzle Inlet (from OGV) ===\n');
fprintf('To5 (K)                : %.2f\n', nozzle.To5);
fprintf('Po5 (kPa)              : %.2f\n', nozzle.Po5/1e3);
fprintf('gamma_n, eta_n         : %.3f, %.3f\n', nozzle.gamma, nozzle.eta_n);

% ====== Nozzle performance outputs ======
fprintf('\n=== Nozzle Performance (nozz.m) ===\n');
fprintf('Thrust (N)             : %.3f\n', conv_nozzle.thrust);
fprintf('Throat area Ath (m^2)  : %.6f\n', conv_nozzle.Ath);
fprintf('Exit area   Ae  (m^2)  : %.6f\n', conv_nozzle.Ae);
fprintf('Exit pressure Pe (Pa)  : %.2f\n', conv_nozzle.Pe);
fprintf('Exit temperature Te(K) : %.2f\n', conv_nozzle.Te);
fprintf('Exit velocity ue (m/s) : %.2f\n', conv_nozzle.ue);
fprintf('Exit Mach Me           : %.3f\n', conv_nozzle.Me);
fprintf('\n');
