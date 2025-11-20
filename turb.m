function turbine = turb(a1_turb, Nstages_turb, N, givens, properties, turbine)
    % -------- required inputs --------
    if ~isfield(turbine,'W_total')
        error('turb.m: turbine.W_total is required.');
    end
    if ~isfield(turbine,'Nstages') || turbine.Nstages ~= Nstages_turb
        warning('turb.m: turbine.Nstages not set or inconsistent; using Nstages_turb argument.');
        turbine.Nstages = Nstages_turb;
    end
    Nst = Nstages_turb;
    % -------- linear "work ratio" distribution (internal) --------
    if isfield(turbine,'Wfrac')
        % explicit fractions from user
        Wfrac_use = turbine.Wfrac(:).';
        if numel(Wfrac_use) ~= Nst
            error('turb.m: length(turbine.Wfrac) must equal Nstages_turb.');
        end
        if any(Wfrac_use <= 0) || ~all(isfinite(Wfrac_use))
            error('turb.m: turbine.Wfrac must be positive and finite.');
        end
        Wfrac_use = Wfrac_use / sum(Wfrac_use);
        Wrs_use   = Wfrac_use;
    else
        % Build Wr(k) = Wr1 - d*(k-1) with PRODUCT(Wr) = 1 (product-type condition)
        Wr1 = getfield_def(turbine,'Wr1', 1.2);    % analogous to Pr1_comp
        if Nst == 1
            Wrs_use   = 1.0;
            Wfrac_use = 1.0;
        else
            % --- NEW: handle Wr1 <= 1 by falling back to equal work split ---
            if Wr1 <= 1.0
                % equal work ratios if the first-stage bias is not > 1
                Wrs_use   = ones(1, Nst);
                Wfrac_use = ones(1, Nst) / Nst;
            else
                idx = 0:(Nst-1);
                fun = @(d) prod(Wr1 - d*idx) - 1.0;
                d_low  = 0.0;
                d_high = Wr1/(Nst-1);  % ensures last Wr >= 0 at upper bound
                d      = fzero(fun, [d_low, d_high]);
                Wrs_use = Wr1 - d*idx;
                if any(Wrs_use <= 0)
                    error('turb.m: constructed work ratios Wrs became non-positive; adjust Wr1.');
                end
                Wfrac_use = Wrs_use / sum(Wrs_use);
            end
        end
    end
    % -------- store distribution in turbine struct --------
    turbine.Wrs   = Wrs_use;
    turbine.Wfrac = Wfrac_use;
    tmpl = struct( ...
    'Pr',NaN, 'Po3',NaN, 'To3',NaN, 'Po2',NaN, 'To2',NaN, ...
    'rh',NaN, 'rtr',NaN, 'rtn',NaN, 'rmn',NaN, 'rmr',NaN, ...
    'psi',NaN, 'phi',NaN, 'eta',NaN, 'eta_poly',NaN, ...
    'a3',NaN, 'b3',NaN, 'sigmaR',NaN, 'sigmaN',NaN, ...
    'Nr',NaN, 'Nn',NaN, 'Reor',NaN, 'Reon',NaN, 'AN2',NaN, ...
    'Po2_Po1',NaN, 'Po3_Po2',NaN, 'a2',NaN, 'chi2',NaN, 'cz',NaN, ...
    'ANsq',NaN, 'sigmaCent',NaN, 'sigmaBend',NaN, ...
    'M2',NaN, ...
    'Wdot',NaN, 'Wdot_req',NaN );
    turbineStages(1:Nst+1) = tmpl;
    % -------- seed "stage 0" (inlet to turbine) --------
    turbineStages(1).Pr  = 1.0;
    turbineStages(1).Po3 = properties.Po4;
    turbineStages(1).To3 = properties.To4;
    turbineStages(1).a3  = a1_turb;   % inlet absolute angle
    % -------- loop over stages --------
    eta_sum = 0.0;
    for k = 1:Nst
        turbine.WperStage = turbine.W_total * Wfrac_use(k);
        % Optional: force last-stage a3 (e.g., 0 deg)
        if k == Nst && isfield(turbine,'last_stage_a3')
            turbine.force_a3 = turbine.last_stage_a3;
        else
            if isfield(turbine,'force_a3')
                turbine = rmfield(turbine,'force_a3');
            end
        end
        prev = turbineStages(k);
        if k == 1
            % First stage: cz is a design variable inside first_turb_stage
            stage_k = first_turb_stage(N, givens, properties, turbine, prev);
        else
            % Remaining stages: use cz from previous stage
            stage_k = turb_stage(N, turbineStages(k).cz, givens, properties, turbine, prev);
        end
        turbineStages(k+1) = stage_k;
        if isfield(stage_k,'eta') && isfinite(stage_k.eta)
            eta_sum = eta_sum + stage_k.eta;
        end
    end
        turbine.allStages = turbineStages(2:end); % wipe out first dummy stage

    % keep the stage-average efficiency if you still want it
    %stageAvgEfficiency = eta_sum / Nst;

    % overall total-pressure ratio Po_out/Po_in
    Pr_total = 1.0;
    for k = 2:(Nst+1)
        if isfield(turbineStages(k),'Pr') && isfinite(turbineStages(k).Pr)
            Pr_total = Pr_total * turbineStages(k).Pr;
        end
    end
    turbine.Pr = Pr_total;

    % overall polytropic efficiency for the entire turbine
    Po_in  = properties.Po4;                 % turbine inlet total pressure
    Po_out = Po_in * Pr_total;               % Pr_tot = Po_out/Po_in

    To_in  = properties.To4;                 % turbine inlet total temperature
    To_out = turbineStages(end).To3;         % last stage outlet total T

    gamma  = properties.gamma_comb;

    % Prtot_total = Po_in / Po_out
    Prtot_total = Po_in / Po_out;

    % etaST_total from overall T and P ratios
    num_TT = 1 - To_out/To_in;
    den_TT = 1 - Prtot_total^((1-gamma)/gamma);
    etaST_total = num_TT / den_TT;

    % num_poly = log(To1/To3);
    % den_poly = log(1 + (1/etaST)*(To1/To3 - 1));
    num_poly  = log(To_in/To_out);
    den_poly  = log(1 + (1/etaST_total)*(To_in/To_out - 1));
    eta_poly_total = num_poly / den_poly;

    % store polytropic efficiency as the "average" efficiency reported
    turbine.polyEfficiency      = eta_poly_total;
    %turbine.stageAvgEfficiency = stageAvgEfficiency;

    turbine.Nstages = Nst;
    turbine.W_total = turbine.W_total;
    turbine.Wfrac   = Wfrac_use;
    turbine.Wrs     = Wrs_use;
end

function val = getfield_def(s, name, default)
    if isstruct(s) && isfield(s,name)
        val = s.(name);
    else
        val = default;
    end
end
 