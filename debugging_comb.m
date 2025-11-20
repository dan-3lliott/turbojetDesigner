clear all
clc

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
compressor.Pr_igv = 0.9920; %stagnation pressure ratio for inlet guide vane
compressor.Pr_ogv = 0.9941; %same thing for ogv
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
combustor.Ar_prediff = 1.3;
combustor.theta = 10; %deg
combustor.rmr = 0.95; %mean radius ratio from compressor -> combustor, manually select
combustor.ALr = 0.375; %liner area ratio, manually select

%assumptions - turbine
turbine.taperratio = 0.5;
turbine.thicknessmax = 0.2; %proportion of chord b
turbine.zweif = 0.8;
turbine.Retmin = 1e5;
turbine.Retmax = 5e5;
turbine.etaShaft = 0.995;
turbine.aFinal = 0; %deg
turbine.Nstages = 3; %PICK THIS!!!!!

%assumptions - nozzle
nozzle.gamma = 1.35;
nozzle.etaAdia = 0.95;

%calculate stagnation conditions at inlet
properties.To1 = properties.Ta*(1+((properties.gamma_air-1)/2)*givens.M^2);
properties.To2 = properties.To1; %assuming adiabatic diffuser
properties.Po1 = properties.Pa*(1+((properties.gamma_air-1)/2)*givens.M^2)^(properties.gamma_air/(properties.gamma_air-1));
properties.Po2 = properties.Po1*diffuser.Rd;

%calculate diffuser velocity
properties.T2 = properties.To2*(1+((properties.gamma_air-1)/2)*diffuser.M2max^2)^-1;
diffuser.cmax = diffuser.M2max * sqrt(properties.gamma_air*properties.R_air*properties.T2)*cosd(compressor.ai);

%DEBUG INPUTS MANUALLY SET
properties.Po3 = 738095; %Pa
properties.To3 = 590.67; %K
compressor.rt = 0.1304; %m
compressor.stages(1).rhs = 0.1197; %m
compressor.stages(1).M1rel = 0.3;

combustor_out = comb(givens, properties, combustor, compressor)
disp('Po4/Po3:');
disp(combustor_out.Po4/properties.Po3);
disp('L:');
disp(combustor_out.L);