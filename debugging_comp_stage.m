clear all
clc

%debugging variables
a1_comp = 35; %deg
M2 = 0.5;
Nstages_comp = 1;
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
properties.gamma_comb = 1.28;
properties.MW_comb = 28.7; %g/mol
properties.R_comb = (properties.Rbar/properties.MW_comb);
properties.To4 = 1680; %K

%assumptions - diffuser
diffuser.Rd = 0.985; %pressure recovery coefficient
diffuser.M2max = 0.45;

%assumptions - compressor
compressor.Dfmax = 0.6;
compressor.Rebmin = 5e5;
compressor.Mrelmax = 0.79;
compressor.bladeARmin = 1.1;
compressor.bladeARmax = 1.6;
compressor.wbarfactor = 2.3;
compressor.taperRatio = 0.8;
compressor.thicknessMax = 0.1; %proportion of chord b
compressor.aFinal = 0; %deg
compressor.a3StageMax = 45; %deg
compressor.a3a1offset = 5; %deg
compressor.Prreq = 10.9;
compressor.Pr = 1; %initial value only

%assumptions - combustor
combustor.etaComb = 0.99;
combustor.dHR = 43.8e6; %J/kg
combustor.tr = 2.1e-3; %sec
combustor.mSnoutRatio = 0.65;
combustor.dPoLoverPo3 = 0.041;

%assumptions - turbine
turbine.taperratio = 0.5;
turbine.thicknessmax = 0.2; %proportion of chord b
turbine.zweif = 0.8;
turbine.Retmin = 1e5;
turbine.Retmax = 5e5;
turbine.etaShaft = 0.995;
turbine.aFinal = 0; %deg

%assumptions - nozzle
nozzle.gamma = 1.35;
nozzle.etaAdia = 0.95;

%calculate stagnation conditions at inlet
properties.To1 = properties.Ta*(1+((properties.gamma_air-1)/2)*givens.M^2);
properties.To2 = properties.To1; %assuming adiabatic diffuser
properties.Po1 = properties.Pa*(1+((properties.gamma_air-1)/2)*givens.M^2)^(properties.gamma_air/(properties.gamma_air-1));
properties.Po2 = properties.Po1*diffuser.Rd;

%determine cz given choice of M2 and a1 - relates to diffuser and IGV
T2 = properties.To2*(1+((properties.gamma_air-1)/2)*M2^2)^-1;
%compressor.cz = M2*sqrt(properties.gamma_air*properties.R_air*T2)*cosd(a1_comp);
%prepare "zeroth" stage structure - inlet guide vane. this needs to
%change
compressorStages(1).Pr = 1;
compressorStages(1).Po3 = properties.Po2;
compressorStages(1).To3 = properties.To2;
compressorStages(1).rhr = 1;
compressorStages(1).rhs = 1;
compressorStages(1).rmr = 1;
compressorStages(1).rms = 1;
compressorStages(1).psi = 1;
compressorStages(1).phi = 1;
compressorStages(1).eta = 1;
compressorStages(1).M1rel = 1;
compressorStages(1).a3 = a1_comp;
compressorStages(1).chi2r = 1;
compressorStages(1).chi3 = 1;
compressorStages(1).Dfr = 1;
compressorStages(1).Dfs = 1;
compressorStages(1).sigmaR = 1;
compressorStages(1).sigmaS = 1;
compressorStages(1).Po2 = 1;

amax = 28; %deg
[compressor_stage, rt_comp, N_rpm, a1_comp, cz_comp] = first_comp_stage(1.53, diffuser, givens, properties, compressor, compressorStages(1), amax)
%calculate other necessary values