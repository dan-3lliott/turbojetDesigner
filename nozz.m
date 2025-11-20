function conv_nozzle = nozz(givens, properties, nozzle, combustor)

% givens
mdot = givens.mdot_air;
R_n = properties.R_comb;
cp_n = nozzle.cp_n;
To5 = nozzle.To5;
Po5 = nozzle.Po5;
g_n = nozzle.gamma;
eta_n = nozzle.eta_n;
Pa = properties.Pa;
Minf = givens.M;
Ta = properties.Ta;
g_d = 1.4;
R_d = 287;
To4 = properties.To4;
To3 = properties.To3;
cp_c = properties.Cp_comb;

f = combustor.f;
Pe = Po5 * (2/(g_n + 1))^(g_n/(g_n - 1)); % Pa
gamma_exp = (g_n - 1)/g_n;
Te = To5*(1 - eta_n*(1 - (Pe/Po5)^gamma_exp)); % K
%Te = To5*(1 - (1 - (Pe/Po5)^gamma_exp)); % K
ue = sqrt(2*cp_n*(To5 - Te)); % m/s

% Choked flow at throat
choked_gammaexp = -(g_n + 1)/(2*(g_n - 1));
choked_factor = (Po5/sqrt(To5))*sqrt(g_n/R_n)*((g_n + 1)/2)^choked_gammaexp;
Ath = mdot*(1 + f)*(choked_factor)^-1; % m^2
gamma_factor = ((g_n + 1)/2)^((g_n + 1)/(2*(g_n - 1)));

Ae = (mdot * (1 + f) * sqrt(R_n*To5) / Po5) * (1/sqrt(g_n)) * gamma_factor; % m^2 (should be same as Ath but just a check)
u = Minf*sqrt(g_d*R_d*Ta); % m/s
Me = ue/sqrt(g_n*R_n*Te);
thrust = mdot*((1 + f)*ue - u) + (Pe - Pa)*Ae; % N

conv_nozzle = nozzle;
conv_nozzle.thrust = thrust;
conv_nozzle.Ath      = Ath;          
conv_nozzle.Ae       = Ae;          
conv_nozzle.Pe       = Pe;           
conv_nozzle.Te       = Te;           
conv_nozzle.ue       = ue;           
conv_nozzle.Me       = Me;           


end