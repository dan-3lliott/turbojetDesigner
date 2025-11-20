function od = offDesign(givens, properties, diffuser, compressor, combustor, turbine, nozzle)
    %=====ASSUMPTIONS AND QUANTITIES=====

    %define needed quantities
    To4_min = 1500; %K
    To4_max = 1700; %K
    To4_sweep = To4_min:To4_max;
    Tref = 288.15; %K
    Pref = 101325; %Pa

    %define offdesign object
    od.To4 = To4_sweep;
    od.Po2 = properties.Po2;
    od.Ta = 216.7; %K
    od.Pa = 17.9e3; %Pa
    od.M = 0.92;

    %calculate new off design properties
    od.To1 = od.Ta*(1+((properties.gamma_air-1)/2)*od.M^2);
    od.To2 = od.To1; %assuming adiabatic diffuser
    od.Po1 = od.Pa*(1+((properties.gamma_air-1)/2)*od.M^2)^(properties.gamma_air/(properties.gamma_air-1));
    od.Po2 = od.Po1*diffuser.Rd;
    od.u = od.M*sqrt(properties.gamma_air*properties.R_air*od.Ta);

    %=====CALCULATIONS=====

    %off-design compressor pressure ratio
    od.To3 = od.To2.*(1 + ((od.To4./od.To2)./(properties.To4/properties.To2)).*(compressor.To3/properties.To2 - 1)); %assumes that 1+f is essentially the same
    od.Pr_comp = (od.To3./od.To2).^((properties.gamma_air*compressor.eta)/(properties.gamma_air-1));
    
    %off-design corrected mass flow
    mc2_d = givens.mdot_air*sqrt(properties.To2/Tref)/(properties.Po2/Pref);
    term1 = sqrt(properties.To4/properties.To2)./sqrt(od.To4./od.To2);
    od.mc2 = mc2_d.*(od.Pr_comp./compressor.Pr).*term1;

    %off-design mass flow
    od.m2 = od.mc2.*(od.Po2/Pref)./sqrt(od.To2/Tref);

    %on-design corrected rpm
    Nc_d = compressor.N/sqrt(properties.To2/Tref);

    %off-design corrected rpm
    od.Nc = Nc_d*sqrt(od.To4./od.To2)/sqrt(properties.To4/properties.To2);

    %off-design rpm
    od.N = od.Nc.*sqrt(od.To2./Tref);

    %off-design To5, assuming To5/To4 ratio is ~= const
    od.To5 = od.To4 .* (turbine.allStages(end).To3/properties.To4);

    %diffuser inlet mach
    M1 = sym('M1');
    od.M1 = zeros(1,length(od.m2));
    for i = 1:length(od.m2)
        symexp = (od.m2(i) == ((diffuser.A1*properties.Po1)/sqrt(properties.To1))*sqrt(properties.gamma_air/properties.R_air)*(M1*(1+((properties.gamma_air-1)/2)*M1^2)^(-(properties.gamma_air+1)/(2*(properties.gamma_air-1)))));
        od.M1(i) = vpasolve(symexp, M1, [0 1]); %should be subsonic
    end

    %diffuser exit mach
    M2 = sym('M2');
    od.M2 = zeros(1,length(od.m2));
    for i = 1:length(od.m2)
        symexp = (od.m2(i) == ((diffuser.A2*properties.Po2)/sqrt(properties.To2))*sqrt(properties.gamma_air/properties.R_air)*(M2*(1+((properties.gamma_air-1)/2)*M2^2)^(-(properties.gamma_air+1)/(2*(properties.gamma_air-1)))));
        od.M2(i) = vpasolve(symexp, M2, [0 1]); %should be subsonic
    end
    
    %off-design f - this is an approximation for now calculating nozzle
    %thrust
    od.f = properties.Cp_comb.*(od.To4 - od.To3)./(combustor.etaComb*combustor.dhR - properties.Cp_comb.*od.To4);

    %prt - this just works out to give const prt which we had assumed
    %earlier
    od.Pr_turb = (od.To5./od.To4).^(properties.gamma_comb/(turbine.polyEfficiency*(properties.gamma_comb-1)));

    %off-design Po5
    od.Po5 = (od.Po2 .* od.Pr_comp .* (combustor.Po4/compressor.Po3) .* od.Pr_turb);

    %thrust
    od.Pe = od.Po5.*(1+((properties.gamma_comb-1)/2))^(-properties.gamma_comb/(properties.gamma_comb-1));
    od.Te = od.To5 .*(1 - nozzle.eta_n*(1 - (od.Pe./od.Po5).^((properties.gamma_comb-1)/properties.gamma_comb)));
    od.ue = sqrt(2*properties.Cp_comb.*(od.To5 - od.Te));
    od.Me = od.ue./sqrt(properties.gamma_comb*properties.R_comb.*od.Te);
    od.thrust = od.m2.*((1 + od.f).*od.ue - od.u) + (od.Pe - od.Pa)*nozzle.Ae;

    %st
    od.st = od.thrust./od.m2;

    %sfc
    od.sfc = od.f./od.st;

    %=====OUTPUTS=====

    %generate plots
    figure
    subplot(2,2,1);
    plot(od.To4, od.m2);
    xlabel('To4 (K)');
    ylabel('$\dot{m}_{air}$ (kg/s)', 'Interpreter','latex');
    subplot(2,2,2);
    plot(od.To4, od.N);
    xlabel('To4 (K)');
    ylabel('RPM');
    subplot(2,2,3);
    plot(od.To4, od.Pr_comp);
    xlabel('To4 (K)');
    ylabel('Pr_{c}');
    subplot(2,2,4);
    plot(od.To4, od.thrust);
    xlabel('To4 (K)');
    ylabel('\tau (N)');
    sgtitle('Off-Design Performance');
end

%THINGS I NEED