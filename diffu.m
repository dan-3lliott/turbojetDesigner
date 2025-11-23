function diffuser = diffu(diffuser, givens, properties, compressor)
    %define inlet mach as freestream mach
    diffuser.M1 = givens.M;
    %define diffuser exit area
    diffuser.A2 = pi*(compressor.rt^2 - compressor.stages(1).rhs^2); %must equal IGV area
    %define M2
    diffuser.M2 = compressor.M2_diff;
    %define Po2
    diffuser.Po2 = properties.Po2; %already calculated in engine.m using given Rd
    %calculate diffuser inlet area
    term1 = (diffuser.M2/diffuser.M1);
    num = (1+((properties.gamma_air-1)/2)*diffuser.M1^2);
    den = (1+((properties.gamma_air-1)/2)*diffuser.M2^2);
    term2 = (num/den)^((properties.gamma_air+1)/(2*(properties.gamma_air-1)));
    diffuser.A1 = diffuser.A2*diffuser.Rd*term1*term2;
    %calculate diffuser inlet diameter
    diffuser.D1 = 2*sqrt(diffuser.A1/pi);
    %define diffuser cz2
    diffuser.cz2 = compressor.cz; %since const cz across IGV
end