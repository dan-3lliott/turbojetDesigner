function compressor = comp(Pr1_comp, diffuser, givens, properties, compressor)
    %prepare dummy stage
    compressorStages(1) = dummy_stage(properties, compressor);

    %prepare compressor structure for stages
    compressor.Pr = 1;
    compressor.W_total = 0;

    %assume linear pressure ratio distribution
    fun = @(d) prod(Pr1_comp - d*(0:compressor.Nstages-1)) - (compressor.Prreq/(compressor.Pr_igv*compressor.Pr_ogv));
    d = fzero(fun, [0, Pr1_comp/(compressor.Nstages-1)]);
    Prs_comp = Pr1_comp - d*(0:compressor.Nstages-1);

    amax = (compressor.a3StageMax) - (compressor.a3a1offset * compressor.Nstages);

    %optimize first stage
    [compressorStages(2), compressor.rt, compressor.N, a1_before_first_stage, compressor.cz, compressor.M2_diff] = first_comp_stage(Pr1_comp, diffuser, givens, properties, compressor, compressorStages(1), amax);
    
    %go back and update first stage to be the inlet guide vane
    compressorStages(1) = igv_comp(compressor.Pr_igv, givens, properties, compressor, a1_before_first_stage, compressor.cz, properties.To2, properties.Po2, compressor.ai);
    compressor.Pr = compressor.Pr * compressorStages(1).Pr * compressorStages(2).Pr;
    compressor.W_total = compressor.W_total + compressorStages(2).Wdot;

    for i = 3:(compressor.Nstages+1)
        %optimize individual stage
        compressorStages(i) = comp_stage(Prs_comp(i-1), givens, properties, compressor, compressorStages(i-1));
        %update total compressor parameters
        compressor.Pr = compressor.Pr * compressorStages(i).Pr;
        compressor.W_total = compressor.W_total + compressorStages(i).Wdot;
    end
    %optimize outlet guide vane
    compressorStages(end+1) = ogv_comp(compressor.Pr_ogv, givens, properties, compressor, compressor.af, compressor.cz, 1, compressorStages(end).To3, compressorStages(end).Po3, compressorStages(end).a3);
    compressor.Pr = compressor.Pr * compressorStages(end).Pr;

    %final properties for output
    compressor.Po3 = compressorStages(end).Po3;
    compressor.To3 = compressorStages(end).To3;

    %total efficiency for output
    compressor.eta = ((properties.gamma_air-1)/properties.gamma_air)*log(compressor.Pr)/log(compressor.To3/properties.To2);

    %store compressor stages in compressor structure
    compressor.stages = compressorStages;
end