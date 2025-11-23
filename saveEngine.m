function engine = saveEngine(givens, properties, diffuser, compressor, combustor, turbine, nozzle)
    %{
    DEBUGGING START
    figure
    hold on
    plot(1:length(compressor.stages), [compressor.stages.rhr]);
    plot(1:length(compressor.stages), ones(length(compressor.stages),1)*compressor.rt);
    hold off
    xlabel('Stage');
    ylabel('Radial Distance');
    legend('Hub', 'Tip', 'location', 'southeast');
    ylim([0 compressor.rt+0.1]);

    figure
    hold on
    yyaxis left;
    plot(0:length(compressor.stages), [properties.Po1 compressor.stages.Po3]);
    yyaxis right;
    plot(0:length(compressor.stages), [properties.To1 compressor.stages.To3]);
    hold off
    legend('Po', 'To');

    figure
    hold on
    plot(1:length(compressor.stages), [compressor.stages.a3]);
    plot(1:length(compressor.stages), [compressor.stages.a2]);
    hold off
    legend('a3', 'a2');


    disp('Compressor cz:');
    disp(compressor.cz);
    disp('Compressor RPM:');
    disp(compressor.N);
    disp('Compressor Power:');
    disp(compressor.W_total);
    disp('Compressor Exit M:');
    disp(compressor.stages(end).M1rel);
    %}

    %off-design calculations
    od = offDesign(givens, properties, diffuser, compressor, combustor, turbine, nozzle);
    writetable(struct2table(od), 'outputs/offdesign.xlsx');

    disp('Diffuser Structure:');
    disp(diffuser);

    disp('Compressor Structure:');
    disp(compressor);

    disp('Combustor Structure:');
    disp(combustor);

    disp('Turbine Structure:');
    disp(turbine);

    disp('Nozzle Structure:');
    disp(nozzle);

    disp('Off-Design Structure:');
    disp(od);
    %DEBUGGING STOP

    %graph the entire engine
    %diffuser
    casing_radii = [diffuser.D1/2];
    hub_radii = [0];
    x_pos = [0 (compressor.rt - diffuser.D1/2)/tand(10)];
    %compressor casing - super simple
    casing_radii = [casing_radii ones(1,(2*length(compressor.stages))).*compressor.rt];
    %igv
    hub_radii = [hub_radii compressor.stages(1).rhs];
    nextX = x_pos(end) + ((compressor.rt - compressor.stages(1).rhs)/compressor.stages(1).bladeARs)*abs(cosd(compressor.stages(1).a3 - compressor.stages(1).a2));
    x_pos = [x_pos nextX];
    hub_radii = [hub_radii compressor.stages(1).rhs];
    nextX = x_pos(end) + ((compressor.rt - compressor.stages(1).rhs)/compressor.stages(1).bladeARs)*abs(cosd(compressor.stages(1).a3 - compressor.stages(1).a2));
    x_pos = [x_pos nextX];
    %compressor
    for i = 2:(length(compressor.stages)-1)
        nextX = x_pos(end) + ((compressor.rt - compressor.stages(i).rhr)/compressor.stages(i).bladeARr)*abs(cosd(compressor.stages(i).a2 - compressor.stages(i).a1));
        x_pos = [x_pos nextX];
        hub_radii = [hub_radii compressor.stages(i).rhr compressor.stages(i).rhs];
        nextX = x_pos(end) + ((compressor.rt - compressor.stages(i).rhs)/compressor.stages(i).bladeARs)*abs(cosd(compressor.stages(i).a3 - compressor.stages(i).a2));
        x_pos = [x_pos nextX];
    end
    %ogv
    hub_radii = [hub_radii compressor.stages(end).rhs];
    nextX = x_pos(end) + ((compressor.rt - compressor.stages(end).rhs)/compressor.stages(end).bladeARs)*abs(cosd(compressor.stages(end).a3 - compressor.stages(end).a2));
    x_pos = [x_pos nextX];
    hub_radii = [hub_radii compressor.stages(end).rhs];
    nextX = x_pos(end) + ((compressor.rt - compressor.stages(end).rhs)/compressor.stages(end).bladeARs)*abs(cosd(compressor.stages(end).a3 - compressor.stages(end).a2));
    x_pos = [x_pos nextX];
    %combustor
    x_pos(end) = x_pos(end) + combustor.L_diff;
    casing_radii = [casing_radii combustor.ro3e combustor.ro3e];
    hub_radii = [hub_radii combustor.ri3e combustor.ri3e];
    x_pos = [x_pos x_pos(end) + (combustor.L - combustor.L_diff) x_pos(end) + (combustor.L - combustor.L_diff + 0.015)];
    %turbine
    for i = 1:length(turbine.allStages)
        casing_radii = [casing_radii turbine.allStages(i).rtn turbine.allStages(i).rtr];
        hub_radii = [hub_radii turbine.allStages(i).rh turbine.allStages(i).rh];
        nextX = x_pos(end) + ((turbine.allStages(i).rtn - turbine.allStages(i).rh)/turbine.allStages(i).bladeARn)*abs(cosd(turbine.allStages(i).a2 - turbine.allStages(i).a1));
        x_pos = [x_pos nextX];
        nextX = x_pos(end) + ((turbine.allStages(i).rtr - turbine.allStages(i).rh)/turbine.allStages(i).bladeARr)*abs(cosd(turbine.allStages(i).a3 - turbine.allStages(i).a2));
        x_pos = [x_pos nextX];
    end
    %turb ogv
    casing_radii = [casing_radii turbine.ogv.rtn];
    hub_radii = [hub_radii turbine.ogv.rhn];
    nextX = x_pos(end) + ((turbine.ogv.rtn - turbine.ogv.rhn)/turbine.ogv.bladeARn)*abs(cosd(turbine.ogv.a_out - turbine.ogv.a_in));
    x_pos = [x_pos nextX];
    casing_radii = [casing_radii turbine.ogv.rtn];
    hub_radii = [hub_radii turbine.ogv.rhn];
    nextX = x_pos(end) + ((turbine.ogv.rtn - turbine.ogv.rhn)/turbine.ogv.bladeARn)*abs(cosd(turbine.ogv.a_out - turbine.ogv.a_in));
    rnozz = sqrt(nozzle.Ae/pi);
    x_pos = [x_pos (nextX + (turbine.ogv.rtn - rnozz)/tand(20))];
    %nozzle
    casing_radii = [casing_radii sqrt(nozzle.Ae/pi)];
    hub_radii = [hub_radii 0];

    x_stations = (1:length(hub_radii)) - 1;

    figure
    hold on
    plot(x_pos, casing_radii, 'Color', ([138, 43, 226]./255));
    plot(x_pos, -casing_radii, 'Color', ([138, 43, 226]./255));
    x_pos(1) = x_pos(2);
    x_pos(end) = x_pos(end-1);
    plot(x_pos, hub_radii, 'Color', ([34, 154, 109]./255));
    plot(x_pos, -hub_radii, 'Color', ([34, 154, 109]./255));
    axis equal
    hold off
    xlim([0 1.32]);
    xlabel('x (m)');
    ylabel('y (m)');

    %output data files for each structure
    writetable(struct2table(diffuser), 'outputs/diffuser.xlsx');
    writetable(struct2table(compressor), 'outputs/compressor.xlsx');
    for i = 1:length(compressor.stages)
        filename = strcat('outputs/compressor_stage_', num2str(i-1), '.xlsx'); %starts IGV at 0
        writetable(struct2table(compressor.stages(i)), filename);
    end
    writetable(struct2table(combustor), 'outputs/combustor.xlsx');
    writetable(struct2table(turbine), 'outputs/turbine.xlsx');
    for i = 1:length(turbine.allStages)
        filename = strcat('outputs/turbine_stage_', num2str(i), '.xlsx');
        writetable(struct2table(turbine.allStages(i)), filename);
    end
    writetable(struct2table(turbine.ogv), 'outputs/turbine_ogv.xlsx');
    writetable(struct2table(nozzle), 'outputs/nozzle.xlsx');

    %final engine structure
    engine.givens = givens;
    engine.properties = properties;
    engine.diffuser = diffuser;
    engine.compressor = compressor;
    engine.combustor = combustor;
    engine.turbine = turbine;
    engine.nozzle = nozzle;
    engine.od = od;
    general.ST = (nozzle.thrust/givens.mdot_air);
    general.SFC = (combustor.f/general.ST);
    engine.general = general;
end