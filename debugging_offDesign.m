clear all

engine = load('engine_for_od_debugging.mat');
givens = load('givens_for_od_debugging.mat');
properties = load('properties_for_od_debugging.mat');

od = offDesign(engine.finalEngine.givens, engine.finalEngine.properties, engine.finalEngine.diffuser, engine.finalEngine.compressor, engine.finalEngine.combustor, engine.finalEngine.turbine, engine.finalEngine.nozzle)