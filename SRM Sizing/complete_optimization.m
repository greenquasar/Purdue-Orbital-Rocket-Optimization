%[length2, width2, inner_width2, newfinal_simulation2, mass2, length1, width1, inner_width1, newfinal_simulation1, mass1] =  complete_optimization(1, 4000, 1.5, 3447378.64659, 'circular', 0.077, 5, 10407, 2.22, ["HTPB"], [298], [920], [1], ["NH4CLO4(I)", "AL"], [298, 200], [1950,2710], [0.88, 0.12], 0.12, 2.75)

%30-70 deltaV or 50-50 deltaV
%solid rocket motor sizing 
function [length2, width2, inner_width2, newfinal_simulation2, mass2, ...
    length1, width1, inner_width1, newfinal_simulation1, mass1] = ...
    complete_optimization(dt, delta_V, TWR, ...
    maxPres, shape, f_inert, payloadMass, atmoPressure, ...
    OF, fuels, f_temps, f_densities, f_fracs, oxidizers, ...
    o_temps, o_densities, o_fracs, diaU, lenU)
%% Inputs
    %dt: time step (s)
    %deltaV: desired deltaV (m/s)
    %TWR: thrust to weight ratio (1)
    %shape: circular/square shape: "circular"/"square" 
    %maxPres: maximum chamber pressure stage can withstand (Pa)
    %f_inert: inert mass fraction
    %shape: circular/square shape: "circular"/"square"
    %atmoPressure: atmospheric pressure (Pa)
    %OF: oxidizer to fuel ratio
    %fuel: CEA fuel name, string
    %f_t: fuel inlet temp (K)
    %f_dens: fuel density in solid state (kg/m^3)
    %oxidizer: CEA oxidizer name, string
    %o_t: oxidizer inlet temp (K)
    %o_dens: fuel density in solid state (kg/m^3)
    %o_fracs: mass fractions for each oxidizer (same length as oxidizer, sum to 1)
%% Outputs
    %length: length of stage (m)
    %width: width of stage (m)
    %inner_width: inner width of stage (m)
    %final_simulation: Simulate_Reverse final output as a list
%% Program
clear
clc
%Optimize Stage 2
[length2, width2, inner_width2, final_simulation2, mass2, deltaVWorking, massWorking, diameterWorking, lengthWorking, inradWorking] = optimize(1, 4000, 1.5, 3447378.64659, 'circular', 0.077, 1.33, 10407, 2.22, ["HTPB"], [298], [920], [1], ["NH4CLO4(I)", "AL"], [298, 200], [1950,2710], [0.88, 0.12], 0.12, 2.75)
fprintf('Stage 2 length is %0.4f\n', length2);
fprintf('Stage 2 width is %0.4f\n', width2);
fprintf('Stage 2 inner width is %0.4f\n', inner_width2);
fprintf('Stage 2 mass is %0.4f\n', mass2);
stage1payloadmass = mass2 + 1.33;
%Optimize Stage 1
[length1, width1, inner_width1, final_simulation1, mass1, newdeltaVWorking, newmassWorking, newdiameterWorking, newlengthWorking, newinradWorking] = optimize(1, 4000, 1.5, 3447378.64659, 'circular', 0.077, stage1payloadmass, 10407, 2.22, ["HTPB"], [298], [920], [1], ["NH4CLO4(I)", "AL"], [298, 200], [1950,2710], [0.88, 0.12], 0.15, 10)
fprintf('Stage 1 length is %0.4f\n', length1);
fprintf('Stage 1 width is %0.4f\n', width1);
fprintf('Stage 1 inner width is %0.4f\n', inner_width1);
fprintf('Stage 1 mass is %0.4f\n', mass1);
end