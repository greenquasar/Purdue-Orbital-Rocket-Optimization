%[length, width, inner_width, final_simulation] = optimize(1, 4000, 1.5, 3447378.64659, 'circular', 0.077, 5, 8161.4, 2.22, ["HTPB"], [298], [920], [1], ["NH4CLO4(I)", "AL"], [298, 200], [1950,2710], [0.88, 0.12], 0.2032, 10)

%solid rocket motor sizing 
function [length, width, inner_width, final_simulation] = ...
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
%Optimize Stage 2
[length, width, inner_width, final_simulation] = optimize(1, 4000, 1.5, 3447378.64659, 'circular', 0.077, 5, 8161.4, 2.22, ["HTPB"], [298], [920], [1], ["NH4CLO4(I)", "AL"], [298, 200], [1950,2710], [0.88, 0.12], 0.2032, 10)
