%% solid rocket motor sizing code
function [length, width, inner_width, final_simulation] = ...
    optimize(dt, delta_V, TWR, ...
    maxPres, shape, f_inert, payloadMass, atmoPressure, ...
    OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs, diaU, lenU)
    %% Inputs
    %dt: time step (s)
    % test
    %error_tolerance: acceptable error in deltaV and TWR (%)
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
    %passthrough_args = [maxPres, shape, f_inert, payloadMass, atmoPressure, OF, fuel, f_temp, f_dens, oxidizer, o_temp, o_dens]
    
    max_iterations = 10;
    %deltaV should be controlled by length??
    %TWR should be conrolled by width and inner radius??

    %start with a guess for length and width
    %loop if iter less than max_iter
    diaL = 0.1 * diaU;
    lenL = 0.1 * lenU;
    index = 0;
    NewEntry = [0, 0, 0, 0, 0];
    for dia = linspace(diaL, diaU, max_iterations)
        inradU = (dia / 2) * 0.8;
        inradL = 0.03 * (dia / 2);
        for len = linspace(lenU,lenL, max_iterations)
            for inrad = linspace(inradL, inradU, max_iterations)
                index = index + 1;
                [T, W, P_c, Thrust, R_b, burn_time, M, Mdot, A_t, deltaV, specificImpulse, propMass, C_t, C_star] = Simulate_Reverse(dt, len, dia, inrad, maxPres, shape, f_inert, payloadMass, atmoPressure, OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs);
                NewEntry = [dia, len, inrad, M, deltaV];
                LoopResults(index, :) = NewEntry;
            end
        end
    end
    
    indexx = find(LoopResults(:,4) >= delta_V)
    massWorking = min(LoopResults(:,4));
    deltaVWorking = LoopResults(min(LoopResults(:,4)), 5);
    diameterWorking = LoopResults(min(LoopResults(:,4)), 1);
    lengthWorking = LoopResults(min(LoopResults(:,4)), 2);
    inradWorking = LoopResults(min(LoopResults(:,4)), 3);
    
    mass = min(massWorking);
    width = diameterWorking(mass == massWorking);
    length = lengthWorking(mass == massWorking);
    inner_width = inradWorking(mass == massWorking);

end