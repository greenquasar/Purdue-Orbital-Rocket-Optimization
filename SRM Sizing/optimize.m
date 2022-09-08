%[length, width, inner_width, final_simulation, mass, deltaVWorking, massWorking, diameterWorking, lengthWorking, inradWorking] = optimize(1, 4000, 1.5, 3447378.64659, 'circular', 0.077, 5, 10407, 2.22, ["HTPB"], [298], [920], [1], ["NH4CLO4(I)", "AL"], [298, 200], [1950,2710], [0.88, 0.12], 0.12, 2.75, 187.147)

%% solid rocket motor sizing code
function [length, width, inner_width, final_simulation, mass, ...
    deltaVWorking, massWorking, diameterWorking, lengthWorking, ...
    inradWorking] = optimize(dt, delta_V, TWR, ...
    maxPres, shape, f_inert, payloadMass, atmoPressure, ...
    OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, ...
    o_densities, o_fracs, diaU, lenU, startingalt)
    %% Inputs
    %dt: time step (s)
    % test
    %error_tolerance: acceptable error in deltaV and TWR (%)
    %deltaV: desired deltaV (m/s)
    %TWR: thrust to weight ratio (1)
    %shape: circular/square shape: "circular"/"square" 
    %maxPres: maximum chamber pressure stage can withstand (Pa)
    %f_inert: inert mass fraction
    %payloadMass (kg)
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
    
    max_iterations = 2;
    %deltaV should be controlled by length??
    %TWR should be controlled by width and inner radius??

    %start with a guess for length and width
    %loop if iter less than max_iter
    diaL = 0.1 * diaU; 
    lenL = 0.1 * lenU; 
    
    %upper values needed to be added
    %upper diameter = 8 inches
    %upper length ?
    
    index = 0;
    NewEntry = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    LoopResults = zeros(9,9);
    for dia = linspace(diaL, diaU, max_iterations)
        inradU = (dia / 2) * 0.8;
        inradL = 0.03 * (dia / 2);
        for len = linspace(lenU,lenL, max_iterations)
            for inrad = linspace(inradL, inradU, max_iterations)
                index = index + 1;
                [T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, specificImpulse, propMass, C_t, C_star] = Simulate_Reverse(dt, len, dia, inrad, maxPres, shape, f_inert, payloadMass, atmoPressure, OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs);
                %added altitude analysis here to account for altitude and
                %dynamic pressure
                [alt, drag, velocity, DP] = altitude_analysis(Thrust, M, dt, inrad * 2, startingalt);
                NewEntry = [dia, len, inrad, M(1), deltaV, alt, DP];
                LoopResults(index, :) = NewEntry;
                fprintf('This is iteration #%0.0f.\n', index);
            end
        end
    end
    
    %select best candidate based off the altitude, dynamic pressure, and
    %mass
    
    %use altitude and dynamic pressure instead of deltaV
    %indexDW = LoopResults(:,5) > delta_V; 
    indexDW = find(LoopResults(:,6) > 3048 + startingalt & LoopResults(:,7) < (maxPres/1000));

    %deltaVWorking = LoopResults(indexDW, 5);
    massWorking = LoopResults(indexDW,4);
    diameterWorking = LoopResults(indexDW, 1);
    lengthWorking = LoopResults(indexDW, 2);
    inradWorking = LoopResults(indexDW, 3);
    
    mass = min(massWorking)
    width = diameterWorking(mass == massWorking);
    length = lengthWorking(mass == massWorking);
    inner_width = inradWorking(mass == massWorking);

    final_simulation = 1;
    
    %output thrust to weight ratio
    fprintf('The thrust to weight ratio will be: %.2f/n', TWR);
    
    WritetoExcel(deltaVWorking, massWorking, diameterWorking, lengthWorking, inradWorking);
end