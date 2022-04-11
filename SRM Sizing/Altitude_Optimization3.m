%optimization for single stage sounding rocket
function [stage_length, stage_width, inner_width, ...
    T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, avgSpecificImpulse, propMass, C_t, C_star, Altitude, Drag, Velocity, final_altitude] = ...
    Altitude_Optimization3(target_altitude, starting_altitude, TWRmin, dynamicPressMax, tolerance, ...
    dt, maxPres, shape, f_inert, payloadMass, atmoPressure, OF, ...
    fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs)
%% Inputs
    %target_altitude: target altitude (m)
    %starting_altitude: Starting altitude (m)
    %TWRmin: minimum Thrust to Weight Ratio
    %TWRmax: maximum Thrust to Weight Ratio
    %tolerance: acceptable error in TWR and altitude (fraction)
    %dt: time step (s)
    %maxPres: maximum chamber pressure stage can withstand (Pa)
    %shape: circular/square shape: "circular"/"square"
    %f_inert: inert mass fraction
    %payloadMass: mass of payload (kg)
    %atmoPressure: atmospheric pressure (Pa)
    %OF: oxidizer to fuel ratio
    %fuels: CEA fuel names, string (can be a list must use double quotes)
    %f_temps: fuel inlet temps (K) (same length as fuel)
    %f_densities: fuel densities in solid state (kg/m^3) (same length as fuel)
    %f_fracs: mass fractions for each fuel (same length as fuel, sum to 1)
    %oxidizers: CEA oxidizer names, string (can be a list must use double quotes)
    %o_temps: oxidizer inlet temps (K) (same length as oxidizer)
    %o_densities: fuel densities in solid state (kg/m^3) (same length as oxidizer)
    %o_fracs: mass fractions for each oxidizer (same length as oxidizer, sum to 1)
    
%% Outputs
    %stage_length, length of stage (m)
    %stage_width, diameter of stage (m)
    %inner_width, inner diameter of fuel bore (m)
    %T, array of time steps (s)
    %W, array of web distance (m)
    %P_c, array of chamber pressure (Pa)
    %Thrust, array of Thrust (N)
    %TWR, Thrust to Weight Ratio (1)
    %R_b, array of burn rate ()
    %burn_time, time to expend stage (s)
    %M, array of mass (kg)
    %Mdot, array of mass flow
    %A_t, throat area required to meet max pressure constraint at end of burn(m^2)
    %deltaV, delta V of stage (m/s)
    %avgSpecificImpulse, average of Isp over burn (s)
    %propMass, (kg)
    %C_t, array of coefficient of thrust (1)
    %C_star, array of chamber pressure (Pa)
    %Altitude, array of altitude (m)
    %Drag, array drag forces (N)
    %Velocity, array of velocities (m/s)

    
%% Constants



%% Program

    %passthrough_args = [maxPres, shape, f_inert, payloadMass, atmoPressure, OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs];
    
    %guesses
    stage_length = 0.3007;
    stage_width  = 0.1733;
    inner_width  = 0.0980;
    
    %do while loop
    while (1)
        [T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, avgSpecificImpulse, propMass, C_t, C_star] = ...
    Simulate_Reverse(dt, stage_length, stage_width, inner_width, maxPres, shape, f_inert, payloadMass, atmoPressure, ...
    OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs);
    [Altitude, Drag, Velocity, DynamicPressure] = altitude_analysis(Thrust, M, dt, stage_width, starting_altitude);
        final_altitude = max(Altitude);
        dynamic_pressure = max(DynamicPressure);
        RocketDrawer(stage_length, stage_width, inner_width);
        altitude_pct_error = 100 * abs(final_altitude-target_altitude)/target_altitude;
        maxQ_pct_error = 100 * abs(dynamic_pressure-dynamicPressMax)/dynamicPressMax;
        fprintf("Stage Length = %.4f | Stage Width = %.4f | Inner Width = %.4f\n", stage_length, stage_width, inner_width);
        fprintf("Target Altitude = %.2f | Current Altitude = %.2f | Error = %0.2f%%\n", target_altitude, final_altitude, altitude_pct_error);
        fprintf("Target max Q = %.2f | Current max Q = %.2f | Error = %0.2f%%\n", dynamicPressMax, dynamic_pressure, maxQ_pct_error);
        if (~tolerated(final_altitude, target_altitude, tolerance) || ~tolerated(dynamicPressure, dynamicPressMax, tolerance)) %altitude          
            stage_length_factor = 1;
            stage_width_factor = 1;
            inner_width_factor = 1;
            
            altFactor = (target_altitude/final_altitude)^(1/3);
            qFactor = (dynamicPressMax/dynamicPressure)^(1/3);
            
            
            %%todo the bois
            
            fprintf("Scaling Stage Length by %0.4f\n", stage_length_factor);
            fprintf("Scaling Stage Width by %0.4f\n", stage_width_factor);
            fprintf("Scaling Inner Width by %0.4f\n", inner_width_factor);
            
            stage_length = stage_length * stage_length_factor;
            stage_width = stage_width * stage_width_factor;
            inner_width = inner_width * inner_width_factor;
            
            
            
        else
            break;
        end
    end
    
end

function result = tolerated(value, target, tolerance)
    result = abs(target-value)/target < tolerance;
end
