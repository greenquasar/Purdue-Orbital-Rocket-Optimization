%optimization for single stage sounding rocket
function [T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, avgSpecificImpulse, propMass, C_t, C_star, Altitude, Drag, Velocity, final_altitude] = ...
    Altitude_Optimization(target_altitude, starting_altitude, TWRmin, TWRmax, tolerance, ...
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
    stage_length = 1;
    stage_width = 0.1;
    innerWidth = 0.05;
    
    %do while loop
    while (1) 
        [T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, avgSpecificImpulse, propMass, C_t, C_star] = ...
    Simulate_Reverse(dt, stage_length, stage_width, innerWidth, maxPres, shape, f_inert, payloadMass, atmoPressure, ...
    OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs);
    [Altitude, Drag, Velocity] = altitude_analysis(Thrust, M, dt, stage_width, starting_altitude);
        final_altitude = max(Altitude);
        if (tolerated(final_altitude, target_altitude, tolerance)) %altitude
            %maybe need to calculate how much to change the width or the
            %length to keep volume const
            %stage width
            factor = target_altitude/final_altitude;
            stage_width = stage_width * factor;
        elseif (tolerated(max(TWR), TWRmax, tolerance)) %max TWR
            %stage length
            factor = TWRmax/max(TWR);
            stage_length = stage_length * factor;
        elseif (tolerated(TWR(1), TWRmin, tolerance)) %min TWR
            %inner width
            factor = TWRmin/TWR(1);
            inner_width = inner_width * factor;
        else
            break;
        end
    end
    
end

function ans = tolerated(value, target, tolerance)
    ans = abs(target-value)/target < tolerance;
end
