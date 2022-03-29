%optimization for single stage sounding rocket
function [] = Altitude_Optimization(target_altitude, starting_altitude, TWRmin, TWRmax, ...
    dt, maxPres, shape, f_inert, atmoPressure, OF, ...
    fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs)
%%inputs
    %target_altitude: target altitude (m)
    %starting_altitude: Starting altitude (m)
    %TWRmin: minimum Thrust to Weight Ratio
    %TWRmax: maximum Thrust to Weight Ratio
    %dt: time step (s)
    %maxPres: maximum chamber pressure stage can withstand (Pa)
    %shape: circular/square shape: "circular"/"square"
    %f_inert: inert mass fraction
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
    
%%outputs
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

end
