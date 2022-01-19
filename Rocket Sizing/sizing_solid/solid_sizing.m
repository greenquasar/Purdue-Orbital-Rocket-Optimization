%%INPUTS
    %%STATIC:
        % gravity (accel units)
        % inert mass of solid stage (mass of the hull) (kg)
        % atmospheric pressure (Pa)
        % fuel temp (298 K)
        % oxidizer temp (298 K)
        % time step (s)
        % Fuel Type
        % Oxidizer Type
        % Accelerant
    %%VARYING EACH LOOP:
        % Length (meters)
        % Width (meters)
        % Inner diameter (meters)
        % Maximum Pressure (Pa)
%%OUTPUTS
    % T = Time Array [seconds]
    % W = Web Array (Inner Radius array) [meters] <<<<< look at this >>>>> !! limiting factor !!
    % P_c = Chamber Pressure over time [Pa]
    % Thrust = Thrust Array [Newtons]
    % R_b = Burn Rate [Meters / Seconds]
    % burn_time = Burn Time (End value of Time array) [seconds]
    % M = Mass over time [Kg] <<<< look at initial value of this one for our purposes >>>> !! optimization !!
    % Mdot = Mass flow over time [Kg/s]
    % A_t = Throat Area [meters^2]
    % deltav = Delta V [Km/s] <<<<<<< look at this value for our purposes >>>>>> !! optimization!!
    % specificImpulse = average specific impulse (seconds)
    % propMass = mass of propellant [kg]
    % C_t = Coefficient of Thrust [no units]
    % C_star = Characteristic Temperature [no units?] 
    %


function out = solid_sizing(minTWR, maxTWR, maxPres, deltaV, diameterMin, diameterMax, lengthMin, lengthMax, webMin, webMax, OF)
    addpath('../../SRM Sizing');
    dx = 0.
    1; %diameter step, m
    dy = 0.1; %length step, m
    dw = 0.1; %web step, m
    shape = 'circular';
    fuel = 'HTPB';
    oxidizer = 'NH4CLO4(I)';
    for diameter = diameterMin:dx:diameterMax
        disp(diameter);
        for length = lengthMin:dy:lengthMax
            disp(length);
            for web = webMin:dw:webMax
                disp(web);
                if web < diameter
                    [T, W, P_c, Thrust, R_b, burn_time, M, Mdot, A_t, deltaV, specificImpulse, propMass, C_t, C_star] = Simulate_Reverse(shape, length, diameter, web, maxPres, OF, fuel, oxidizer);
                end
            end
        end
    end
    out = 1;  %todo
    
end