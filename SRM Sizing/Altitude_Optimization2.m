%optimization for single stage sounding rocket
function [stage_length, stage_width, inner_width, ...
    T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, avgSpecificImpulse, propMass, C_t, C_star, Altitude, Drag, Velocity, final_altitude] = ...
    Altitude_Optimization2(target_altitude, starting_altitude, TWRmin, TWRmax, tolerance, ...
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
    n = 5; %grid density

%% Program

    %passthrough_args = [maxPres, shape, f_inert, payloadMass, atmoPressure, OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs];
    
    %guesses
    stage_length = 0.3050;
    stage_width  = 0.1700;
    inner_width  = 0.0950;
    
    range_factor = 1.5;
    stage_length_low  = stage_length/range_factor;
    stage_length_high = stage_length*range_factor;
    stage_width_low   = stage_width/range_factor;
    stage_width_high  = stage_width*range_factor;
    inner_width_low   = inner_width/range_factor;
    inner_width_high  = inner_width*range_factor;
    
    %stage_length_low  = 0.28;
    %stage_length_high = 0.33;
    %stage_width_low   = 0.14;
    %stage_width_high  = 0.18;
    %inner_width_low   = 0.08;
    %inner_width_high  = 0.10;
    
    %do while loop
    my_waitbar = waitbar(0, "Initial Simulation", 'Name', "Altitude_Optimization2");
    while (1)  
        [T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, avgSpecificImpulse, propMass, C_t, C_star] = ...
        Simulate_Reverse(dt, stage_length, stage_width, inner_width, maxPres, shape, f_inert, payloadMass, atmoPressure, ...
        OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs);
        [Altitude, Drag, Velocity] = altitude_analysis(Thrust, M, dt, stage_width, starting_altitude);
        final_altitude = max(Altitude);
        
        altitude_pct_error = 100 * abs(final_altitude-target_altitude)/target_altitude;
        TWRmax_pct_error = 100 * abs(max(TWR)-TWRmax)/TWRmax;
        TWRmin_pct_error = 100 * abs(TWR(1)-TWRmin)/TWRmin;
        
        fprintf("Stage Length = %.4f | Stage Width = %.4f | Inner Width = %.4f\n", stage_length, stage_width, inner_width);
        fprintf("Target Altitude = %.2f | Current Altitude = %.2f | Error = %0.2f%%\n", target_altitude, final_altitude, altitude_pct_error);
        fprintf("Target max TWR = %.2f | Current max TWR = %.2f | Error = %0.2f%%\n", TWRmax, max(TWR), TWRmax_pct_error);
        fprintf("Target min TWR = %.2f | Current min TWR = %.2f | Error = %0.2f%%\n", TWRmin, TWR(1), TWRmin_pct_error);
        
        waitbar(0, my_waitbar, "Current Iteration Progress:");
        
        if (tolerated(final_altitude, target_altitude, tolerance) && tolerated(max(TWR), TWRmax, tolerance) && tolerated(TWR(1), TWRmin, tolerance))      
            break;
            
        else
            %insert codeblock here
            Errors = ones(n,n,n);
            %fill dimensions
            %\/\/ THE TESSERACT OF POSSIBILITIES \/\/
            Dimensions = zeros(n,n,n,3); %val(:,:,:,1) = stage_length | val(:,:,:,2) = stage_width | val(:,:,:,3) = inner_width
            %/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            %stage length
            data = linspace(stage_length_low, stage_length_high, n);
            for x = 1:n
                for y = 1:n
                    for z = 1:n
                        Dimensions(x,y,z,1) = data(x);
                    end
                end
            end
            %stage width
            data = linspace(stage_width_low, stage_width_high, n);
            for x = 1:n
                for y = 1:n
                    for z = 1:n
                        Dimensions(x,y,z,2) = data(y);
                    end
                end
            end
            %inner width
            data = linspace(inner_width_low, inner_width_high, n);
            for x = 1:n
                for y = 1:n
                    for z = 1:n
                        Dimensions(x,y,z,3) = data(z);
                    end
                end
            end
            
            fprintf("Searching in range:\nStage Length: [%.4f, %.4f]\nStage Width:  [%.4f, %.4f]\nInner Width:  [%.4f, %.4f]\n",...
                stage_length_low, stage_length_high, stage_width_low, stage_width_high, inner_width_low, inner_width_high);
            
            %run iterations
            for x = 1:n
                for y = 1:n
                    for z = 1:n
                        %fprintf("(%d,%d,%d)",x,y,z);
                        stage_length = Dimensions(x,y,z,1);
                        stage_width = Dimensions(x,y,z,2);
                        inner_width = Dimensions(x,y,z,3);
                        
                        if (inner_width >= stage_width)
                            error = nan;
                        else 
                            [T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, avgSpecificImpulse, propMass, C_t, C_star] = ...
                            Simulate_Reverse(dt, stage_length, stage_width, inner_width, maxPres, shape, f_inert, payloadMass, atmoPressure, ...
                            OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs);
                            [Altitude, Drag, Velocity] = altitude_analysis(Thrust, M, dt, stage_width, starting_altitude);
                            final_altitude = max(Altitude);
                            
                            waitbar(((x-1)*(n^2)+(y-1)*n+z)/(n^3), my_waitbar, sprintf("Current Iteration Progress: (%d/%d)",((x-1)*(n^2)+(y-1)*n+z), n^3));
                            
                            altitude_pct_error = 100 * abs(final_altitude-target_altitude)/target_altitude;
                            TWRmax_pct_error = 100 * abs(max(TWR)-TWRmax)/TWRmax;
                            TWRmin_pct_error = 100 * abs(TWR(1)-TWRmin)/TWRmin;

                            error = altitude_pct_error^2 + TWRmax_pct_error^2 + TWRmin_pct_error^2;
                        end
                        Errors(x,y,z) = error;
                        
                    end
                end
            end
            fprintf("\n");
            %find lowest error and set new limits
            [~,loc] = min(Errors(:));
            [ii,jj,kk] = ind2sub(size(Errors),loc);
            stage_length = Dimensions(ii,jj,kk,1);
            stage_width  = Dimensions(ii,jj,kk,2);
            inner_width  = Dimensions(ii,jj,kk,3);
            if (ii == 1 || ii == n)
                disp("Stage Length range too restrictive");
            elseif (jj == 1 || jj == n)
                disp("Stage Width range too restrictive");
            elseif (kk == 1 || kk == n)
                disp("Inner Width range too restrictive");
            end
            stage_length_low  = Dimensions(ii-1,jj,kk,1);
            stage_length_high = Dimensions(ii+1,jj,kk,1);
            stage_width_low   = Dimensions(ii,jj-1,kk,2);
            stage_width_high  = Dimensions(ii,jj+1,kk,2);
            inner_width_low   = Dimensions(ii,jj,kk-1,3);
            inner_width_high  = Dimensions(ii,jj,kk+1,3);
        end
    end
    close(my_waitbar);
    
end

function result = tolerated(value, target, tolerance)
    result = abs(target-value)/target < tolerance;
end
