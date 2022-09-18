%% solid rocket motor simulation code
function [T, W, P_c, Thrust, TWR, R_b, burn_time, M, Mdot, A_t, deltaV, avgSpecificImpulse, propMass, C_t, C_star] = ...
    Simulate_Reverse(dt, stage_length, stage_width, innerWidth, maxPres, shape, f_inert, payloadMass, atmoPressure, ...
    OF, fuels, f_temps, f_densities, f_fracs, oxidizers, o_temps, o_densities, o_fracs) 
    %% Inputs
    
    %dt: time step (s)
    %stage_length: chamber length (m)
    %stage_width: outer chamber width (m)
    %innerWidth: inner chamber (m)
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
    
    %% Program
    
    % Constants (SI)
    g = 9.80665;    %[m/s^2]
    


    %GUESSTIMATED DUMMY VALUES FROM CEARUN, get better ones from CEA

    n = 0.098; %n is dimensionless
    C = 6.2; %C has units of [m^(2n+1)N^(-n)] at pressure of 1000 psi (6.895 MPA)

    C_t = [];
    C_star = [];
    Isp = [];
   
    [C_t(1), C_star(1), Isp(1)] = SRM_CEA(maxPres,OF,fuels,f_temps,f_fracs,oxidizers,o_temps,o_fracs,atmoPressure);
    disp(C_star(1));
    %error checking
    if (stage_width <= innerWidth)
       error("Width is less than or equal to inner width. Impossible!");
    end
    %todo if fuel and ox lists aren't the same length
    
    %convenience terms by manipulating inputs
    r_max = stage_width / 2;  %[m]
    r_min = innerWidth / 2; %[m]
    
    f_dens = 0;
    for i = 1:length(fuels)
        f_dens = f_dens + f_densities(i) * f_fracs(i);
    end
    o_dens = 0;
    for i = 1:length(fuels)
        o_dens = o_dens + o_densities(i) * o_fracs(i);
    end
        
    %propDens = (OF * o_dens + f_dens) / (OF + 1); %density of combined propellant (kg/m^3)
    propDens = 1758; %kg/m^3, from prop
    %Finding throat area for end of burn
    % / 1000 is to make all constant terms in kg*s/mm, 1st 1-e6 is to make maxPres in Mpa, last 1e6 is to make At in m^2
    A_t = ((Surface_Area(shape, r_max, stage_length) * propDens * C * C_star(1) / 1000) / (g * (maxPres*1e-6)^(1-n)))*1e-6;  %[m^2]
    throatDiameter = 2*sqrt(A_t/pi);    %[m]
    
    
    
    %Set initial values
    T(1) = 0;
    W(1) = r_max;
    P_c(1) = maxPres;
    
    
    %Rocket Optimized for average Chamber Pressure instead of Maximum PSI
    %P_c(1) = 2530375.93;
    
    R_b(1) = (C * (P_c(1) / (1 * 10^6)) ^ n) * 0.001; %Change Pressure unit to MPA and Burn rate to m/s
    propVol = (Area(shape, r_max)-Area(shape, r_min))*stage_length;
    propMass = propVol*propDens;
    f_prop = 1-f_inert;
    totalMass = propMass/(f_prop) + payloadMass;
    inertMass = totalMass-propMass;
    M(1)=inertMass;
    i = 2;
    aIhateTheFrench = 0.0174; %in/s
    nTheEnglishAsWell = 0.04;
    my_waitbar1 = waitbar(0, "Simulation Progress:");
    while W(i-1) > r_min
        %disp(string(W(i-1)/r_max)-innerWidth)
%         disp("Progress: "+string(2*(r_max-W(i-1))*100/r_max)+"%")
        waitbar(1-(W(i-1)-r_min)/(r_max-r_min), my_waitbar1);
        %time step
        T(i) = T(i-1) + dt;
        %CEA Call
        [C_t(i), C_star(i), Isp(i)] = SRM_CEA(P_c(i-1),OF,fuels,f_temps,f_fracs,oxidizers,o_temps,o_fracs,atmoPressure);
        %C_t = [C_t, c_t];
        %C_star = [C_star, c_star];
        %end
        %Chamber pressure
        P_c = [P_c, 0];
        % / 1000 is to make all constant terms in kg*s/mm, 1st 1e6 is to make At in mm^2, last 1e6 is to make P_c in Pa
        P_c(i) = ((Surface_Area(shape, W(i-1), stage_length) * propDens * C * C_star(i) / 1000) / (g * A_t*1e6))^(1/ (1 - n))*1e6; %Pa
        P_c_inStupidImperialUnits(i) = P_c(i) * 0.000145038;
        %Thrust (N)
        %Thrust(i) = C_t(i) * A_t * P_c(i);
        %Burn Rate
        R_b(i) = (C * (P_c(i)*1e-6) ^ n) * 0.001;   %Change Pressure unit to MPA and Burn rate to m/s
        R_b_dumb(i) = (aIhateTheFrench * P_c_inStupidImperialUnits(i)^nTheEnglishAsWell); %in in/s
        %New web distance
        W(i) = W(i-1) - R_b(i) * dt;
 
        %calculate propellant mass
        %Mdot(i) = (Area(shape, W(i-1)) - Area(shape, W(i)))*stage_length*propDens;
        Mdot(i) = Surface_Area(shape, W(i), stage_length)*R_b(i)*propDens;
        M(i) = M(i-1) + Mdot(i)*dt;
        Thrust(i) = Isp(i) * g * Mdot(i);
        
        %increment index
        i = i + 1;
    end
    %disp("Progress: 100%")
    close(my_waitbar1);
    burn_time = T(end);
    accel = Thrust(2:end)./M(2:end);
    

    % Specific Impusle (Isp) [sec]
    avgSpecificImpulse = sum(Isp, 'all')/max(size(Isp));   %max(size(X)) is the same as length(X)
    %%%% Old code (using Riemman Sum): 
    %specificImpulse = trapz(dt,Thrust)/(propMass*g);
  
    % Delta-V [m/s]
    v_e = avgSpecificImpulse * g;  %effective exhaust velocity

    deltaV = v_e * log(M(end) / M(1));
    
    TWR = Thrust./(M.*g);
    %%%% Old code (using Riemman Sum):
    %deltaV = trapz(dt, Thrust(2:end)./M(2:end));
    %deltaV = trapz(dt, Thrust(2:end)) / trapz(dt, M(2:end));
    
    %%Formatting
    %graphs
    figure(1);
    subplot(2,4,1);
    plot(flip(T),Thrust);
    
    title('Thrust');
    grid on;
    
    subplot(2,4,2);
    plot(flip(T),M);
    title('Mass');
    grid on;
    
    subplot(2,4,3);
    plot(flip(T),W);
    title('Web Distance');
    grid on;
    
    subplot(2,4,4);
    plot(flip(T),P_c);
    title('Chamber Pressure');
    grid on;
    
    subplot(2,4,5);
    plot(flip(T),Mdot);
    title('Mass Flow');
    grid on;
    
    subplot(2,4,6);
    plot(flip(T),R_b);
    title('Burn Rate');
    grid on;
    
    subplot(2,4,7);
    plot(flip(T),TWR);
    title('Thrust to Weight Ratio');
    grid on;
    
    subplot(2,4,8);
    plot(flip(T),C_t);
    title('Thrust Coefficient');
    grid on;
    movegui('northwest');
    
    figure(2)
    plot(flip(T), R_b_dumb);
    title('Burn Rate in silly Units');
    xlabel('time (s)');
    ylabel('in/s');
    
    %%flip things for output
    M = flip(M);
    Thrust = flip(Thrust);
    W = flip(W);
    P_c = flip(P_c);
    R_b = flip(R_b);
    Mdot = flip(Mdot);
    TWR = flip(TWR);
    disp(propMass);
    disp(propMass+inertMass);
end

function area = Surface_Area(shape, w, l)
% INPUTS
%shape, the cross sectional shape of the grain
%w, the current web distance
%l, the length of the rocket
    switch(shape)
        %case 'aft finocyl'
        %    area = ;%todo
        case 'circular'
            area = 2 * pi * w * l;
        case 'square'
            area = 8 * w * l;
        otherwise
            fprintf("Invalid shape parameter, see Surface_Area function for valid choices");
            area = 0;
    end
end

function area = Area(shape, w)
% INPUTS
%shape, the cross sectional shape of the grain
%w, the current web distance
    switch(shape)
        %case 'aft finocyl'
        %    area = ;%todo
        case 'circular'
            area = pi * w^2;
        case 'square'
            area = w*w;
        otherwise
            fprintf("Invalid shape parameter, see Surface_Area function for valid choices");
            area = 0;
    end
end

