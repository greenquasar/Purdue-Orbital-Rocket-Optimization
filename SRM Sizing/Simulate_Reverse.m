%THIS IS A TEST OF THE NEW WAY WE WILL SIMULATE OUR SRM
%%

%% solid rocket motor sizing code
function [T, W, P_c, Thrust, R_b, burn_time, M, Mdot, A_t, deltaV, specificImpulse, propMass, C_t, C_star] = ...
    Simulate_Reverse(dt, shape, length, width, innerWidth, maxPres, f_inert, atmoPressure, ...
    OF, fuel, f_temp, f_dens, oxidizer, o_temp, o_dens)
    %% Sample Function call
    %    
    
    %% Inputs
    
    %dt: time step (s)
    %shape: circular/square shape: "circular"/"square"
    %length: chamber length (m)
    %width: outer chamber width (m)
    %innerWidth: inner chamber (m)
    %maxPres: maximum chamber pressure stage can withstand (Pa)
    %f_inert: inert mass fraction
    %atmoPressure: atmospheric pressure (Pa)
    %OF: oxidizer to fuel ratio
    %fuel: CEA fuel name, string
    %f_t: fuel inlet temp (K)
    %f_dens: fuel density in solid state (kg/m^3)
    %oxidizer: CEA oxidizer name, string
    %o_t: oxidizer inlet temp (K)
    %o_dens: fuel density in solid state (kg/m^3)
    
    
    %% Outputs
    %T, array of time steps (s)
    %W, array of web distance (m)
    %P_c, array of chamber pressure (Pa)
    %Thrust, array of Thrust (N)
    %R_b, array of burn rate ()
    %burn_time, time to expend stage (s)
    %M, array of mass (kg)
    %Mdot, array of mass flow
    %A_t, throat area required to meet max pressure constraint at end of burn(m^2)
    %deltaV, delta V of stage (m/s)
    %specificImpulse, average of Isp over burn (s)
    %propMass, (kg)
    %C_t, array of coefficient of thrust (1)
    %C_star, array of chamber pressure (Pa)
    
    %% Program
    
    % Constants (SI)
    clc;
    g = 9.80665;    %[m/s^2]

    %GUESSTIMATED DUMMY VALUES FROM CEARUN, get better ones from CEA
    C = 6.2; %C at pressure of 1000 psi (6.895 MPA)
    n = 0.098; %Burn rate exponent at 1000 psi (6.895 MPA)

    C_t = [];
    C_star = [];
    Isp = [];
   
    [C_t(1), C_star(1), Isp(1)] = SRM_CEA(maxPres,OF,fuel,f_temp,oxidizer,o_temp,atmoPressure);
    
    if (width <= innerWidth)
       error("Width is less than or equal to inner width. Impossible!");
    end
    
    %convenience terms by manipulating inputs
    r_max = width / 2;  %[m]
    r_min = innerWidth / 2; %[m]
    propDens = (OF * o_dens + f_dens) / (OF + 1); %density of combined propellant (kg/m^3)
    
    %Finding throat area for end of burn
    A_t = ((Surface_Area(shape, r_max, length) * propDens * C * C_star(1)) / (g * maxPres^(1-n)));  %[m]
    throatDiameter = 2*sqrt(A_t/pi);    %[m]
    
    %Set initial values
    T(1) = 0;
    W(1) = r_max;
    P_c(1) = maxPres;
    R_b(1) = (C * (P_c(1) / (1 * 10^6)) ^ n) * 0.001;   %Change Pressure unit to MPA and Burn rate to m/s
    propVol = (Area(shape, r_max)-Area(shape, r_min))*length;
    propMass = propVol*propDens;
    totalMass = propMass/(1-f_inert);
    inertMass = totalMass-propMass;
    M(1)=inertMass;
    i = 2;
    while W(i-1) > r_min
        disp(string(W(i-1)/r_max))
        %time step
        T(i) = T(i-1) + dt;
        %CEA Call
        [C_t(i), C_star(i), Isp(i)] = SRM_CEA(P_c(i-1),OF,fuel,f_temp,oxidizer,o_temp, atmoPressure);
        %C_t = [C_t, c_t];
        %C_star = [C_star, c_star];
        %end
        %Chamber pressure
        P_c = [P_c, 0];

        P_c(i) = ((Surface_Area(shape, W(i-1), length) * propDens * C * C_star(i)) / (g * A_t))^(1/ (1 - n));

        %Thrust (N)
        Thrust(i) = C_t(i) * A_t * P_c(i); %N %c_t changes over time! this eqn doesn't apply %T=mdot*ve (this assumes ideal nozzle design)
        %Burn Rate
        R_b(i) = (C * (P_c(i) / (1 * 10^6)) ^ n) * 0.001;   %Change Pressure unit to MPA and Burn rate to m/s
        %New web distance
        W(i) = W(i-1) - R_b(i) * dt;
        
       
        %calculate propellant volume
        Mdot(i) = (Area(shape, W(i-1)) - Area(shape, W(i)))*length*propDens;
        M(i) = M(i-1) + Mdot(i);
        %increment index
        i = i + 1;
    end
    
    burn_time = T(end);
    accel = Thrust(2:end)./M(2:end);
    
    % Specific Impusle (Isp) [sec]
    Isp = Isp / g;  %Convert [m/sec] to [sec]
    specificImpulse = sum(Isp, 'all')/max(size(Isp));   %max(size(X)) is the same as length(X)
    %%%% Old code (using Riemman Sum): 
    %specificImpulse = trapz(dt,Thrust)/(propMass*g);
    
    % Delta-V [m/s]
    v_e = specificImpulse * g;  %effective exhaust velocity
    deltaV = v_e * log(M(end) / M(1));
    %%%% Old code (using Riemman Sum):
    %deltaV = trapz(dt, Thrust(2:end)./M(2:end));
    %deltaV = trapz(dt, Thrust(2:end)) / trapz(dt, M(2:end));
    
    %%Formatting
    %graphs
    figure(1);
    subplot(2,4,1);
    plot(-T,Thrust);
    title('Thrust');
    grid on;
    
    subplot(2,4,2);
    plot(-T,M);
    title('Mass');
    grid on;
    
    subplot(2,4,3);
    plot(-T,W);
    title('Web Distance');
    grid on;
    
    subplot(2,4,4);
    plot(-T,P_c);
    title('Chamber Pressure');
    grid on;
    
    subplot(2,4,5);
    plot(-T,Mdot);
    title('Mass Flow');
    grid on;
    
    subplot(2,4,6);
    plot(-T,R_b);
    title('Burn Rate');
    grid on;
    
    subplot(2,4,7);
    plot(-T,C_star);
    title('C Star');
    grid on;
    
    subplot(2,4,8);
    plot(-T,C_t);
    title('Thrust Coefficient');
    grid on;
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
