%THIS IS A TEST OF THE NEW WAY WE WILL SIMULATE OUR SRM
%%

%% solid rocket motor sizing code
function [T, W, P_c, Thrust, R_b, burn_time, M, Mdot, A_t, deltaV, specificImpulse, propMass] = Simulate_Reverse(dt, shape, length, width, innerWidth, maxPres, f_inert, atmoPressure,OF,fuel,f_t,oxidizer,o_t)
    %% Constants (SI)
    g = 9.80665; 
    %% Inputs
    %time step (s)
    %shape: circular/square shape: "circular"/"square"
    %length: chamber length (m)
    %width: outer chamber width (m)
    %innerWidth: inner chamber (m)
    %maxPres: maximum chamber pressure desired (Pa)
    %f_inert: inert mass fraction
    %atmoPressure: atmospheric pressure (Pa)
    %OF: oxygen to fuel ratio
    %fuel: CEA fuel name
    %f_t: fuel inlet temp (K)
    %oxidizer: CEA oxidizer name
    %o_t: oxygen inlet temp (K)
    
    %rocket geometry parameters: Width, inner width, shape, and length [all in meters]
    %maximum pressure desired [Pa]

    %GUESSTIMATED DUMMY VALUES FROM CEARUN, get better ones from CEA
    C = 6.2; %C at pressure of 1000 psi (6.895 MPA)
    n = 0.098; %Burn rate exponent at 1000 psi (6.895 MPA)
    propDens = 1500; %Average density for 70% AP-HTPB (kg/m^3) %maybe use rho from cea?
    
    C_t = [];
    C_star = [];
   
    [c_t, c_star] = SRM_CEA(maxPres,OF,fuel,f_t,oxidizer,o_t, atmoPressure);
    C_t = [C_t, c_t];
    C_star = [C_star, c_star];

    if (width <= innerWidth)
       error("Width is less than or equal to inner width. Impossible!");
    end

    r_max = width / 2;
    r_min = innerWidth / 2;
    %finding throat area for end of burn
    A_t = ((Surface_Area(shape, r_max, length) * propDens * C * c_star) / (g * maxPres^(1-n)));
    throatDiameter = 2*sqrt(A_t/pi); %m
    %set initial values
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
        %if(we gon call it?)
        [c_t, c_star] = SRM_CEA(P_c(i-1),OF,fuel,f_t,oxidizer,o_t, atmoPressure);
        C_t = [C_t, c_t];
        C_star = [C_star, c_star];
        %end
        %Chamber pressure
        P_c = [P_c, 0];

        P_c(i) = ((Surface_Area(shape, W(i-1), length) * propDens * C * c_star) / (g * A_t))^(1/ (1 - n));

        %Thrust (N)
        Thrust(i) = c_t * A_t * P_c(i); %N %c_t changes over time! this eqn doesn't apply %T=mdot*ve (this assumes ideal nozzle design)
        %Burn Rate
        R_b(i) = (C * (P_c(i) / (1 * 10^6)) ^ n) * 0.001;   %Change Pressure unit to MPA and Burn rate to m/s
        %New web distance
        W(i) = W(i-1) - R_b(i) * dt;
        
        
        %calculate propellant volume
        Mdot(i) = (Area(shape, W(i)) + Area(shape, W(i-1)))*length*propDens;
        M(i) = M(i-1) + Mdot(i);
        %increment index
        i = i + 1;
    end
    burn_time = T(end);
    accel = Thrust(2:end)./M(2:end);
    deltaV = trapz(dt, Thrust(2:end)./M(2:end));
    
    v_e = 2000;
    deltaVcheck = v_e*log(M(1)-M(end)); %need v_e
    specificImpulse = trapz(dt,Thrust)/(propMass*g);
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
