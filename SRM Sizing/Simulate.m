%% solid rocket motor sizing code
function [T, W, P_c, Thrust, R_b, burn_time, A_t] = Simulate(dt, shape, length, width, innerWidth, maxPres)
%% Constants (SI)
g = 9.81; 
%% Inputs
%time step [seconds]
%rocket geometry parameters: Width, inner width, shape, and length [all in meters]
%maximum pressure desired [Pa?]

%GUESSTIMATED DUMMY VALUES FROM CEARUN, get better ones from CEA
C = 6.2; %C at pressure of 1000 psi (6.895 MPA)
n = 0.098; %Burn rate exponent at 1000 psi (6.895 MPA)
C_star = 1077.8; %C* for 70% AP-HTPB 
propDens = 1500; %Average density for 70% AP-HTPB (kg/m^3)
c_t = 0.5; %units? who needs those

    r_max = width / 2;
    r_min = innerWidth / 2;
    %finding throat area for end of burn
    A_t = (Surface_Area(shape, r_min, length) * propDens * C * C_star / (g * maxPres^(1-n)));
    %set initial values
    T(1) = 0;
    W(1) = r_min;
    P_c(1) = ((Surface_Area(shape, W(1), length) * propDens * C * C_star) / (g * A_t))^(1/ 1 - n);
    R_b(1) = C * P_c(1)^n;
    
    i = 2;
    while W(i-1) < r_max
        %time step
        T(i) = T(i-1) + dt;
        %Chamber pressure
        P_c = [P_c, 0];
        P_c(i) = ((Surface_Area(shape, W(i-1), length) * propDens * C * C_star) / (g * A_t))^(1/ 1 - n);
        %Thrust
        Thrust(i) = c_t * A_t * P_c(i);
        %Burn Rate
        R_b(i) = C * P_c(i)^n;
        %New web distance
        W(i) = W(i-1) + R_b(i) * dt;
        %increment index
        i = i + 1;
    end
    burn_time = T(i-1);
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