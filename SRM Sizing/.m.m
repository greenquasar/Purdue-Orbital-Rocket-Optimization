%% solid rocket motor sizing code

%% Constants (SI)
g = 9.81; 
%% Inputs
%time step [seconds]
%grain shape = string ex. 'aft finocyl' must be matched by a case in
%  surface area function
%rocket geometry parameters
% length, case length [m]
% diameter, case diameter [m]
% A_t, throatArea [m]
%propellant characteristics: isp, rho, c_star
%a (burn rate coefficient), 
%n (burn rate exponent)(n<1 (heister232))

%% Functions
function [] = Simulate(dt, grainShape, length, diameter, A_t, rho, c_star, a, n)
    r_max = diameter/2;
    %=======

    %initialized arrays for burn surface area, chamber pressure, rate of burn,
    %web distance and time
    A_b=[]; P_c=[]; R_b=[]; 
    T=[]; W=[];
    %time = 0 and web distance = 0
    T(1) = 0;
    W(1) = 0;

    %calculate initial burn surface area
    A_b(1) = Surface_Area(grainShape, W(1), length);

    %throat area
    %A_t = (s*rho_p*a*c_star)/(g*p_c^(1-n));


    %compute initial chamber pressure, p_c
    P_c(1) = (a*rho*A_b(1)*c_star/(g*A_t))^(1/(1-n)); %todo, fix

    %compute initial burn rate, r_b
    R_b(1) = A_b(1)*P_c(1)^n;

    i = 2;
    while r < r_max
        %time step
        T(i) = T(i-1) + dt;
        %current burn rate
        R_b(i) = a * P_c(i-1)^n;
        %compute new web distance
        W(i) = W(i-1) + R_b(i)*dt;
        %compute updated burn surface area
        A_b(i) = Surface_Area(grainShape, W(i), length);
        %compute updated chamber pressure
        P_c(i) = (a*rho*A_b(i)*c_star/(g*A_t))^(1/(1-n));
        %(optional) compute mass flow, exit mach, exit pressure, exit velocity, thrust

        %(optional) calculate current mass, Δm, ΔV

        %increment index
        i = i + 1;
    end
    %burn time = t
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
            area = 2*pi*w*l;
        case 'square'
            area = 8*w*l;
        otherwise
            fprintf("Invalid shape parameter, see Surface_Area function for valid choices");
            area = 0;
    end
end
function burnTime = burnTime(width, bRate, shape)
%INPUTS
%width, width of propellant grain
%bRate, burn rate of propellant
%shape, geometry of propellant grain
    switch(shape)
        case 'circular'
            burnTime = width / (2 * bRate);
        otherwise
            fprintf("Invalid shape parameter, see burnTime function for valid choices");
            burnTime = 0;
    end
end
function vol = volume(length, width, shape)
%INPUTS
%length, length of propellant grain
%width, width of propellant grain
%shape, geometry of propellant grain
    %get inner diameter from another function
    switch(shape)
        case 'circular'
            vol = (width - innerDiam)^2 * pi * length;
        otherwise
            fprintf("Invalid shape parameter, see volume function for valid choices");
            vol = 0;
    end
end
function [optLength, optWidth] = optimizeDims(minLength, maxLength, ratio, deltaV, bRate, shape, thrust, propDens, dt)
%INPUTS
%minLength and maxLength, the range of lengths we want to optimize for
%ratio, ratio between length and width we want to maintain
%deltaV, the deltaV we require
%bRate, burn rate of propellant
%shape, geometry of propellant grain
%thrust, array of thrust over time
%propDens, density of propellant
%dt, time step
    %dl, amount we increase length width each time
    dl = 0.01;
    %iterate from smallest size to largest size
    for length = minLength:dl:maxLength
        %retrieving width from length and dimensional ratio
        width = length / ratio;
        %getting total burntime from width, burn rate and grain shape
        bTime = burnTime(width, bRate, shape);
        %initializing value of impulse
        impulse = 0;
        %loop to calculate total impulse from thrust array and burntime
        for t = 0:dt:bTime
            %summing impulse from thrust array and timestep
            impulse = impulse + thrust(time) * dt;
        end
        %retrieving volume using length, width, and shape)
        vol = volume(length, width, shape);
        %finding initial mass, m = dv
        mass = vol * propDens;
        %finding deltaV from impulse and total mass
        calc_deltaV = impulse / mass;
        %checking to see if deltaV has been achieved
        if calc_deltaV >= deltaV
            %if achieved, return length and width
            optLength = length;
            optWidth = width;
            return;
        end
    end
end
        