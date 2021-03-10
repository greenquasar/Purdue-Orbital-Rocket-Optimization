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
        