function [mass, volume] = liquid_sizing(isp, propMassFraction, payloadMass, deltaV, finenessRatio, propDensity, containerThickness);
    %mass calculations
    massRatio = exp(deltaV/(g*isp));
    propMass = payloadMass * ((massRatio - 1) / (((1 - massRatio) / propMassFraction) + massRatio)); %[kg]
    inertMass = (propMass / propMassFraction) - propMass; %[kg]
    totalMass = propMass + inertMass + payloadMass; %[kg]

    %dimensioning
    propVolume = totalMass / propDensity; %[m^3]
    diameter = (4 * propVolume / (pi * finenessRatio))^(1/3); %[m]
    length = diameter * finenessRatio; %[m]
    diameter = diameter + 2 * containerThickness;
    crossSectionalArea = (1/4) * pi * (diameter^2); %[m^2]
end