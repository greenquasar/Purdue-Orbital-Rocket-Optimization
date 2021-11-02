function [mass, volume] = liquid_sizing(isp, propMassFraction, payloadMass, deltaV, propDensity);
    %mass calculations
    massRatio = exp(deltaV/(g*isp));
    propMass = payloadMass * ((massRatio - 1) / (((1 - massRatio) / propMassFraction) + massRatio)); %[kg]
    inertMass = (propMass / propMassFraction) - propMass; %[kg]
    totalMass = propMass + inertMass + payloadMass; %[kg]

    %dimensioning
    propVolume = totalMass / propDensity; %[m^3]
end