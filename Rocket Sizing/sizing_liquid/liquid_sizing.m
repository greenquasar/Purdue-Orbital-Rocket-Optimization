function [mass, volume] = liquid_sizing(isp, propMassFraction, payloadMass, deltaV, propDensity)
    %example call: [mass, volume] = liquid_sizing(270, 0.9, 1, 1000, 1)

    %%DESCRIPTION
    % Calculate rocket sizing (mass and volume)
    % 
    % Purdue Orbital Fall 2021
    % Rocket Sizing Code
    % Adapted from python rocket-sizing code
    % Authors: Blake Lowe, Alex Hanna
    %
    %%INPUTS
    % 1. isp [s]
    % 2. propMassFraction [1], fraction of stage mass which is propellant
    % varies from [0,1]
    % 3. payloadMass [kg]
    % 4. deltaV [m/s]
    % 5. propDensity [kg/m^3], weighted average of oxidizer and fuel based
    % on oxidizer:fuel ratio
    %
    %%OUTPUTS
    % 1. mass [kg], total mass including given payload
    % 2. volume [m^3], propellant volume
    
    %constants
    g = 9.81;

    %mass calculations
    massRatio = exp(deltaV/(g*isp)); %[dimensionless]
    propMass = payloadMass * ((massRatio - 1) / (((1 - massRatio) / propMassFraction) + massRatio)); %[kg]
    inertMass = (propMass / propMassFraction) - propMass; %[kg]
    totalMass = propMass + inertMass + payloadMass; %[kg]

    %dimensioning
    propVolume = totalMass / propDensity; %[m^3]
    
    %outputs
    mass = totalMass;
    volume = propVolume;
end