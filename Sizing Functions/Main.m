%Purdue Orbital Fall 2020
%Rocket Sizing Code
%Adapted from python rocket-sizing code
%Authors: Blake Lowe

%Calculate rocket sizing (mass and dimensions) for given inputs
%single stage version

clearvars;
%CONSTANTS
g = 9.81; %[m/s^2]
%INPUTS
isp = 300; %[s]
propMassFraction = 0.8; %between 0 and 1
payloadMass = 2.5; %[kg]
deltaV = 9804.9; %[m/s]for a single stage
finenessRatio = 15; %[dimensionless] length/diameter
propDensity = 1; %[kg/m^3]
containerThickness = 0; %[m] set to 0 to ignore
%use below two to calculate full rocket length
%engineLength = 0.5; %[m]
%payloadLength = 0.5; %[m]

%CALCULATIONS
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