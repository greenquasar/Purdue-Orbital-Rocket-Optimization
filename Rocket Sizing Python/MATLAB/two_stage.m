%Purdue Orbital Fall 2020
%Rocket Sizing Code
%Adapted from python rocket-sizing code
%Authors: Blake Lowe

%Calculate rocket sizing (mass and dimensions) for given inputs
%two stage version

clearvars;
%%CONSTANTS
g = 9.81; %[m/s^2]
%%INPUTS
isp1 = 276; %[s]
isp2 = 276; %[s]
propMassFraction1 = 0.8; %between 0 and 1
propMassFraction2 = 0.8; %between 0 and 1
payloadMass = 15; %[kg]
deltaVtotal = 3711.3; %[m/s]for a single stage
deltaVsplit = 0.5; %between 0 and 1 (1 is all first stage, 0 is all second stage)
finenessRatio = 10; %[dimensionless] length/diameter
propDensity1 = 1252; %[kg/m^3]
propDensity2 = 1252; %[kg/m^3]
%use below two to calculate full rocket length
containerThickness = 0.01; %[m]
engine1Length = 0.5; %[m]
engine2Length = 0.5; %[m]
payloadLength = 0.5; %[m]

%%CALCULATIONS
%mass calculations
deltaV1 = deltaVtotal*deltaVsplit;
deltaV2 = deltaVtotal-deltaV1;
%second stage
massRatio2 = exp(deltaV2/(g*isp2));
propMass2 = payloadMass * ((massRatio2 - 1) / (((1 - massRatio2) / propMassFraction2) + massRatio2)); %[kg]
inertMass2 = (propMass2 / propMassFraction2) - propMass2; %[kg]
totalMass2 = propMass2 + inertMass2 + payloadMass; %[kg]
%first stage
massRatio1 = exp(deltaV1/(g*isp1));
propMass1 = totalMass2 * ((massRatio1 - 1) / (((1 - massRatio1) / propMassFraction1) + massRatio1)); %[kg]
inertMass1 = (propMass1 / propMassFraction1) - propMass1; %[kg]
totalMass = propMass1 + inertMass1 + totalMass2; %[kg]
%dimensioning
propVolume1 = propMass1 / propDensity1; %[m^3]
propVolume2 = propMass2 / propDensity2; %[m^3]
propVolume = propVolume1 + propVolume2;
propDiameter = (4 * propVolume / (pi * finenessRatio))^(1/3); %[m]
propLength = finenessRatio * propDiameter;
diameter = propDiameter + 2 * containerThickness;
length = propLength + engine1Length + engine2Length + payloadLength;
crossSectionalArea = (1/4) * pi * (diameter^2); %[m^2]

%%OUTPUTS
out = [totalMass, totalMass2, length, diameter, crossSectionalArea];
fprintf("Rocket mass: %.2f [kg]\n", totalMass);
fprintf("Second stage mass: %.2f [kg]\n", totalMass2);
fprintf("Rocket Length: %.2f [m]\n", length);
fprintf("Rocket diameter: %.2f [m]\n", diameter);
fprintf("Cross sectional area: %.2f [m^2]\n", crossSectionalArea);
