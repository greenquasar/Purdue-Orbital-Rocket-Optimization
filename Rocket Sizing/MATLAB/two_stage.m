%Purdue Orbital Fall 2020
%Rocket Sizing Code
%Adapted from python rocket-sizing code
%Authors: Blake Lowe

%Calculate rocket sizing (mass and dimensions) for given inputs
%two stage version
function [out] = two_stage(isp1, isp2, propMassFraction1, propMassFraction2, payloadMass, deltaVtotal, deltaVsplit, finenessRatio, propDensity1, propDensity2, containerThickness, engine1Length, engine2Length, payloadLength)  
% sample call -- [] = two_stage(276, 276, 0.8, 0.8, 1.25, 8000, 0.5, 10, 1252, 1252, 0.01, 0.5, 0.5, 0.5)
clearvars;

%%INPUTS
% 1. isp1 (s)
% 2. isp2 (s)
% 3. propMassFraction1
% 4. propMassFraction2
% 5. payloadMass (kg)
% 6. deltaVtotal (m/s)
% 7. deltaVsplit
% 8. finenessRatio
% 9. propDensity1 (kg/m^3)
% 10. propDensity2 (kg/m^3
% 11. containerThickness (m)
% 12. engine1Length (m)
% 13. engine2Length (m)
% 14. payloadLength (m)

%%OUTPUTS
% 1. totalMass = total rocket mass (kg)
% 2. totalMass2 = second stage mass (kg)
% 3. length = rocket length (m)
% 4. diamter = rocket diamter (m)
% 5. crossSectionalArea = cross sctional area (m^2)

%%CONSTANTS
g = 9.81; %[m/s^2]

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
