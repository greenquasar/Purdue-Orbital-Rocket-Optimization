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
    % 1. Thrust (Newtons)
    % 2. Fuel Type 
    % 3. Fuel Temp (Kelvin)
    % 4. Oxidizer Type
    % 5. Oxidizer Temp (Kelvin)
    % 6. Chamber Pressure (Pa)
    % 7. Gravity
    % 8. Atmospheric Pressure (Pa)
    %
    %%Varying INPUTS
    % 1. Chamber Diameter (Inner Diameter - Meters)
    % 2. Oxidizer Fuel Ratio (low, high; no units)
    % 3. Fineness Ratio (Length/Width)
    %%OUTPUTS
    % 1. OF - Oxidizer Fuel Ratio
    % 2. cstar - Characteristic Temperature
    % 3. isp - Impulse(seconds)
    % 4. m_dot_total - change of total mass (kg/s)
    % 5. m_dot_fuel - change of mass of fuel (kg/s)
    % 6. m_dot_ox - change of mass of oxidizer (kg/s)
    % 7. Dt - time step(s)
    % 8. eps = expansion ratio (chamber pressure/atmospheric pressure)
    % 9. De - exit diameter (inches)
    % 10. Dc - chamber diameter (inches)
    % 11. Lstar - chamber characteristic length (inches)
    % 12. Lc - Chamber Length (inches)
    % 13. qc - heat flux in chamber (Btu/h-ft^2)
    % 14. qt - heat flux a throat of nozzle (Btu/f-t^2)
    
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