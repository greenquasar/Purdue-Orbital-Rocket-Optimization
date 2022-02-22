%function [mass, volume] = liquid_sizing(isp, propMassFraction, payloadMass, deltaV, propDensity)
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
    % 1. Thrust (N)- tbd random value added for now, needs to be changed
    % from lbf to N
    % 2. Fuel Type 
    % 3. Fuel Temp (Kelvin)- can be accessed from thermo.inp
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
    
    %Paths
    %addpath('purdue-orbital-mission-design\LRE Sizing\CEA');
   
    

    %mass calculations
    %massRatio = exp(deltaV/(g*isp)); %[dimensionless]
    %propMass = payloadMass * ((massRatio - 1) / (((1 - massRatio) / propMassFraction) + massRatio)); %[kg]
    %inertMass = (propMass / propMassFraction) - propMass; %[kg]
    %totalMass = propMass + inertMass + payloadMass; %[kg]
 
    %dimensioning
    %propVolume = totalMass / propDensity; %[m^3]
     
    %thrust calculations
    %delta v = 2km/s
    %mass of rocket = total mass
    %totmass*dv = totmass*accel = thrust
    
    %initializing variables
    Thrust = 1000.0; %constant
     Fuel_Types = ['H','BH'];
     Fuel_Temps = [298.150, 200];
%     Fuels = cell(2,2); %defines empty matrix for fuels
%     Fuels{1,1} = 'H'; %defines first fuel (hydrogen)
%     Fuels{1,2} = 298.150; %defines first fuel temperature
%     Fuels{2,1} = 'BH'; %defines second fuel (boron monohydride)
%     Fuels{2,2} = 200; %defines second fuel temperature
     
    %define loop variable for fuels
    i = 1;
    
    %David Tweak
    Ox_Temp = 298.150;
    
    Oxidizers = ['O2-'];
    Ox_Temps = [298.150];
    Chamber_Pressure = 1000; %constant*
    dDiameter = .1;
    Diameter_max = 0.2;
    Diameter_min = .1;
    OF1 = 0;
    OF2 = 1;
    
    %main call loop
    for Oxidizer = Oxidizers
        
        for Fuel = Fuel_Types
            Fuel_Temp = Fuel_Temps(i);
            for Chamber_Diameter = Diameter_min:dDiameter:Diameter_max
                
                [out] = liquid_rocket_design(Thrust, Fuel, Fuel_Temp, Oxidizer, Ox_Temp, Chamber_Pressure, Chamber_Diameter, OF1, OF2);
               
            end
            
            i = i + 1;
          
        end
    end
    % OF1 is set to 0 and OF2 is set to 1 see liquid_rocket _engine_design
    % for more info
    % potential design: have an array containing all fuel types, storage
    % temps, and other data related to fuel and cycle through it storing
    % the results in another array.
    
    %outputs
    %mass = totalMass;
    %volume = propVolume;
%end