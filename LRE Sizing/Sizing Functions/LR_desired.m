fclose('all');

%% To Do
% 1)Make it so that depending on what fuel and oxidizers are given, script
% finds correct CEA Temp to pass through
% 2)iterate through chamber diameter to get maximum thrust (to get desired
% delta v
% 3)output mass
% 4) make Lstar function

%% Provided Bounds
% 1)dv = 2 km/s
% 2)Pc = [1000 3000] psi
% 3)fineness ratio 10


% input values: (Thrust, Fuel, Fuel_Temp, Oxidizer, Ox_Temp, Chamber_Pressure, Chamber_Diameter, OF1, OF2)
% provides [OF cstar Isp m_dot_total m_dot_fuel m_dot_ox Dt eps De Dc Lstar Lc q_c q_t]

% CEA call: [cstar, isp, M, gamma, T, rho, mu, Pr, Mw, k] = EngineCEA(Pc,OF,fuel,f_t,oxidizer,o_t)
% CEA first call initializations

Thrust = 5000; %[N]
Fuel = 'RP-1'; %'RP-1'
Fuel_Temp = 298.150; %[K] 298.150 for RP-1
Oxidixer = 'O2(L)'; %'O2(L)'
Ox_Temp = 90.170; %[K] 90.170 for O2(L)
Chamber_Pressure = 14000000; %[Pa]
Chamber_Diameter = .5; %[m]
OF1 = 0.1; %OF lower bound
OF2 = 5; %OF upper bound

[out] = liquid_rocket_engine_design(Thrust, Fuel, Fuel_Temp, Oxidizer, Ox_Temp, Chamber_Pressure, Chamber_Diameter, OF1, OF2);

%--------------------------------------
function [out] = liquid_rocket_engine_design(Thrust, Fuel, Fuel_Temp, Oxidizer, Ox_Temp, Chamber_Pressure, Chamber_Diameter, OF1, OF2)
%example input:[out] = liquid_rocket_engine_design(5000, 'RP-1', 298.150, 'O2(L)', 90.170, 14000000, 0.5, 0.1, 5) 

%% Liquid Propulsion System Design Code V2(metric units)
%% Header and Information
% Author[s]:Ben Worrell
% Editor[s]: Michael Berthin, Alex Hanna
% Version 1
% Last Update Date: 10/25/2021
% Differences from original: uses only metric units,
% Revision Date
% Purpose: Facilitates liquid rocket engine design based on given 
% thrust, propellant combination, steady state chamber pressure, and O/F
% ratio. Assumes an optimally expanded 80% bell nozzle contour and analyzes
% only steady state parameters. Uses NASA Chemical Equilibrium with
% Applications [CEA] to pull combustion data from which performance
% chacteristics are derived. Engine geometry determination then follows.
% FOR INQUIRIES OR IMPROVEMENTS CONTACT AUTHOR AT bworrell@purdue.edu
%% Sources
% 1. Heister, S. D., Anderson, W. E., Pourpoint Thimote?e, & Cassady, R. J.
%       (2019). Rocket propulsion. Cambridge: Cambridge University Press.
% 2. Huzel, D. K., Arbit, H., & Huang, D. H. (1992). Modern engineering for
%       design of liquid-propellant rocket engines. Washington, D.C.: 
%       American Institute of Aeronautics and Astronautics.
%% Input[s]
% Input Format: Name [data type]: description [unit]
% 1. Thrust [float]: Desired Thrust Level [N]
% 2. Fuel [string]: Fuel Selection
% 3. Fuel_Temp [float]: Fuel Storage Temperature [K]
% 4. Oxidizer [string]: Oxidizer Selection
% 5. Ox_Temp [float]: Oxidizer Storage Temperature [K]
% 6. Chamber_Pressure [float]: Steady State, Average Chamber Pressure [psi]
% 7. Chamber_Diameter [float]: Chamber Diameter [in]
% 8. OF1 [float]: Low initial guess for OF
% 9. OF2 [float]: High initial guess for OF

% IMPORTANT: For any storable propellants, CEA only accepts the temperature
% 298.15 K, this input is for cyrogenic propellants
%% Output[s]
% Output is 1xn vector 'out' components of this vector discussed below
% 1. O/F [float]: Optimized O/F ratio based on performance [N/A]
% 2. cstar [float]: Characteristic velocity [ft/s]
% 3. Isp [float]: Specific impulse [s]
% 4. m_dot_tot [float]: Total mass flow rate [lbm/s]
% 5. m_dot_fuel [float]: Fuel mass flow rate [lbm/s]
% 6. m_dot_ox [float]: Oxidizer mass flow rate [lbm/s]
% 7. Dt [float]: Throat diameter [in]
% 8. eps [float]: Expansion ratio [N/A]
% 9. De [float]: Exit diameter [in]
% 10. Dc [float]: Chamber Diameter [in]
% 11. Lstar [float]: Chamber Characteristic Length [in]
% 12. Lc [float]: Chamber Length [in]
% 13. qc [float]: Heat flux in chamber [Btu/h-ft^2]
% 14. qt [float]: Heat flux at throat of nozzle [Btu/f-ft^2]


%% Add CEA to File Path
% IMPORTANT: If you are using this on a computer for the first time, you must download the CEA files and build the filepath for your own computer
fclose('all');
addpath('X:\Documents\Purdue\Orbital\purdue-orbital-mission-design\LRE Sizing\CEA', '-end');
savepath();


%% Important Notes [this script]
% 1) Deleted unit conversions to maintain SI units - reference original liquid design code for conversions

%% Important Notes [EngineCEA]
% 1) Modify atmoPressure to pressure at desired altitude
% 2) EngineCEA outputs cstar in ft/s (internal conversion), convert it back
% after initial EngineCEA call
% 3) Chamber pressure passed to EngineCEA in psi
% 4) Might want to make inp('fuel_t_unit') = 'K', not sure what default is
% 5) Might want to make inp('ox_t_unit') = 'K', not sure what default is
% 6) All outputs should be in SI units

%% Important Notes [cea_rocket_run]
% 1) Input pressure default bar (changed to psi in EngineCEA)
% 2) cea_rocket_read says all properties in SI units (kg, m, s, K)

%% Initializations
g = 9.81; % acceleleration due to gravity [m/s^2]
Pc = Chamber_Pressure; % [Pa]
pc = Pc / 6895; % [psi] passed to EngineCEA
Dc = Chamber_Diameter; % [m]

pe = 14.7; % exit pressure [psi]
pa = 14.7; % atmospheric pressure [psi]
effcstar = 0.95; % c* efficiency [N/A]
effcf = 0.9; % cf efficiency [N/A]

Ru = 8314; % [J/kmol*K]
F = Thrust; % [N]

%% O/F Ratio Optimization
% This will take a while (anywhere between 7 and 20 minutes)
OF = OF1:0.1:OF2; %vector of O/F 
Ispvec = zeros(1,length(OF)); %vector of corresponding Isp
for j = 1:length(OF)
    [cstar, ~, Mach, gam, T, rho, mu, Pr, Mw, k]= EngineCEA(pc,OF(j),Fuel,Fuel_Temp,Oxidizer,Ox_Temp);
    % Formatting
    % 1 -> Chamber Condition
    % 2 -> Throat Condition
    % 3 -> Nozzle Condition
    epsilon = 1./Mach(3) .* ((2+(gam(3)-1).*Mach(3).^2)./(gam(3)+1)).^((gam(3)+1)./(2.*(gam(3)-1)));
    cf = sqrt(((2.*gam(3,1).^2)./(gam(3,1)-1)) .* (2./(gam(3,1)+1)).^((gam(3,1)+1)./(gam(3,1)-1)) .* (1-(pe/pc).^((gam(3)-1)./gam(3)))) + (pe/pc - pa/pc).*epsilon; % Dimensionless thrust
    cstar = cstar / 3.281; %convert to [m/s]
    Ispvec(j) = (cf.*effcf.*cstar.*effcstar)./g; %[s]
end
OFopt = mean(OF(Ispvec == max(Ispvec))); %collects O/F that yields highest Isp
OF = OFopt; %simplifies O/F argument

%% Combustion Performance Analysis
[~, ~, Mach, gam, T, rho, mu, Pr, Mw, k] = EngineCEA(pc,OF,Fuel,Fuel_Temp,Oxidizer,Ox_Temp);
cstar = sqrt((Ru*T(1)/(gam(1)*Mw(1)))*((2/(gam(1)+1))^(-(gam(1)+1)/(gam(1)-1)))); %characteristic velocity [m/s]
eps = (1/Mach(3))*((2+(gam(3)-1)*(Mach(3)^2))/(gam(3)+1))^((gam(3)+1)/(2*(gam(3)-1))); %nozzle expansion ratio [N/A]
cf = sqrt(((2*gam(3)^2)/(gam(3)-1))*((2/(gam(3)+1))^((gam(3)+1)/(gam(3)-1)))*(1-(pe/pc)^((gam(3)-1)/gam(3))))+(pe/pc - pa/pc)*eps; %thrust coefficient [N/A]
Isp = (cf*cstar*effcstar*effcf)/g; %specific impulse [s]
ve = sqrt((2*gam(1)*Ru*T(1))/(Mw(1)*(gam(1)-1))*(1-(pe/pc)^((gam(1)-1)/gam(1)))); %nozzle exit velocity [m/s]
m_dot_tot = F/ve; %steady total mass flow [kg/s]
m_dot_f = m_dot_tot*(1/(OF+1)); %steady fuel mass flow [kg/s]
m_dot_ox = m_dot_f*OF; %steady oxidizer mass flow [kg/s]

%% Engine Geometry Determination
% Chamber
At = m_dot_tot*cstar*effcstar/(pc*6894.76); %throat area [m^2]
Dt = sqrt(4*At/pi); %throat diameter [m]
De = sqrt((Dt^2)*eps); %exit area [m^2]
chi = Dc^2/Dt^2; %contraction ratio [N/A]

Lstar = 1.143; %characteristic length [m](45 in)

% Uncomment this line once getlstar.m function is fully built at this point
% 45 in is a relatively good estimate for a lot of these propellants
%Lstar = getlstar(Fuel, Oxidizer); %chamber characteristic length [m]
Vc = Lstar*At; %chamber volume [m^3]
Lc = Vc/((pi/4)*Dc^2); %chamber length [m]

%% Nozzle Contour
[x,y] = nozzle(Dt*39.3701/2,eps,chi,Lc*39.3701); %plots nozzle contour [m] %39.3701 in/m

%% Thermal Analysis
% Chamber
q_c = heatflux(T(1),T,Mach(1),rho,mu,gam,k,Pr(1),Ru/Mw(1),Dc);
% Throat
q_t = heatflux(T(2),T,Mach(2),rho,mu,gam,k,Pr(2),Ru/Mw(2),Dt);

%% Build Output Vector
out = [OF cstar Isp m_dot_total m_dot_fuel m_dot_ox Dt eps De Dc Lstar Lc q_c q_t];
end
