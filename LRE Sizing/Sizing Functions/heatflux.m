function [q] = heatflux(T,Tvec,M,rhovec,muvec,gamvec,kvec,Pr,R,D)
%% Heat Flux Calculator
%% Header and Information
% Author[s]: Ben Worrell
% Revision: 1
% Revision Date: 4/28/2020
% Purpose: Calculates the heat flux for an engine using some engine
% parameters and a NASA CEA output. Calculates the recovery temperature and
% arithmetic mean values using the topics from Heister to evaluate the
% gas phase convective heat flux.
% FOR INQUIRIES OR IMPROVEMENTS CONTACT AUTHOR AT bworrell@purdue.edu
%% Source[s]
% 1. Heister, S. D., Anderson, W. E., Pourpoint Thimote?e, & Cassady, R. J.
%       (2019). Rocket propulsion. Cambridge: Cambridge University Press.
%% Input[s]
% Input Format: Name [data type]: description [unit]
% 1. T [float]: Wall Temperature Adjacent to Hot Gas [K]
% 2. Tvec [float]: 3x1 Vector of Temperatures [K]
% 3. M [float]: Mach Number at Evaluation Point [N/A]
% 4. rhovec [float]: 3x1 Vector of Densities [kg/m^3]
% 5. muvec [float]: 3x1 Vector of Kinematic Viscosities [Pa-s]
% 6. gamvec [float]: 3x1 Vector of Speific Heat Ratios [N/A]
% 7. kvec [float]: 3x1 Vector of Thermal Conductivities [W/m-K]
% 8. Pr [float]: Prandtl Number [N/A] 
% 9. R [float]: Gas Constant of Flow of Interest [J/kg-K]
% 10. D [float]: Diameter of Region of Interest [m]

% IMPORTANT: This code only works for one CEA output where indices 1, 2 and
% 3 are the chamber, throat, and nozzle of the engine respectively
%% Output[s]
% 1. q [float]: Heat flux at point of interest [W/m^2]
%% Initializations
% Check if chamber conditions are being evaluated, mach number output for
% CEA is 0, because M << 1. 0.1 Assumed mach number in chamber.
if M  == 0
    M = 0.1;
end
% Evaluate parameters of interest
index = find(Tvec == T); %Find index in CEA outputs where the interest is 
rho = rhovec(index); %Density of interest [kg/m^3]
mu = muvec(index); %Kinematic viscosity of interest [Pa-s]
gam = gamvec(index); %Ratio of specific heats of interest [N/A]
k = kvec(index); %Thermal conductivity of interest [W/m-K]
%% Recovery Temperature and Mean Parameters Calculation
r = Pr^(1/3); %Recovery Factor for a turbulent free boundary layers [N/A]
Tr = T*(1+((gam-1)/2)*r*M^2); %Recovery Temperature [K]
T_am = 0.5*(Tr + T); %Mean Temperature [K]
rho_am = interp1([Tvec(3),Tvec(2),Tvec(1)],[rhovec(3),rhovec(2),rhovec(1)],T_am,'linear','extrap'); %mean density [kg/m^3]
mu_am = interp1([Tvec(3),Tvec(2),Tvec(1)],[muvec(3),muvec(2),muvec(1)],T_am,'linear','extrap'); %mean viscosity [Pa-s]
%% Dimensional Analysis and Heat Transfer Calculation
Re = rho*sqrt(gam*R*T)*D/(mu); % Reynolds Number [N/A]
Nu = 0.026*Re^0.8*Pr^0.8*(rho_am/rho)^0.8*(mu_am/mu)^0.2; %Nusselt Number from Bartz Equation
q = Nu*k/(D) * (Tr-T); %Heat flux [W/m^2]
end