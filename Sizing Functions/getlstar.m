function [Lstar] = getlstar(fuel, oxidizer)
%% L* DATABASE
%% Header and Information
% Author[s]: Ben Worrell
% Revision: 1
% Revision Date: 4/15/2020
% Purpose: Pulls inputted L* data for a given propellant combination
% sourced from various areas.
% FOR INQUIRIES OR IMPROVEMENTS CONTACT AUTHOR AT bworrell@purdue.edu
%% Sources
% Huzel, D. K., Arbit, H., & Huang, D. H. (1992). Modern engineering for 
% design of liquid-propellant rocket engines. Washington, D.C.: American 
% Institute of Aeronautics and Astronautics.
%% Input[s]
% Input Format: Name [data type]: description [unit]
% 1. Fuel [string]: fuel selection [N/A]
% 2. Oxidizer [string]: oxidizer selection [N/A]
%% Output[s]
% 1. Lstar [float]: Approximated L* value [m]
%% Build Data
propellant_combination = ['RP-1O2(L)' 'CH4(L)O2(L)' 'C2H6O2(L)']; %fuel and oxidizer combination
lstar_data = [1.143 1.5240 1.016]; 
%% Initializations
comb = strcat(fuel,oxidizer);
index = find(propellant_combination == comb);

if ~index
    error('Invalid Fuel/Oxidizer Combination')
else
    Lstar = lstar_data(index);
end
end