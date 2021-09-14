% Purdue Orbital: HTPB and H2O2

CEA_RUN = true;
CEA_SAVE_FILE = 'cea_PurdueOrbital.mat';

inp = containers.Map;
inp('type') = 'eq fr';              % Sets the type of CEA calculation
inp('p') = 20;                % Chamber pressure
inp('p_unit') = 'bar';              % Chamber pressure units
inp('o/f') = 7.6;              % Mixture ratio
%inp('sup') = 180;               % Supersonic area ratios
inp('pip') = 20;                 % Pressure ratios
inp('fuel') = 'HTPB';             % Fuel name from thermo.inp
inp('fuel_t') = 298;                % Fuel inlet temperature
inp('ox') = ["H2O2(L)" "H2O(L)"];     % Ox name from thermo.inpj
inp('ox_wt%') = [0.9 0.1];             % Oxidator weight percentages
inp('ox_t') = [298 298];               % Ox inlet temperature
inp('ox_t_unit') = 'K';
inp('file_name') = 'HTPBH2O2_Percent90_OF_7.6_.inp';% Input/output file name

if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end