function [c_t, C_star, Isp] = SRM_CEA(Pc,OF,fuel,f_t,f_fracs,oxidizer,o_t,o_fracs,atmoPressure)

fclose('all');

addpath('CEA');
% CEA_ROCKET_EXAMPLE: Uses MATLAB CEA wrapper. For in-depth
% documentation read the headers of cea_rocket_run.m,
% cea_rocket_run_single.m, and cea_rocket_read.m

% Change this variable to true to rerun CEA instead of using saved values
CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

atmoPressurePSI = atmoPressure/6895;
PcPSI = Pc/6895;

% The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% hash in C) as input. The dictionary uses MATLAB character arrays as the
% keys, and the value data type varies by which key is used. Details of
% each key are listed in cea_rocket_run.m
% For example: inp('key') = value.
inp = containers.Map;
inp('type') = 'eq';                   % Sets the type of CEA calculation
inp('p') = PcPSI;                        % Chamber pressure
inp('p_unit') = 'psi';                % Chamber pressure units
inp('o/f') = OF;                      % Mixture ratio
%inp('sup') = 70;                    % Supersonic area ratios
%14.7 psi is standard 
%1.45 is 18km
atmoPressurePSI = atmoPressure/6895;
PcPSI = Pc/6895;
inp('pip') = PcPSI/atmoPressurePSI;           % Pressure ratios
inp('fuel') = fuel;               % Fuel name from thermo.inp
inp('fuel_t') = f_t;               % Fuel inlet temperature
inp('fuel_wt%') = f_fracs;         % Fuel mass fraction
inp('ox') = oxidizer;                  % Ox name from thermo.inp (O2(L))
inp('ox_t') = o_t;                  % Ox inlet temperature
inp('ox_wt%') = o_fracs;            % Ox weight percentage
inp('file_name') = 'SRM_CEA.inp';   % Input/output file name

if CEA_RUN
    data = cea_rocket_run(inp);       % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

% The output data structure, called 'data' in this case, is also a MATLAB
% map. 'data' contains a single entry for each of the CEA calculation types
% listed ('eq' and 'fr'). For instance, if only 'fr' is listed, then 'data'
% will only contain a single entry under data('fr').
data_eq = data('eq');

c_f_vec = squeeze(data_eq('cf')); %cf is coefficient of thrust at the exit pressure (not desired)
c_f = c_f_vec(2);

c_t = c_f; % wrong for now

%look at line 84 of lre_design.m
%set pa from atmo pressure, set pe as 
%get Mach, gammas
%calculate epsilon
%calculate ct

C_star_vec = squeeze(data_eq('cstar'));
C_star = C_star_vec(2);
isp_vec = squeeze(data_eq('isp'));
Isp = isp_vec(2)/9.80665; % convert m/s to s

%uncomment below to see list of all terms cea outputs
%keys(data_eq)

%%%%%%%%%%%%%%% keys(data_eq) %THIS LISTS ALL OUTPUTS %%%%%%%%%%%%%%%%%%%%%
% Use keys(data_eq) or keys(data_fr) to see the contents of each map
% respectively. Every output of CEA is contained in these keys, including
% molar concentrations. Most keys contain a 3D array with columns
% corresponding to the pressure, O/F, and area/pressure ratio inputs
% respectively. If only a single value is given for one of these inputs,
% the output will still be a 3D array. The squeeze() MATLAB function must
% be used to reduce the number of dimensions appropriately. Read the notes
% at the top of cea_rocket_read.m for more details.

end