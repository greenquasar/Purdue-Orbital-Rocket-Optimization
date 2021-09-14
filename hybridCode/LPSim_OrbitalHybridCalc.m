clc
clear

i=1;                                                                       %iterative counter
m_dot_ox = 3;                                                              %usually at 1.75: mass flow rate of oxidizer
r(i) = 1;                                                                  %radius of fuel port
L = 6;                                                                     %length of fuel grain
r_outer = 2;                                                               %maximum radius of fuel port
Ap(i) = pi*r(i)^2;                                                         %area of fuel port as initial radius increases through burn time
p_c(i) = 250;                                                              %chamber pressure (initial set at 375 psi)
o_f(1) = 6; %was 5                                                         %???
N=1;                                                                       %number of fuel ports

g0 = 32.174*12;                                                            %acceleration due to gravity but in in/s^2    
p_a = 14.696;                                                              %psi atmospheric pressure (0 if vacuum, 14.696 at sea level)
R = 4920/32.2;


fuel =  'Paraffin_Wax';                                                    % Fuel choice (match name to CEA inp files!!!)
rho_f = 0.03251456;                                                        % fuel density [lb/in^3]
oxidizer = 'N2O';                                                  % Oxidizer choice (match name to CEA inp files!!!)

[t_c(i), o_f_cea(i), gamma(i), c_star(i), MW] = callCEA1(o_f(i), fuel, oxidizer);

a = .472;                                                                  %[m^(1+2n) kg^(-n) s^(n-1)] regression rate coefficient r=a*Gox^n (empirical data from research)
n = .555;                                                                   %regression rate coefficient r=a*Gox^n (empirical data from research)
a = a * ((3.28*12)^(1+2*n)) * (2.2^(-n))/1000;                             %-> in^(1+2n) lbm^(-n) s^n-1

N=1;
delta_t = .1;

while(r(i)<r_outer)
    per(i) = 2*pi*r(i);                                                    %calculates inner perimeter of fuel grain
    Ab(i) = per(i) *L*N;                                                   %finds total inner surface area of fuel ports
    m_flux_ox(i) = m_dot_ox/(N*(pi*r(i)^2));                               %Oxidizer mass flux lbm/in^2-s
    r_dot(i) = a*(m_flux_ox(i))^n;                                         
    m_dot_fuel(i) = r_dot(i)*rho_f*Ab(i);                                  %regression rate of the fuel grain in/s
    m_dot_total(i) = m_dot_ox + m_dot_fuel(i);                             %total mass flow rate (fuel and oxidizer)
    o_f(i) =m_dot_ox / m_dot_fuel(i);                                      %oxidizer to fuel ratio, changes bc hybrid engine
    
    if(i==1)
        M_e = sqrt((2/(gamma(i) - 1)).*((p_c(1)./p_a).^((gamma(i) - 1)/gamma(i)) - 1));     %mach number
        A_t = (c_star.*m_dot_total(i))./(p_c(i)*g0);    %does this calculate nozzle throat area?
        A_e = (A_t./M_e).*((1 + ((gamma - 1)./2)*M_e.^2)./((gamma + 1)./2)).^((gamma + 1)./(2*(gamma - 1))); %exit throat area?
        epsilon = A_e/A_t;  %expansion (contraction?) ratio
    else
        p_c(i)= (m_dot_total(i-1)*c_star(i-1)/(A_t))/g0; %chamber pressure (I assume)
    end
    
    [CF(i), t_c(i), ISP(i), o_f_cea(i), gamma(i), c_star(i)] = callCEA(p_c(i), o_f(i), fuel, oxidizer, epsilon);
    
    Isp(i) = CF(i)*c_star(i)/g0;
    thrust(i) = CF(i)*p_c(i)*A_t;
    
    i=i+1;
    r(i) = r(i-1) + (r_dot(i-1)) * delta_t;
end
avg_thrust = mean(thrust)
avg_ISP = mean(ISP)
avg_p_c = mean(p_c)

t = 0:delta_t:delta_t*(i-2);

% figure(1)
% plot(t,m_dot_ox);
% plot(t,m_dot_fuel,'linewidth',2.5);
% title('Fuel Mass flows')
% xlabel('Time [s]')
% ylabel('Mass flow [kg/s]')
% grid on
% 
% figure(2)
% plot(t,m_dot_total);
% title('Total Mass flows')
% xlabel('Time [s]')
% ylabel('Mass flow [kg/s]')
% grid on
% 
% 
% figure(3)
% plot(t,thrust,'linewidth',2.5);
% title('Thrust','fontsize',18)
% xlabel('Time [s]','fontsize',12)
% ylabel('Thrust [lbf]','fontsize',12)
% grid on
% 
% figure(4)
% plot(t,Isp,'linewidth',2.5);
% title('ISP','fontsize',18)
% xlabel('Time [s]','fontsize',12)
% ylabel('ISP [s]','fontsize',12)
% grid on
% 
% figure(5)
% plot(t,p_c);
% title('Pressure')
% xlabel('Time [s]')
% ylabel('Pressure [psi]')
% grid on
% 
% figure(7)
% plot(t,r_dot, 'linewidth',2.5);
% title('Regression rate','fontsize',18)
% xlabel('Time [s]','fontsize',12)
% ylabel('Regression rate [in/s]','fontsize',12)
% grid on
% 
% figure(8)
% plot(t,o_f);
% title('O/F')
% xlabel('Time [s]')
% ylabel('O/F')
% grid on

function [t_c, o_f, gamma, c_star, MW] = callCEA1(OtoF, fuel, oxidizer)
    CEA_RUN = true;                                                        % initializes program
    CEA_SAVE_FILE = 'cea_PurdueOrbital.mat';                               % save file

    inp = containers.Map;                                                  % format
    inp('type') = 'eq';                                                    % Sets the type of CEA calculation
    inp('p') = 300:10:400;                                                        % Chamber pressure
    inp('p_unit') = 'psi';                                                 % Chamber pressure units
    inp('o/f') = OtoF;                                                     % Mixture ratio
    inp('fuel') = fuel;                                                    % Fuel name from thermo.inp
    inp('fuel_t') = 298;
    inp('ox') = oxidizer;                                                  % Fuel name from thermo.inp
    inp('ox_t') = 298;                                                     % Ox inlet temperature
    inp('file_name') = 'orbital.inp';                                      % Input/output file name

    %Call the CEA MATLAB code
    if CEA_RUN
        data = cea_rocket_run(inp);
        save(CEA_SAVE_FILE, 'data');
    end

    data_eq = data('eq');
    %combustion temperature
    t_c = squeeze(data_eq('t'));
    t_c = t_c (1,1);
    t_c = t_c *9/5;                                                        %K -> R
    %oxidizer to fuel ratio
    o_f = squeeze(data_eq('o/f'));
    %molecular weight
    MW=squeeze(data_eq('m'));
    MW = MW * 2.205;                                                       %kg->lbs
    %specific heat ratio
    gamma=squeeze(data_eq('gammas'));
    gamma=gamma(1, 1);
    %characteristic velocity
    c_star = squeeze(data_eq('cstar'));
    c_star = c_star(1, 1);
    c_star = c_star * 39.37;                                               % m/s -> in/s
    %chamber pressure
    p_c=squeeze(data_eq('p'))*14.5;
end

function [CF,t_c, ISP, o_f, gamma, c_star] = callCEA(p_c, OtoF, fuel, oxidizer, epsilon)
    CEA_RUN = true;                                                        % initializes program
    CEA_SAVE_FILE = 'cea_PurdueOrbital.mat';                               % save file

    inp = containers.Map;                                                  % format
    inp('type') = 'eq';                                                    % Sets the type of CEA calculation
    inp('p') = p_c;                                                        % Chamber pressure
    inp('p_unit') = 'psi';                                                 % Chamber pressure units
    inp('o/f') = OtoF;                                                     % Mixture ratio
    inp('sup') = epsilon;                                                  % Supersonic area ratios (use sup because we aren't perfectly expanded)
    inp('fuel') = fuel;                                                    % Fuel name from thermo.inp
    inp('fuel_t') = 298;
    inp('ox') = oxidizer;                                                  % Fuel name from thermo.inp
    inp('ox_t') = 298;                                                     % Ox inlet temperature
    inp('file_name') = 'LPSim.inp';                                      % Input/output file name

    %Call the CEA MATLAB code
    if CEA_RUN
        data = cea_rocket_run(inp);
        save(CEA_SAVE_FILE, 'data');
    end

    data_eq = data('eq');
    %combustion temperature
    t_c = squeeze(data_eq('t'));
    t_c = t_c (1,1);
    t_c = t_c *9/5;                                                        %K -> R
    %oxidizer to fuel ratio
    o_f = squeeze(data_eq('o/f'));
    %molecular weight
    %MW=squeeze(data_eq('m'));
    %specific heat ratio
    gamma=squeeze(data_eq('gammas'));
    gamma=gamma(1, 1);
    %characteristic velocity
    c_star = squeeze(data_eq('cstar'));
    c_star = c_star(1, 1);
    c_star = c_star * 3.281*12;                                            % m/s -> in/s
    %Isp
    ISP=squeeze(data_eq('isp'));
    ISP=ISP(2);
    ISP = ISP / 9.81;                                                      % m/s -> s
    %Thrust Coefficient
    CF = squeeze(data_eq('cf'));
    CF = CF(2);
end