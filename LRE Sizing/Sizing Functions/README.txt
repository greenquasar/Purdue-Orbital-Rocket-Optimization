%LIQUID ROCKET ENGINE SIZING CODE
%Author: Ben Worrell [bworrell@purdue.edu]
%Date: 4/28/2020
%Revision: 1
%Requires NASA CEA Matlab Wrapper
%Function: liquid_rocket_engine_design.m
%Sub-Functions: EngineCEA.m, nozzle.m, and heatflux.m

%% DESCRIPTION OF ENGINE DESIGN FUNCTION
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
% 1. Thrust [float]: Desired Thrust Level
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
% IMPORTANT: If you are using this on a computer for the first time,
% you must download the CEA files and build the filepath for your own
% computer

%% DESCRIPTION OF NOZZLE DESIGN FUNCTION (nozzle.m)
% Purpose: Uses a quadratic bezier curve to generate a parabolic profile
% for an 80% bell nozzle contour.
% FOR INQUIRIES OR IMPROVEMENTS CONTACT AUTHOR AT bworrell@purdue.edu
%% Sources
% 1. Heister, S. D., Anderson, W. E., Pourpoint Thimote?e, & Cassady, R. J.
%       (2019). Rocket propulsion. Cambridge: Cambridge University Press.
% 2. Huzel, D. K., Arbit, H., & Huang, D. H. (1992). Modern engineering for 
%       design of liquid-propellant rocket engines. Washington, D.C.: American 
%       Institute of Aeronautics and Astronautics.
% 3. Newlands, R. (n.d.). The Thrust Optimized Parabolic Nozzle.
% 4. Sutton, G. P., & Biblarz, O. (2017). Rocket propulsion elements. Hoboken, NJ: Wiley.
%% Input[s]
% Input Format: Name [data type]: description [unit]
% 1. Rt [float]: throat radius [m]/[in]
% 2. eps [float]: expansion ratio [N/A]
% 3. chi [float]: contraction ratio [N/A]
% 4. Lc [float]: chamber length [m]/[in]
%% Output[s]
% 1. x: generated contour profile x data
% 2. y: generated contour profile y data

%% DESCRIPTION OF HEAT FLUX FUNCTION (heatflux.m)
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