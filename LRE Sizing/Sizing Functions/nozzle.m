function [x, y] = nozzle(Rt, eps, chi,Lc)
%% Nozzle Profile Generation Code
%% Header and Information
% Author[s]: Ben Worrell
% Revision: 1
% Revision Date: 4/15/2020
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
%% Initializations
thd = 37.5; %divergence angle [deg]
the = 7; %convergence angle [deg]
%% From Chamber to Throat Enterance Profile
% Initialize x and y plotting
theta = linspace(-130,-90,100);
x_conv = 1.5.*Rt.*cosd(theta);
y_conv = 1.5.*Rt.*sind(theta)+1.5.*Rt+Rt;
%% From Throat Enterance Profile to Beginning of Nozzle Contour
theta = linspace(-90,(thd-90),100);
x_div = 0.382.*Rt.*cosd(theta);
y_div = 0.382.*Rt.*sind(theta)+0.382.*Rt+Rt;
%% Determination of Nozzle Contour Bezier Curve Properties
Nx = x_div(end); %x coordinate for beginning of nozzle
Ny = y_div(end); %y coordinate for beginning of nozzle
Ex = 0.8*((sqrt(eps)-1)*Rt/tand(15)); %x coordinate of exit plane
Ey = sqrt(eps)*Rt; %y coordinate of exit plane
m1 = tand(thd); %slope of line tangent to nozzle start
m2 = tand(the); %slope of line tangent to exit plane
c1 = Ny - m1*Nx; %equation for line tangent to nozzle start
c2 = Ey - m2*Ex; %equation for line tangent to nozzle exit
Qx = (c2-c1)/(m1-m2); %x coordinate of point where two lines intersect
Qy = (m1*c2-m2*c1)/(m1-m2); %y coordinate of point where two line smeet
%% Nozzle Contour
t = linspace(0,1);
%quadratic bezier curves for nozzle contour
x_noz = (((1-t).^2).*Nx) + (2.*(1-t).*t.*Qx) + ((t.^2).*Ex);
y_noz = (((1-t).^2).*Ny) + (2.*(1-t).*t.*Qy) + ((t.^2).*Ey);
%% Determination of Chamber Contour Bezier Curve Properties
Cx = x_conv(1); %x coordinate of where throat contour begins
Cy = y_conv(1); %y coordinate of where throat contour begins
L_conv = 0.8*((((sqrt(chi)-1)*Rt))/tand(45)); %length of the converging portion of the nozzle
Chamx = Cx - L_conv; %x coordinate of chamber beginning
Chamy = sqrt(chi)*Rt; %y coordinate of chamber beginning
m1 = tand(-45);
m2 = tand(180);
c1 = Cy - m1*Cx;
c2 = Chamy - m2*Chamx;
Px = (c2 - c1)/(m1 - m2);
Py = (m1*c2 - m2*c1)/(m1-m2);
%% Chamber Contour
x_cham = (((1-t).^2).*Cx) + (2.*(1-t).*t.*Px) + ((t.^2).*Chamx);
y_cham = (((1-t).^2).*Cy) + (2.*(1-t).*t.*Py) + ((t.^2).*Chamy);
x_chamber = [(Chamx - Lc) Chamx]; %creating x values of chamber line
Rc = sqrt(chi)*Rt;
y_chamber = [Rc Rc]; %creating y values of chamber line
%% Build Contour Vector
% x = [x_chamber x_cham x_conv x_div x_noz];
% y = [y_chamber y_cham y_conv y_div y_noz];
x = [x_chamber x_conv x_div x_noz];
y = [y_chamber y_conv y_div y_noz];
%% Plot
figure
plot(x,y,'k',x,-y,'k')
title('80% Bell Nozzle Contour')
xlabel('x')
ylabel('y')
grid on
axis equal
end