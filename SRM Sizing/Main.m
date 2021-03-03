%% solid rocket motor sizing code

%% Constants (SI)
g = 9.81; 








%% Inputs
%time step [seconds]
%grain shape = string ex. 'aft finocyl' must be matched by a case in
%  surface area function
%case length [m]
%case diameter [m]
%r_max = case diameter / 2
%propellant characteristics (density probably)

%% Functions
function [] = SRM(dt, grainShape, length, diameter, isp)
r_max = diameter/2;

%start program execution
%create empty arrays for A, p_c, r_b
<<<<<<< Updated upstream
t=0; dt=.0001;
=======
%time = 0 and web distance = 0
T=[]; W=[];
T(1) = 0;
W(1) = 0;

>>>>>>> Stashed changes

%initialized arrays for burn surface area,chamber pressure,
%rate of burn
A=[]; P_c=[]; R_b=[]; 

A(1) = Surface_Area('circular', W(1), length);
%calculate initial burn surface area, A

%throat area
A_t = (s*rho_p*a*c_star)/(g*p_c^(1-n));


%compute initial chamber pressure, p_c
p_c = [a*rho_p*A_b*c_star/(g*A_t)]^(1/(1-n));

%compute initial burn rate, r_b
r_b = a*p_c^n;

while r < r_max
    t+=dt;           %time step
    r_b=a*p_c^n;     %current burn rate
    w+=r_b*dt;       %new web distance
    A_b

%compute new web distance

%compute updated burn surface area

%compute updated chamber pressure

%compute mass flow, exit mach, exit pressure, exit velocity, thrust

%calculate current mass, Δm, ΔV 

%t+=dt

%burn time = t



function area = Surface_Area(shape, w, l)
%shape, the cross sectional shape of the 
    switch(shape)
        %case 'aft finocyl'
        %    area = ;%todo
        case 'circular'
            area = 2*pi*w*l;
        case 'square'
            area = 8*w*l;
        otherwise
            fprintf("Invalid shape parameter, see Surface_Area function for valid choices");
            area = 0;
    end
