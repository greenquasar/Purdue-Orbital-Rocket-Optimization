%solid rocket motor sizing code

%inputs
%time step [seconds]
%grain shape = string ex. 'aft finocyl' must be matched by a case in
%  surface area function
%case length [m]
%case diameter [m]
%r_max = case diameter / 2
%propellant characteristics (density probably)

function [] = SRM(dt, etc)

%start program execution
%create empty arrays for A, p_c, r_b
%time = 0 and web distance = 0


%initialized arrays for burn surface area,chamber pressure,
%rate of burn
A=[];p_c=[];r_b=[]; 

%calculate initial burn surface area, A
%assumed to be circular


%compute initial chamber pressure, p_c

%compute initial burn rate, r_b

%while r < r_max

%compute new web distance

%compute updated burn surface area

%compute updated chamber pressure

%compute mass flow, exit mach, exit pressure, exit velocity, thrust

%calculate current mass, Δm, ΔV 

%t+=dt

%burn time = t



function [area] = Surface_Area(w, shape)
    switch(shape)
        case 'aft finocyl'
            area = 1;%todo
        case 'circular'
            area = 1;%todo
    end
    area = 0;
