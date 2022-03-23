function [] = WritetoExcel(fuel,ox,fuel_oxRatio,length,width,innerwidth,totmass,mass1st,mass2nd,thrust_time,mass_time,thrustweight_time,pressure_time,dV,isp)
% This function outputs the parameters to an excel file named 'Output'
% 
% This is poorly written and needs to be fixed but it should work
% MAKE SURE YOU HAVE EXCEL ON YOURE SYSTEM
% thows an error if you have Output.xlsx open on your system

filename = 'Output.xlsx';

scalars = [length, width, innerwidth, totmass, mass1st, mass2nd, dV, isp, fuel, ox, fuel_oxRatio];

xlswrite(filename, scalars, 1, 'D2:D12');
%xlswrite(filename, thrust_time, 2, '')
end

