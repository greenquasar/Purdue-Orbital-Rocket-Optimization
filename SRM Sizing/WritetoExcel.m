function [] = WritetoExcel(fuel,ox,fuel_oxRatio,length,width,innerwidth,totmass,mass1st,mass2nd,thrust_time,mass_time,thrustweight_time,pressure_time,dV,isp)
% This function outputs the parameters to an excel file named 'Output'
% 
% This is poorly written and needs to be fixed but it should work
% MAKE SURE YOU HAVE EXCEL ON YOURE SYSTEM
filename = 'Output.xlsx';
xlswrite(filename, length, 1, 'D2');
xlswrite(filename, width, 1, 'D3');
xlswrite(filename, innerwidth, 1, 'D4');
xlswrite(filename, totmass, 1, 'D5');
xlswrite(filename, mass1st, 1, 'D6');
xlswrite(filename, mass2nd, 1, 'D7');
xlswrite(filename, dV, 1, 'D8');
xlswrite(filename, isp, 1, 'D9');
xlswrite(filename, fuel, 1, 'D10');
xlswrite(filename, ox, 1, 'D11');
xlswrite(filename, fuel_oxRatio, 1, 'D12');
xlswrite(filename, thrust_time, 2, '')
end

