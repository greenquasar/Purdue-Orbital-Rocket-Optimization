function [] = WritetoExcel(fuel,ox,fuel_oxRatio,length,width,innerwidth,totmass,mass1st,mass2nd,thrust_time,mass_time,thrustweight_time,pressure_time,time,dV,isp)
    % This function outputs the parameters to an excel file named 'InputOutput.xlsx'
    % 
    % This is poorly written and needs to be fixed but it should work
    % MAKE SURE YOU HAVE EXCEL ON YOUR SYSTEM
    % throws an error if you have Output.xlsx open on your system
    %
    % Make sure that any vector input is a VERTICAL vector, otherwise it will
    % not print properly.

    filename = 'InputOutput.xlsx';

    scalars = [length, width, innerwidth, totmass, mass1st, mass2nd, dV, isp, fuel, ox, fuel_oxRatio]';

    %prints out the values in the scalar vector above
    writematrix(scalars,filename,'Sheet',1,'range','D2:D12')

    %prints out the data for the thrust over time graph
    writematrix(time,filename,'Sheet',2,'range','A3:A77')
    writematrix(thrust_time,filename,'Sheet',2,'range','B3:B77')

    %prints out the data for the mass over time graph
    writematrix(time,filename,'Sheet',2,'range','D3:D77')
    writematrix(mass_time,filename,'Sheet',2,'range','E3:E77')

    %prints out the data for the thrust to weight ratio over time graph
    writematrix(time,filename,'Sheet',2,'range','G3:G77')
    writematrix(thrustweight_time,filename,'Sheet',2,'range','H3:H77')

    %prints out the data for the pressure over time graph
    writematrix(time,filename,'Sheet',2,'range','J3:J77')
    writematrix(pressure_time,filename,'Sheet',2,'range','K3:K77')

end