function [] = WritetoExcel(DeltaV, Mass, ChbrDiameter, ChbrLength, InnerRad)
    % This function outputs the parameters to an excel file named 'InputOutput.xlsx'
    % 
    % This is poorly written and needs to be fixed but it should work
    % MAKE SURE YOU HAVE EXCEL ON YOUR SYSTEM
    % throws an error if you have Output.xlsx open on your system
    %
    % Make sure that any vector input is a VERTICAL vector, otherwise it will
    % not print properly.

    filename = 'InputOutput.xlsx';

    %prints out the values in the scalar vector above
    writematrix(DeltaV,filename,'Sheet',1,'range','B2:B31')
    writematrix(Mass,filename,'Sheet',1,'range','D2:D31')
    writematrix(ChbrDiameter,filename,'Sheet',1,'range','F2:F31')
    writematrix(ChbrLength,filename,'Sheet',1,'range','H2:H31')
    writematrix(InnerRad,filename,'Sheet',1,'range','J2:J31')

    %prints out the data for the thrust over time graph
    %writematrix(time,filename,'Sheet',2,'range','A3:A77')
    %writematrix(thrust_time,filename,'Sheet',2,'range','B3:B77')

    %prints out the data for the mass over time graph
    %writematrix(time,filename,'Sheet',2,'range','D3:D77')
    %writematrix(mass_time,filename,'Sheet',2,'range','E3:E77')

    %prints out the data for the thrust to weight ratio over time graph
    %writematrix(time,filename,'Sheet',2,'range','G3:G77')
    %writematrix(thrustweight_time,filename,'Sheet',2,'range','H3:H77')

    %prints out the data for the pressure over time graph
    %writematrix(time,filename,'Sheet',2,'range','J3:J77')
    %writematrix(pressure_time,filename,'Sheet',2,'range','K3:K77')

end