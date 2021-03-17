function [out] = liquid_rocket_engine_optimization(Thrust, Fuel, Fuel_Temp, Oxidizer, Ox_Temp, Chamber_Pressure, OF1, OF2, minChamberDiameter, maxChamberDiameter, stepSize)
    engines = [];
    for Chamber_Diameter = minChamberDiameter:stepSize:maxChamberDiameter
        fclose all;
        fprintf(string(Chamber_Diameter)+'\n');
        engines = [engines; liquid_rocket_engine_design(Thrust, Fuel, Fuel_Temp, Oxidizer, Ox_Temp, Chamber_Pressure, Chamber_Diameter, OF1, OF2)];
    end
    isps = engines(:,3);
    bestIndex = find(max(isps));
    out = engines(bestIndex);
end