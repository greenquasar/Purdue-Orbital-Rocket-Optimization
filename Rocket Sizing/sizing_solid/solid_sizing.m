function [mass] = solid_sizing(max_thrust, chamb_pressure, deltaV, chamb_dia_min, chamb_dia_max, chamb_len_min, chamb_len_max, fuel_ox)
    chamb_dia = chamb_dia_min;
    while chamb_dia <= chamb_dia_max
        %call function
        chamb_dia = chamb_dia + .1;
    end
    