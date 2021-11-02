function out = solid_sizing(minTWR, maxTWR, maxPres, deltaV, diameterMin, diameterMax, lengthMin, lengthMax, webMin, webMax, OF)
    addpath('../../SRM Sizing');
    dx = 0.1; %diameter step, m
    dy = 0.1; %length step, m
    dw = 0.1; %web step, m
    shape = 'circular';
    fuel = 'HTPB';
    oxidizer = 'NH4CLO4(I)';
    for diameter = diameterMin:dx:diameterMax
        disp(diameter);
        for length = lengthMin:dy:lengthMax
            disp(length);
            for web = webMin:dw:webMax
                disp(web);
                if web < diameter
                    [T, W, P_c, Thrust, R_b, burn_time, M, Mdot, A_t, deltaV, specificImpulse, propMass, C_t, C_star] = Simulate_Reverse(shape, length, diameter, web, maxPres, OF, fuel, oxidizer);
                end
            end
        end
    end
    out = 1;  %todo
    
end