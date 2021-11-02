function out = solid_sizing(TWRtarget, TWRtol, deltaVtarget, deltaVtol, maxPres, diameterMin, diameterMax, OF)
    addpath('../../SRM Sizing');
    %define what we know
    shape = 'circular';
    fuel = 'HTPB';
    oxidizer = 'NH4CLO4(I)';
    
    %guess at everything else
    
    
    
    while true %do while
        [T, W, P_c, Thrust, R_b, burn_time, M, Mdot, A_t, deltaV, specificImpulse, propMass, C_t, C_star] = Simulate_Reverse(shape, length, diameter, web, maxPres, OF, fuel, oxidizer);
        if ((abs(deltaV-deltaVtarget)/deltaVtarget < deltaVtol) && (abs((T/W)-TWRtarget)/TWRtarget < TWRtol)) %condition
            break; 
        if ()%dv to low or high
          %todo
        if () %twr too low or high
          %todo
      end
    end
    out = 1;
end