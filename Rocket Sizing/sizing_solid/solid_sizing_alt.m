function out = solid_sizing_alt(TWRtarget, TWRtol, deltaVtarget, deltaVtol, maxPres)
    addpath('../../SRM Sizing');
    %define what we know
    shape = 'circular';
    fuel = 'HTPB';
    oxidizer = 'NH4CLO4(I)';
    OF = 1;
    %todo make these inputs ^^
    
    %guess at everything else
    length = 1;
    diameter = 1;
    web = 0.25 * diameter;
    
    while true %do while
        fprintf("Simulating with length=%.2fm, diameter=%.2fm, web=%.2fm\n", length, diameter, web);
        [T, W, P_c, Thrust, R_b, burn_time, M, Mdot, A_t, deltaV, specificImpulse, propMass, C_t, C_star] = Simulate_Reverse(shape, length, diameter, web, maxPres, OF, fuel, oxidizer);
        deltaVdiff = deltaV-deltaVtarget;
        deltaVpercentdiff = deltaVdiff/deltaVtarget;
        TWRdiff = (T/W)-TWRtarget;
        TWRpercentdiff = TWRdiff/TWRtarget;
        fprintf("Current rocket has ΔV=%.0fm/s, TWR=%0.2f\n", deltaV, T/W);
        fprintf("Target rocket has  ΔV=%.0fm/s, TWR=%0.2f\n", deltaVtarget, TWRtarget);
        
        if ((abs(deltaVpercentdiff) < deltaVtol) && (abs(TWRpercentdiff) < TWRtol)) %condition
            break; 
        if ()%dv to low or high
          %todo change l, d, w
        if () %twr too low or high
          %todo change l, d, w
      end
    end
    out = 1;
end