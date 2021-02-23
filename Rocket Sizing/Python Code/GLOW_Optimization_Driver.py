# Optimization Function
# Iterates through Velocity split range
#
# Inputs:
# 1. Number of Stages
# 2. Array of velocity split range
# 3. Velocity split step size
# 4. Array of Isp Range
# 5. Isp step size
# 6. Array of propellant mass fraction range
# 7. Propellant mass fraction step size
# 5. Altitude
# 6. PayloadMass
#
# Output
# 1. GLOW for each configuration
# 2. Array of all tabulated data

import Rocket
import numpy as np
from tqdm import tqdm


def GLOW_Optimization_Driver(Num_Stages, VelocitySplitRange, VelocitySplitStep,
                                IspRange, IspStep, PropMassFractRange, PropMassFractStep,
                                Altitude, PayloadMass):

    # Array to hold data
    Data = []

    # Determines Ideal Velocity
    Videal = 1.2 * np.sqrt(3.9857E14/(Altitude + 6378000))

    # Loop through range of Velocity Split
    # Split is the velocity percentage of the second stage
    for split in tqdm(range(VelocitySplitRange[0], VelocitySplitRange[1] + 1, VelocitySplitStep)):
        DeltaVsplit = [Videal*split/100, Videal*(100 - split)/100]

        # Loops through Isp of first stage
        for Isp1 in range(IspRange[0], IspRange[1] + 1, IspStep):

            # Loops through Isp of second stage
            for Isp2 in range(IspRange[0], IspRange[1] + 1, IspStep):

                # Loops through Propellant Mass Fraction of First Stage
                for PropMassFract1 in range(PropMassFractRange[0], PropMassFractRange[1] + 1, PropMassFractStep):

                    # Loops through propellant Mass Fraction of Second Stage
                    for PropMassFract2 in range(PropMassFractRange[0], PropMassFractRange[1] + 1, PropMassFractStep):
                        Isp = [Isp2, Isp1]
                        PropMassFract = [PropMassFract2, PropMassFract1]
                        rocket = Rocket.Rocket(Num_Stages, Isp, DeltaVsplit, PayloadMass, PropMassFract)
                        rocket.createStages()
                        rocket.calcStageMass()
                        Data.append([100 - split, split, Isp1, Isp2, PropMassFract1, PropMassFract2, rocket.GLOW])
    return Data
