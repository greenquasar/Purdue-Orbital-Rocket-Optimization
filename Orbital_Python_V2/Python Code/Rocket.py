# Rocket Class
# Holds array of stages and other information about the rocket
#
# Inputs:
# 1. Number of stages
# 2. Array of Isp
# 3. Array of Delta V
# 4. Payload Mass
# 5. Array of propellant mass fraction
#
# Output
# 1. GLOW
#
# **IMPORTANT** All arrays are organized from [Stage 2, Stage 1]

import Stages


class Rocket:
    # Constructor
    def __init__(self, num_Stages, Isp, DeltaVsplit, PayloadMass, PropMassFract):
        self.num_Stages = num_Stages
        self.Isp = Isp
        self.DeltaVsplit = DeltaVsplit
        self.PropMassFract = []
        self.PropMassFract[:] = [x / 100 for x in PropMassFract]
        self.PayloadMass = PayloadMass
        self.stage = []
        self.GLOW = None

    # Add stage objects to list of stages
    def createStages(self):
        for i in range(0, self.num_Stages):
            self.stage.append(Stages.Stages(self.Isp[i], self.PropMassFract[i], None, self.DeltaVsplit[i]))

    # Calculates the mass values for each stage
    def calcStageMass(self):
        for i in range(0, self.num_Stages):
            if i == 0:
                self.stage[i].PayloadMass = self.PayloadMass

            self.stage[i].StageFunction()
            if i < self.num_Stages - 1:
                if self.stage[i].TotalMass < 0:
                    self.GLOW = "Not Possible"
                    break
                else:
                    self.stage[i + 1].PayloadMass = self.stage[i].TotalMass
            else:
                if self.stage[i].TotalMass < 0:
                    self.GLOW = "Not Possible"
                else:
                    self.GLOW = self.stage[i].TotalMass
