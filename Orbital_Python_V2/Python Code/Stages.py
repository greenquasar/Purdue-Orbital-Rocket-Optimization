# Stages Class
# Calculates mass of an individual stage
#
# Inputs:
# 1. Isp
# 2. Propellant Mass Fraction
# 3. Payload Mass
# 4. Velocity
#
# Calculates:
# 1. Propellant mass
# 2. Inert Mass
# 3. Mass Ratio
# 4. Total Mass (Includes mass of payload)

GRAV = 9.81

import math as m

class Stages:
    # Constructor
    def __init__(self, Isp, PropellantMassFract, PayloadMass, VelocitySplit):
        self.Isp = Isp
        self.PropellantMassFract = PropellantMassFract
        self.PayloadMass = PayloadMass
        self.VelocitySplit = VelocitySplit
        self.PropellantMass = None
        self.InertMass = None
        self.MassRatio = None
        self.TotalMass = None

    def MassRatioFun(self):
        self.MassRatio = m.exp((self.VelocitySplit / (GRAV * self.Isp)))
    def PropellantMassFun(self):
        self.PropellantMass = self.PayloadMass * ((self.MassRatio - 1) / (((1 - self.MassRatio) / self.PropellantMassFract) + self.MassRatio))
    def InertMassFun(self):
        self.InertMass = (self.PropellantMass / self.PropellantMassFract) - self.PropellantMass
    def TotalMassFun(self):
        self.TotalMass = self.PropellantMass + self.InertMass + self.PayloadMass

    def StageFunction(self):
        self.MassRatioFun()
        self.PropellantMassFun()
        self.InertMassFun()
        self.TotalMassFun()
