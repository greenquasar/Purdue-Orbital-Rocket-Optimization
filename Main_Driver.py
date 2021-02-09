# Main Driver for Optimization Software

import GLOW_Optimization_Driver
import xlsxwriter
#import tqdm

# Initializes number of stages
Num_Stages = 2

# Velocity split range and the step size
VelocitySplitRange = [25, 75]
VelocitySplitStep = 1

# Isp range and step size
# Applied to all stages
IspRange = [150, 350]
IspStep = 25

# Propellant mass fraction range and step size
# Applied to all stages
PropMassFractRange = [25, 90]
PropMassFractStep = 5

# Orbit Altitude
Altitude = 150000

# Payload Mass
PayloadMass = 1.33

print("Simulation Started")
Data = GLOW_Optimization_Driver(Num_Stages, VelocitySplitRange, VelocitySplitStep,
                                IspRange, IspStep, PropMassFractRange, PropMassFractStep,
                                Altitude, PayloadMass)

# For printing data from simulation to excel sheet
workbook = xlsxwriter.Workbook('Data_Stages_' + str(Num_Stages) + '.xlsx')
worksheet = workbook.add_worksheet("Sheet 1")

print("Saving Data")

# Excel sheet column titles
title = ["Velocity Split of Stage 1",
            "Velocity Split of Stage 2",
            "Isp of Stage 1",
            "Isp of Stage 2",
            "Propellant Mass Fraction of Stage 1",
            "Propellant Mass Fraction of Stage 2",
            "GLOW"]
worksheet.write_row(0, 0, title)

# prints data to excel sheet
sheet = 1
col = 0
for row in tqdm(range(len(Data))):
    worksheet.write_row((row+sheet) % 1048576, col, Data[row])

    if (((row+sheet+1) % 1048576) == 0) & (row > 2):
        sheet = sheet + 1
        worksheet = workbook.add_worksheet("Sheet " + str(sheet))
        worksheet.write_row(0, 0, title)
workbook.close()
print("Simulation complete")
