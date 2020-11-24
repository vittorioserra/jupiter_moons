# -*- coding: utf-8 -*-

from pymeeus.Epoch import Epoch
from ..pymeeus.JupiterMoons import JupiterMoons

epoch = Epoch()
epoch.set(2020, 1, 2, 12, 0, 0, utc = True) # this should be in terrestrial time, not utc


# 1 s in jd = 1.157401129603386e-05
time_step = 1.157401129603386e-05 #time-step -->60 seconds

phenomena_matrix = JupiterMoons.check_phenomena(epoch)
print(phenomena_matrix)

print(JupiterMoons.rectangular_positions(epoch))

JupiterMoons.crossing_point(epoch, 3600)

for a in range(1, 3600):

    epoch += time_step
    phenomena_matrix = JupiterMoons.rectangular_positions(epoch, solar=True)

    for i in range(0,1):
        if ((phenomena_matrix[i][0] < 1) and (phenomena_matrix[i][0] > 0)):
            print(epoch.get_full_date())
            print ("Phenomena of satellite " +str(i))
            print(phenomena_matrix[i])

        if ((phenomena_matrix[i][1] < 1) and (phenomena_matrix[i][1] > 0)):
            print(epoch.get_full_date())
            print("Phenomena of satellite " + str(i))
            print(phenomena_matrix[i])

        if ((phenomena_matrix[i][0] > -1) and (phenomena_matrix[i][0] < 0)):
            print(epoch.get_full_date())
            print("Phenomena of satellite " + str(i))
            print(phenomena_matrix[i])

        if ((phenomena_matrix[i][1] > -1) and (phenomena_matrix[i][1] < 0)):
            print(epoch.get_full_date())
            print("Phenomena of satellite " + str(i))
            print(phenomena_matrix[i])


