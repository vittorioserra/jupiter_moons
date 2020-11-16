# -*- coding: utf-8 -*-

from pymeeus.Epoch import Epoch
from pymeeus.JupiterMoons import JupiterMoons

epoch = Epoch()
epoch.set(1992, 12, 19, 5 + 59 / 3600)

# 1 s in jd = 1.157401129603386e-05
time_step = 60 * 5 * 1.157401129603386e-05

for i in range(100) :

    print("Alpha, Central Point, Beta, Spread" + str(JupiterMoons.round_base_cone_param(epoch + i*time_step)))

epoch.set(1997, 12, 19, 5 + 59 / 3600)

print("On the other half of the orbit, the results should be almost the same")

for i in range(100) :

    print("Alpha, Central Point, Beta, Spread" + str(JupiterMoons.round_base_cone_param(epoch + i*time_step)))

epoch.set(2003, 12, 19, 5 + 59 / 3600)

print("On the other half of the orbit, the results should be almost the same")

for i in range(100) :

    print("Alpha, Central Point, Beta, Spread" + str(JupiterMoons.round_base_cone_param(epoch + i*time_step)))

epoch.set(2008, 12, 19, 5 + 59 / 3600)

print("On the other half of the orbit, the results should be almost the same")

for i in range(100) :

    print("Alpha, Central Point, Beta, Spread" + str(JupiterMoons.round_base_cone_param(epoch + i*time_step)))


epoch = Epoch()
epoch.set(1992, 12, 19, 5 + 59 / 3600)

# 1 s in jd = 1.157401129603386e-05
time_step = 60 * 5 * 1.157401129603386e-05

for i in range(0, 100):

    epoch += time_step
    print(epoch.get_full_date())
    jupiter_moons_tuple = JupiterMoons.rectangular_positions(epoch, solar = True)
    moon_tuple = jupiter_moons_tuple[0][:]
    print(str(moon_tuple))
    cone = JupiterMoons.round_base_cone_param(epoch)
    print("Alpha, Central Point, Beta, Spread" + str(JupiterMoons.round_base_cone_param(epoch + i*time_step)))
    cone_size_moon = JupiterMoons.check_umbra_eclipse(moon_tuple[0], moon_tuple[1], moon_tuple[2], cone[0])
    print("Specific cone size :")
    print(str(cone_size_moon))
