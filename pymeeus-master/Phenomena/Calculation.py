# -*- coding: utf-8 -*-

from pymeeus.Epoch import Epoch
from pymeeus.JupiterMoons import JupiterMoons

epoch = Epoch()
epoch.set(1992, 12, 16, 0 + 59 / 3600)

time_step = 3600 * 1.157401129603386e-05

for i in range(1, 100):
    epoch += time_step
    print(epoch.get_full_date())
    print(JupiterMoons.check_phenomena(epoch))
