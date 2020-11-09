from pymeeus.Epoch import Epoch
from pymeeus.JupiterMoons import JupiterMoons

epoch = Epoch()
# epoch.set(1992, 12, 16, 0 + Epoch.tt2ut(1992, 12) / 3600)
epoch.set(1992, 12, 16, 0 + 59 / 3600)
# epoch.set(1992, 10, 13, 0)

print(epoch.jde())
print(JupiterMoons.check_phenomena(epoch))
