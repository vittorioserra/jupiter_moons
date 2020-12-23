from pymeeus.Epoch import Epoch
from pymeeus.Angle import Angle
from pymeeus.Jupiter import Jupiter
from pymeeus.Earth import Earth
import pymeeus.Coordinates
from pymeeus.Sun import Sun
import numpy as np

year = 2020
month = 11
day = 28
epoch = Epoch(year, month, day)
"""        
longitude(Angle) – Geodetic longitude, as an Angle object.It is measured positively west from Greenwich, and negatively
to the east.
latitude(Angle) – Geodetic latitude, as an Angle object
"""
# pos of Friedrichshafen
# NOTE: according to Meeus Chapter 15, Page 101, longitudes eastward from greenwhich should have a negative sign
L = Angle(-9.47554)
phi = Angle(47.65689)

#preliminary calculation to DeltaT
t = (year-2000)/100
#after the year 2000, relevant polynomial
DeltaT = 102+102*t+25.3*(t**2)
#add correction to avoid discontinuity
DeltaT += 0.37*(year-2100)

#for stars and planets
h_0_deg = -0.34

true_obliquity: Angle = pymeeus.Coordinates.true_obliquity(epoch)
nutation_longitude = pymeeus.Coordinates.nutation_longitude(epoch)
# Apparent sidereal time
theta: float = epoch.apparent_sidereal_time(true_obliquity, nutation_longitude)  # in days
theta0 = theta * 360.985647  # in degrees
#theta0 = Angle(theta, ra=True)  # as an angle

alpha1, delta1, elongation = Jupiter.geocentric_position(epoch.__isub__(1))
alpha2, delta2, elongation = Jupiter.geocentric_position(epoch)
alpha3, delta3, elongation = Jupiter.geocentric_position(epoch.__iadd__(1))

#approximate time

cos_H0= (np.sin(np.deg2rad(h_0_deg))-np.sin(np.deg2rad(phi))*np.sin(np.deg2rad(delta2)))/(np.cos(np.deg2rad(phi))*np.cos(np.deg2rad(delta2)))
H0_rad =np.arccos(cos_H0)
H0_deg = np.rad2deg(H0_rad)
#the body is not circumpolar, itcannot remain the whole day abovethe horizon

#times as fractions of a day
m_0 = (alpha2+L-theta0)/360
m_1 = m_0-H0_deg/360
m_2 = m_0 + H0_deg/360

theta_greenwich_0= theta0 + 360.985647*m_0
theta_greenwich_1= theta0 + 360.985647*m_1
theta_greenwich_2= theta0 + 360.985647*m_2