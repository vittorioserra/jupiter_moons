from pymeeus.Epoch import Epoch
from pymeeus.Angle import Angle
from pymeeus.Jupiter import Jupiter
from pymeeus.Earth import Earth
import pymeeus.Coordinates
from pymeeus.Sun import Sun


def rising_setting_jupiter(epoch: Epoch, longitude, latitude):
    # Time of interest
    year, month, day = epoch.get_date()
    """
    alpha1(Angle) – Apparent right ascension the previous day at 0h TT, as an Angle object
    delta1(Angle) – Apparent declination the previous day at 0 h TT, as an Angle object

    alpha2(Angle) – Apparent right ascension the current day at 0 h TT, as an Angle object
    delta2(Angle) – Apparent declination the current day at 0 h TT, as an Angle object

    alpha3(Angle) – Apparent right ascension the following day at 0h TT, as an Angle object
    delta3(Angle) – Apparent declination the following day at 0 h TT, as an Angle object
    """
    alpha1, delta1, elongation = Jupiter.geocentric_position(epoch.__isub__(1))
    alpha2, delta2, elongation = Jupiter.geocentric_position(epoch)
    alpha3, delta3, elongation = Jupiter.geocentric_position(epoch.__iadd__(1))

    """
    h0(Angle) – ‘Standard’ altitude: the geometric altitude of the center of the body at the time of apparent rising
    or setting, as degrees in an Angle object.It should be - 0.5667 deg for stars and planets, -0.8333 deg for the
    Sun, and 0.125 deg for the Moon.
    """
    h0 = Angle(-0.5667)
    """
    delta_t(float) – The difference between Terrestrial Time and Universal Time(TT - UT) in seconds of time
    """
    delta_t: float = Epoch.tt2ut(year, month)
    """
    theta0(Angle) – Apparent sidereal time at 0 h TT on the current day for the meridian of Greenwich, as degrees 
    in an Angle object
    """
    true_obliquity: Angle = pymeeus.Coordinates.true_obliquity(epoch)
    nutation_longitude = pymeeus.Coordinates.nutation_longitude(epoch)
    # Apparent sidereal time
    theta: float = epoch.apparent_sidereal_time(true_obliquity, nutation_longitude)  # in days
    theta = theta / 24.0  # in hours
    theta0 = Angle(theta, ra=True)  # as an angle

    rising, transit, setting = pymeeus.Coordinates.times_rise_transit_set(longitude, latitude, alpha1, delta1,
                                                                          alpha2, delta2, alpha3, delta3, h0,
                                                                          delta_t, theta0)
    epoch_rising = Epoch(year, month, day, rising)
    epoch_transit = Epoch(year, month, day, transit)
    epoch_setting = Epoch(year, month, day, setting)

    y_rising, m_rising, d_rising, h_rising, mi_rising, s_rising = epoch_rising.get_full_date()
    y_transit, m_transit, d_transit, h_transit, mi_transit, s_transit = epoch_transit.get_full_date()
    y_setting, m_setting, d_setting, h_setting, mi_setting, s_setting = epoch_setting.get_full_date()

    if h_transit < h_rising:
        d_transit += 1
    if h_setting < h_rising:
        d_setting += 1

    print("Rising time:  " + str(y_rising) + "/" + str(m_rising) + "/" + str(d_rising) + ", " +
          str(h_rising) + ":" + str(mi_rising) + ":" + str(int(s_rising)))

    print("Transit time: " + str(y_transit) + "/" + str(m_transit) + "/" + str(d_transit) + ", " +
          str(h_transit) + ":" + str(mi_transit) + ":" + str(int(s_transit)))

    print("Setting time: " + str(y_setting) + "/" + str(m_setting) + "/" + str(d_setting) + ", " +
          str(h_setting) + ":" + str(mi_setting) + ":" + str(int(s_setting)))

    print("\nSoll (für 28.11.2020):\n"
          "\tRising = 11:19\n\tTransit = 15:40\n\tSetting = 20:02")


def rising_setting_sun(epoch, longitude, latitude):
    # Time of interest
    year, month, day = epoch.get_date()
    """
    alpha1(Angle) – Apparent right ascension the previous day at 0h TT, as an Angle object
    delta1(Angle) – Apparent declination the previous day at 0 h TT, as an Angle object

    alpha2(Angle) – Apparent right ascension the current day at 0 h TT, as an Angle object
    delta2(Angle) – Apparent declination the current day at 0 h TT, as an Angle object

    alpha3(Angle) – Apparent right ascension the following day at 0h TT, as an Angle object
    delta3(Angle) – Apparent declination the following day at 0 h TT, as an Angle object
    """
    alpha1, delta1, r = Sun.apparent_rightascension_declination_coarse(epoch.__isub__(1))
    alpha2, delta2, r = Sun.apparent_rightascension_declination_coarse(epoch)
    alpha3, delta3, r = Sun.apparent_rightascension_declination_coarse(epoch.__iadd__(1))

    """
    h0(Angle) – ‘Standard’ altitude: the geometric altitude of the center of the body at the time of apparent rising
    or setting, as degrees in an Angle object.It should be - 0.5667 deg for stars and planets, -0.8333 deg for the
    Sun, and 0.125 deg for the Moon.
    """
    h0 = Angle(-0.8333)
    """
    delta_t(float) – The difference between Terrestrial Time and Universal Time(TT - UT) in seconds of time
    """
    delta_t: float = Epoch.tt2ut(year, month)
    """
    theta0(Angle) – Apparent sidereal time at 0 h TT on the current day for the meridian of Greenwich, as degrees 
    in an Angle object
    """
    true_obliquity: Angle = pymeeus.Coordinates.true_obliquity(epoch)
    nutation_longitude = pymeeus.Coordinates.nutation_longitude(epoch)
    # Apparent sidereal time
    theta: float = epoch.apparent_sidereal_time(true_obliquity, nutation_longitude)  # in days
    theta = theta / 24.0  # in hours
    theta0 = Angle(theta, ra=True)  # as an angle

    rising, transit, setting = pymeeus.Coordinates.times_rise_transit_set(longitude, latitude, alpha1, delta1,
                                                                          alpha2, delta2, alpha3, delta3, h0,
                                                                          delta_t, theta0)
    epoch_rising = Epoch(year, month, day, rising)
    epoch_transit = Epoch(year, month, day, transit)
    epoch_setting = Epoch(year, month, day, setting)

    y_rising, m_rising, d_rising, h_rising, mi_rising, s_rising = epoch_rising.get_full_date()
    y_transit, m_transit, d_transit, h_transit, mi_transit, s_transit = epoch_transit.get_full_date()
    y_setting, m_setting, d_setting, h_setting, mi_setting, s_setting = epoch_setting.get_full_date()

    if h_transit < h_rising:
        d_transit += 1
    if h_setting < h_rising:
        d_setting += 1

    print("Rising time:  " + str(y_rising) + "/" + str(m_rising) + "/" + str(d_rising) + ", " +
          str(h_rising) + ":" + str(mi_rising) + ":" + str(int(s_rising)))

    print("Transit time: " + str(y_transit) + "/" + str(m_transit) + "/" + str(d_transit) + ", " +
          str(h_transit) + ":" + str(mi_transit) + ":" + str(int(s_transit)))

    print("Setting time: " + str(y_setting) + "/" + str(m_setting) + "/" + str(d_setting) + ", " +
          str(h_setting) + ":" + str(mi_setting) + ":" + str(int(s_setting)))


class VisibilityJupiter:
    def __init__(self, epoch, longitude, latitude):
        print("\nJUPITER")
        rising_setting_jupiter(epoch, longitude, latitude)
        print("\nSUN")
        rising_setting_sun(epoch, longitude, latitude)


if __name__ == "__main__":
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
    longitude = Angle(9.47554)
    latitude = Angle(47.65689)

    VisibilityJupiter(epoch, longitude, latitude)
