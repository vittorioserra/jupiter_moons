import functools
from enum import unique, IntEnum
from typing import List

from pymeeus.Epoch import Epoch
from pymeeus.Jupiter import Jupiter
from pymeeus.JupiterMoons import JupiterMoons
from pymeeus.Earth import Earth
from pymeeus.Coordinates import mean_obliquity
import numpy as np
from numpy import sin, cos, sqrt, deg2rad, tan, arctan2, arcsin


from math import atan, tan, sqrt

# Radius of Jupiter moons in Jupiter's radii
# moons_radii = [3643.2 / (2 * 69911), 3121.6 / (2 * 69911), 5262.4 / (2 * 69911), 4820.6 / (2 * 69911)]
moons_radii = [3643.2 / (2 * 71492), 3121.6 / (2 * 71492), 5262.4 / (2 * 71492), 4820.6 / (2 * 71492)]


@functools.total_ordering
class Result:
    def __init__(self, epoch: Epoch):
        """ Constructor for Results class

        :param epoch: Epoch the Result was calculated for
        :type epoch: Epoch
        """

        # Initiate variables
        self.epoch = epoch

    def __lt__(self, other):
        """ Compares Result to an other Result depending
        on the Epoch the calculation was made for

        :param other: Other Result
        :type other: Any

        :returns: Whether other Result has higher Epoch or not
        :rtype: bool

        :raises: TypeError, if other is no instance of Result
        """

        # Check type
        if not isinstance(other, Result):
            raise TypeError("Can't compare different types")

        # Compare Epochs
        return self.epoch.jde() < other.epoch.jde()

    def __eq__(self, other):
        """ Checks if Result is equal to an other Result depending
        on the Epoch the calculation was made for

        :param other: Other Result
        :type other: Any

        :returns: Whether other Result has same Epoch or not
        :rtype: bool

        :raises: TypeError, if other is no instance of Result
        """

        # Check type
        if not isinstance(other, Result):
            raise TypeError("Can't compare different types")

        # Compare Epochs
        return self.epoch.jde() == other.epoch.jde()


class Phenomenon(Result):
    """Class that delivers calculation data for phenomena
    """

    def __init__(self, phenomenon_type: str, sat: int, distance: float, z: float, epoch: Epoch,
                 shadow_radius_umbra: float = 0,
                 shadow_radius_penumbra: float = 0):
        """Constructor for Phenomenon class.
        Sortable by Epoch

        :param phenomenon_type: Type of phenomenon:
            PA: Transit
            OC: Ocultation
            OM: Transit of the umbra of the satellite on the disk of Jupiter
            EC: Eclipse
        :type phenomenon_type: str
        :param sat: Number of satellite - 1
            0: Io
            1: Europa
            2: Ganymede
            3: Callisto
        :type sat: int
        :param distance: Perspective distance to center of Jupiter
        :type distance: float
        :param z: z-Coordinate of satellite
        :type z: float
        :param epoch: Epoch of the Phenomena
        :type epoch: Epoch
        :param shadow_radius_umbra: Radius of the umbra in Jupiter's radii
            (only solar phenomena)
        :type shadow_radius_umbra: float
        :param shadow_radius_penumbra: Radius of the penumbra in Jupiter's radii
            (only solar phenomena)
        :type shadow_radius_penumbra: float
        """

        # Call Result constructor
        super().__init__(epoch)

        # Initialize variables
        self.phenomenon_type = phenomenon_type
        self.sat = sat
        self.distance = distance
        self.z = z
        self.shadow_radius_umbra = shadow_radius_umbra
        self.shadow_radius_penumbra = shadow_radius_penumbra

        self.moon_radius = moons_radii[sat]

        self.phenomenon_occurs = False
        self.shadow_type: str = ""

        # Check shadow phenomenon type
        self._check_shadow_type()

    def get_critical_radius(self) -> float:
        """Method that calculates the perspective radius (distance) to Jupiter
        that is at least necessary to have an extreme shadow phenomenon.

        :returns: Perspective distance (in Jupiter's radii) to Jupiter for
            extreme shadow phenomena
        :rtype: float
        """

        # For eclipse return penumbra radius + moon radius
        if self.phenomenon_type == "EC":
            return self.shadow_radius_penumbra + self.moon_radius
        # For other phenomena return Jupiter's radius (= 1) + moon radius
        else:
            return 1 + self.moon_radius

    def has_right_z_sign(self) -> bool:
        """Determines whether the sign of the z-Coordinate corresponds with the phenomenon
        type. OC and EC only possible for more distant moons, OM and PA for moons that are
        less distant than Jupiter.

        :returns: Whether the sign of the z-Coordinate corresponds with the
            phenomenon type
        :rtype: bool
        """

        return (self.phenomenon_type in ["PA", "OM"] and self.z < 0) or (
                self.phenomenon_type in ["OC", "EC"] and self.z > 0)

    def update(self, epoch: Epoch) -> None:
        """Updates the Phenomenon object to a given Epoch

        :param epoch: Epoch the phenomenon shall be updated
        :type epoch: Epoch

        :rtype: None
        """

        # Calculate updated phenomenon
        updated_phenomenon = Detection.check_single_phenomena(epoch, self.sat, self.phenomenon_type)

        # Copy all class variables from updated phenomenon
        self.phenomenon_type = updated_phenomenon.phenomenon_type
        self.sat = updated_phenomenon.sat
        self.distance = updated_phenomenon.distance
        self.z = updated_phenomenon.z
        self.epoch = updated_phenomenon.epoch
        self.shadow_radius_umbra = updated_phenomenon.shadow_radius_umbra
        self.shadow_radius_penumbra = updated_phenomenon.shadow_radius_penumbra
        self.moon_radius = updated_phenomenon.moon_radius
        self.phenomenon_occurs = updated_phenomenon.phenomenon_occurs
        self.shadow_type = updated_phenomenon.shadow_type

    def _check_shadow_type(self) -> None:
        """Private method that determines shadow type and whether a
        phenomenon occurs

        :rtype: None
        """

        # Check z-Coordinate if Phenomenon is possible
        if self.has_right_z_sign():
            # Set flag
            self.phenomenon_occurs = True

            # Phenomena with EXT and INT
            if self.phenomenon_type in ["PA", "OC", "OM"]:
                # Exterior contact with shadow
                if 1 - self.moon_radius < self.distance <= 1 + self.moon_radius:
                    self.shadow_type = "EXT"
                # Interior contact with shadow
                elif self.distance <= 1 - self.moon_radius:
                    self.shadow_type = "INT"
                else:
                    # Unset flag, no phenomenon detected
                    self.phenomenon_occurs = False

            # Phenomena with PEN
            else:
                # Exterior contact with penumbra and closer
                if self.distance <= self.shadow_radius_penumbra + self.moon_radius:
                    self.shadow_type = "PEN"

                    if self.shadow_radius_umbra - self.moon_radius < self.distance <= self.shadow_radius_umbra + self.moon_radius:
                        # Exterior contact with shadow
                        self.shadow_type = "EXT"
                    elif self.distance <= self.shadow_radius_umbra - self.moon_radius:
                        # Interior contact with shadow
                        self.shadow_type = "INT"
                else:
                    # Unset flag, no phenomenon detected
                    self.phenomenon_occurs = False


class Detection:
    """This class delivers methods for checking for any kind of phenomena
    of Jupiter's four Galilean moons.
    """

    @unique
    class PhenomenaTypes(IntEnum):
        OC = 0
        EC = 1
        PA = 2
        OM = 3

    # Define phenomena types as strings
    phenomena_types_str = ["OC", "EC", "PA", "OM"]

    @staticmethod
    def check_single_phenomena(epoch: Epoch, i_sat: int, phenomenon_type: str) -> Phenomenon:
        """Calculates a Phenomenon object for a given satellite and phenomenon
        type to a given Epoch

        :param epoch: Epoch that shall be checked
        :type epoch: Epoch
        :param i_sat: Index of satellite (i.e. No. of Sat. - 1) that shall
            be checked
        :type i_sat: int
        :param phenomenon_type: Type of the phenomenon being checked
        :type phenomenon_type: str

        :returns: Phenomenon object for the given inputs
        :rtype: Phenomenon

        :raises: TypeError, if no Epoch typed input is given
        """

        # Check input type
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")

        # Calc perspective coordinates of the satellites depending on phenomenon type
        Coords = JupiterMoons.rectangular_positions_jovian_equatorial(epoch, solar=phenomenon_type in ["EC", "OM"])

        # Earth based phenomena
        if phenomenon_type in ["PA", "OC"]:
            # Calc perspective distance to Jupiter for given satellite
            distance = Detection.perspective_distance(Coords[i_sat][0], Coords[i_sat][1])

            # Return Phenomena instance
            return Phenomenon(Detection.phenomena_types_str[Detection.phenomena_types_str.index(phenomenon_type)],
                              i_sat, distance, Coords[i_sat][2], epoch)

        # Solar phenomena
        else:
            # Calc perspective distance to Jupiter for given satellite
            distance = Detection.perspective_distance(Coords[i_sat][0], Coords[i_sat][1])
            # Calc shadow cone parameter
            r_umbra, r_penumbra = Detection.cone_radius(epoch, Coords[i_sat][2])

            # Return Phenomena instance
            return Phenomenon(Detection.phenomena_types_str[Detection.phenomena_types_str.index(phenomenon_type)],
                              i_sat, distance, Coords[i_sat][2], epoch, r_umbra,
                              r_penumbra)

    @staticmethod
    def check_phenomena(epoch: Epoch) -> List[List[Phenomenon]]:
        """Calculates Phenomenon objects for all satellite and phenomena
        types for a given Epoch

        :param epoch: Epoch that shall be checked
        :type epoch: Epoch

        :returns: 2D-Array of Phenomena object, where rows correspond with
            satellites, columns with phenomena types
        :rtype: List[List[Phenomenon]]

        :raises: TypeError, if no Epoch typed input is given
        """

        # Check input type
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")

        # Result matrix, where each rows is for a satellite
        # Column 0: Ocultation (OC)
        # Column 1: Eclipse (EC)
        # Column 2: Transit (PA)
        # Column 3: Satellite umbra transit (OM)
        # noinspection PyTypeChecker
        result_matrix: List[List[Phenomenon]] = [[None, None, None, None],
                                                 [None, None, None, None],
                                                 [None, None, None, None],
                                                 [None, None, None, None]]

        # Calculate coordinates as seen from the Earth
        Coords_Earth = JupiterMoons.rectangular_positions_jovian_equatorial(epoch)
        # Calculate coordinates as seen from the Sun
        Coords_Sun = JupiterMoons.rectangular_positions_jovian_equatorial(epoch, solar=True)

        # Iterate through satellites
        for i_sat in range(0, 4):
            # Iterate through phenomena types
            for phenomenon_type in Detection.PhenomenaTypes:
                if phenomenon_type in [Detection.PhenomenaTypes.PA, Detection.PhenomenaTypes.OC]:
                    # Calc perspective distance to Jupiter
                    distance = Detection.perspective_distance(Coords_Earth[i_sat][0], Coords_Earth[i_sat][1])

                    # Fill result matrix with Phenomena instance
                    result_matrix[i_sat][phenomenon_type] = Phenomenon(Detection.phenomena_types_str[phenomenon_type],
                                                                       i_sat, distance, Coords_Earth[i_sat][2], epoch)
                else:
                    # Calc perspective distance to Jupiter
                    distance = Detection.perspective_distance(Coords_Sun[i_sat][0], Coords_Sun[i_sat][1])
                    # Calc shadow cone parameter
                    r_umbra, r_penumbra = Detection.cone_radius(epoch, Coords_Sun[i_sat][2])

                    # Fill result matrix with Phenomena instance
                    result_matrix[i_sat][phenomenon_type] = Phenomenon(Detection.phenomena_types_str[phenomenon_type],
                                                                       i_sat, distance, Coords_Sun[i_sat][2], epoch,
                                                                       r_umbra,
                                                                       r_penumbra)

        return result_matrix

    @staticmethod
    def perspective_distance(X: float, Y: float, ellipsoid: bool = True):
        """This method calculates the perspective distance to the center of Jupiter
        for the given coordinates. It returns the distance in Jupiter's radii.

        :param X: X-coordinate of the satellite in Jupiter's radii
        :type X: float
        :param Y: Y-coordinate of the satellite in Jupiter's radii
        :type Y: float
        :param ellipsoid: Whether Jupiter shall be modelled elliptical
        :type ellipsoid: bool

        :returns: Perspective distance to Jupiter's center in Jupiter's radii
        :rtype: float
        """

        if ellipsoid:
            # Distort y-coordinate to model elliptical Jupiter
            Y *= 1.071374

        return sqrt(X ** 2 + Y ** 2)

    @staticmethod
    def cone_radius(epoch: Epoch, z: float):
        """Calculates the radius of umbra and penumbra shadow for
        a given z-Coordinate (z > 0 -> more distant than Jupiter)

        :param epoch: Epoch the calculation should be made for
        :type epoch: Epoch
        :param z: Z-Coordinate in Jupiter's radii
        :type z: float

        :returns: Radius of umbra and penumbra shadow in Jupiter's
            radii
        """

        alpha_rad, cone_alpha_vertex, beta_rad, cone_beta_vertex = Detection.round_base_cone_param(epoch)

        r_umbra = (abs(cone_alpha_vertex) - z) / abs(cone_alpha_vertex)
        r_penumbra = (abs(cone_beta_vertex) + z) / abs(cone_beta_vertex)

        return r_umbra, r_penumbra

    @staticmethod
    def planetocentric_declinations_rad(epoch):
        d = epoch - 2433282.5
        T1 = d.jde() / 36525

        alpha_0 = np.deg2rad(268 + 0.1061 * T1)
        delta_0 = np.deg2rad(64.50 - 0.0164 * T1)

        l, b, r = Jupiter.geometric_heliocentric_position(epoch, tofk5=False)

        l = l.rad()
        b = b.rad()

        # Compute the heliocentric position of the Earth
        l0, b0, r0 = Earth.geometric_heliocentric_position(epoch, tofk5=False)

        l0 = l0.rad()
        b0 = b0.rad()

        # jupiter's data
        x = r * np.cos(b) * np.cos(l) - r * cos(l0)
        y = r * np.cos(b) * np.sin(l) - r * sin(l0)
        z = r * np.sin(b) - r * cos(b0)

        Delta = np.sqrt(x ** 2 + y ** 2 + z ** 2)

        l_corr_factor_deg = -0.012990 * Delta / r ** 2

        l_corr_factor_rad = np.deg2rad(l_corr_factor_deg)

        l = l + l_corr_factor_rad

        # calculate x, y, z and Delta anew

        x = r * np.cos(b) * np.cos(l) - r0 * cos(l0)
        y = r * np.cos(b) * np.sin(l) - r0 * sin(l0)
        z = r * np.sin(b) - r0 * sin(b0)

        Delta = np.sqrt(x ** 2 + y ** 2 + z ** 2)

        epsilon_0 = (mean_obliquity(epoch)).rad()

        alpha_s = np.arctan2((np.cos(epsilon_0) * np.sin(l) - np.sin(epsilon_0) * np.tan(b)), np.cos(l))
        delta_s = np.arcsin(cos(epsilon_0) * sin(b) + sin(epsilon_0) * cos(b) * sin(l))

        D_s = -sin(delta_0) * sin(delta_s) - cos(delta_0) * cos(delta_s) * cos(alpha_0 - alpha_s)

        u = y * cos(epsilon_0) - z * sin(epsilon_0)
        v = y * sin(epsilon_0) + z * cos(epsilon_0)

        alpha = arctan2(u, x)
        delta = arctan2(v, sqrt(x ** 2 + u ** 2))
        chsi = arctan2(sin(delta_0) * cos(delta) * cos(alpha_0 - alpha) - sin(delta) * cos(delta_0),
                       cos(delta) * sin(alpha_0 - alpha))

        D_e = -sin(delta_0) * sin(delta) - cos(delta_0) * cos(delta) * cos(alpha_0 - alpha)

        return D_s, D_e


    @staticmethod
    def ellipse_base_cone_param(epoch: Epoch, ellipsoid: bool = True):
        """This method constructs a cone modelling jupiter's shadow, thus here we are considering the rotated ellipsoidal
        shape of jupiter

        :param epoch: Epoch that should be checked
        :type epoch: Epoch
        :param ellipsoid: Whether or not to distort jupiter and the cone
            concurrently
        :type ellipsoid: bool

        :returns:
            alpha_cone_rad : aperture of the cone in radians measured from
                the base,
            cone_vertex_jupiter_radii : distance of the umbral cone's sharpest
                point in jupiter-radii, always behind jupiter if viewed
                from the sun,
            beta_cone_rad : aperture of the penumbral cone in radians measured
                from jupiter,
            cone_beta_vertex_jupiter_radii : distance of the penumbral cone's
                sharpest point in jupiter-radii, always before jupiter if
                viewed from the sun
        :rtype: tuple

        :raises: TypeError, if input parameter has invalid type
        """

        #compute the angles
        D_s, D_e = planetocentric_declinations_rad(epoch)

        # Define radii of sun and Jupiter in AU
        sun_radius_au = 0.00465047
        jupiter_radius_au = 0.10045 * sun_radius_au

        if ellipsoid:
            jupiter_radius_multiplier = 1.071374
        else:
            jupiter_radius_multiplier = 1.0

        jupiter_radius_au = jupiter_radius_au * 1.071374

        # Check if Epoch is given
        if epoch is not None:
            # Check types
            if isinstance(epoch, Epoch):
                # Calculate the position of jupiter in solar-spherical coordinates
                l, b, r = Jupiter.geometric_heliocentric_position(epoch)

                # alpha is the umbra defining angle
                alpha_cone_rad_equatorial = atan(r / (sun_radius_au - jupiter_radius_au))
                alpha_cone_rad_polar = atan(r / (sun_radius_au - jupiter_radius_au*np.cos(D_s)))

                # beta is the penumbra defying angle
                beta_cone_rad_equatorial = atan(r / (sun_radius_au + jupiter_radius_au))
                beta_cone_rad_polar = atan(r / (sun_radius_au + jupiter_radius_au*np.cos(D_s)))

                # Compute distance of the sharpest pint behind jupiter in jupiter radii
                cone_vertex_jupiter_radii_equatorial = tan(alpha_cone_rad_equatorial)
                cone_vertex_jupiter_radii_polar = tan(alpha_cone_rad_polar)

                cone_beta_vertex_jupiter_radii_equatorial = (-1 * jupiter_radius_multiplier * tan(beta_cone_rad_equatorial))
                cone_beta_vertex_jupiter_radii_polar = (-1 * jupiter_radius_multiplier * tan(beta_cone_rad_polar))

                return [alpha_cone_rad_equatorial,alpha_cone_rad_polar], \
                       [cone_vertex_jupiter_radii_equatorial, cone_vertex_jupiter_radii_polar], \
                       [beta_cone_rad_equatorial, beta_cone_rad_polar],\
                       [cone_beta_vertex_jupiter_radii_equatorial, cone_beta_vertex_jupiter_radii_polar]

            else:
                raise TypeError("Invalid input type")

    @staticmethod
    def elliptical_cone_radius(epoch: Epoch, z: float):
        """Calculates the radius of umbra and penumbra shadow for
        a given z-Coordinate (z > 0 -> more distant than Jupiter)

        :param epoch: Epoch the calculation should be made for
        :type epoch: Epoch
        :param z: Z-Coordinate in Jupiter's radii
        :type z: float

        :returns: Radius of umbra and penumbra shadow in Jupiter's
            radii
        """
        # compute the angles

        alpha_cone_rad,cone_alpha_vertex,beta_cone_rad , cone_beta_vertex = Detection.ellipse_base_cone_param(epoch)

        r_umbra_equatorial = (abs(cone_alpha_vertex[0]) - z) / abs(cone_alpha_vertex[0])
        r_umbra_polar = (abs(cone_alpha_vertex[1]) - z) / abs(cone_alpha_vertex[1])
        r_penumbra_equatorial = (abs(cone_beta_vertex[0]) + z) / abs(cone_beta_vertex[0])
        r_penumbra_polar = (abs(cone_beta_vertex[1]) + z) / abs(cone_beta_vertex[1])

        r_umbra =[r_umbra_equatorial, r_umbra_polar]
        r_penumbra = [r_penumbra_polar, r_penumbra_equatorial]

        return r_umbra, r_penumbra

    @staticmethod
    def round_base_cone_param(epoch: Epoch, ellipsoid: bool = True):
        """This method constructs a cone modelling jupiter's shadow

        :param epoch: Epoch that should be checked
        :type epoch: Epoch
        :param ellipsoid: Whether or not to distort jupiter and the cone
            concurrently
        :type ellipsoid: bool

        :returns:
            alpha_cone_rad : aperture of the cone in radians measured from
                the base,
            cone_vertex_jupiter_radii : distance of the umbral cone's sharpest
                point in jupiter-radii, always behind jupiter if viewed
                from the sun,
            beta_cone_rad : aperture of the penumbral cone in radians measured
                from jupiter,
            cone_beta_vertex_jupiter_radii : distance of the penumbral cone's
                sharpest point in jupiter-radii, always before jupiter if
                viewed from the sun
        :rtype: tuple

        :raises: TypeError, if input parameter has invalid type
        """

        # Define radii of sun and Jupiter in AU
        sun_radius_au = 0.00465047
        jupiter_radius_au = 0.10045 * sun_radius_au

        if ellipsoid:
            jupiter_radius_multiplier = 1.071374
        else:
            jupiter_radius_multiplier = 1.0

        jupiter_radius_au = jupiter_radius_au * 1.071374

        # Check if Epoch is given
        if epoch is not None:
            # Check types
            if isinstance(epoch, Epoch):
                # Calculate the position of jupiter in solar-spherical coordinates
                l, b, r = Jupiter.geometric_heliocentric_position(epoch)

                # alpha is the umbra defining angle
                alpha_cone_rad = atan(r / (sun_radius_au - jupiter_radius_au))

                # beta is the penumbra defying angle
                beta_cone_rad = atan(r / (sun_radius_au + jupiter_radius_au))

                # Compute distance of the sharpest pint behind jupiter in jupiter radii
                cone_vertex_jupiter_radii = tan(alpha_cone_rad)
                cone_beta_vertex_jupiter_radii = (-1 * jupiter_radius_multiplier * tan(beta_cone_rad))

                return alpha_cone_rad, cone_vertex_jupiter_radii, beta_cone_rad, cone_beta_vertex_jupiter_radii
            else:
                raise TypeError("Invalid input type")



if __name__ == "__main__":
    test_epoch = Epoch()
    test_epoch.set(2020, 1, 2, 12, 36, 0)

    res = Detection.check_phenomena(test_epoch)

    print()
