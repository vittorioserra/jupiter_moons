import functools
from enum import unique, IntEnum
from typing import List

from pymeeus_optimized.Epoch import Epoch
from pymeeus_optimized.Jupiter import Jupiter
from pymeeus_optimized.JupiterMoons import JupiterMoons

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


class Phenomenom(Result):
    """Class that delivers calculation data for phenomena
    """

    def __init__(self, phenomenom_type: str, sat: int, distance: float, z: float, epoch: Epoch,
                 shadow_radius_umbra: float = 0,
                 shadow_radius_penumbra: float = 0):
        """Constructor for Phenomenom class.
        Sortable by Epoch

        :param phenomenom_type: Type of phenomenom:
            PA: Transit
            OC: Ocultation
            OM: Transit of the umbra of the satellite on the disk of Jupiter
            EC: Eclipse
        :type phenomenom_type: str
        :param sat: Number of satellite - 1
            0: Io
            1: Europa
            2: Ganymed
            3: Callisto
        :type sat: int
        :param distance: Perspectivic distance to center of Jupiter
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
        self.phenomenom_type = phenomenom_type
        self.sat = sat
        self.distance = distance
        self.z = z
        self.shadow_radius_umbra = shadow_radius_umbra
        self.shadow_radius_penumbra = shadow_radius_penumbra

        self.moon_radius = moons_radii[sat]

        self.phenomenom_occurs = False
        self.shadow_type: str = ""

        # Check shadow phenomenom type
        self._check_shadow_type()

    def get_critical_radius(self) -> float:
        """Method that calculates the perspective radius (distance) to Jupiter
        that is at least necessary to have an extreme shadow phenomenom.

        :returns: Perspective distance (in Jupiter's radii) to Jupiter for
            extreme shadow phenomena
        :rtype: float
        """

        # For eclipse return penumbra radius + moon radius
        if self.phenomenom_type == "EC":
            return self.shadow_radius_penumbra + self.moon_radius
        # For other phenomena return Jupiter's radius (= 1) + moon radius
        else:
            return 1 + self.moon_radius

    def has_right_z_sign(self) -> bool:
        """Determines whether the sign of the z-Coordinate corresponds with the phenomenom
        type. OC and EC only possible for more distant moons, OM and PA for moons that are
        less distant than Jupiter.

        :returns: Whether the sign of the z-Coordinate corresponds with the
            phenomenom type
        :rtype: bool
        """

        return (self.phenomenom_type in ["PA", "OM"] and self.z < 0) or (
                self.phenomenom_type in ["OC", "EC"] and self.z > 0)

    def update(self, epoch: Epoch) -> None:
        """Updates the Phenomenom object to a given Epoch

        :param epoch: Epoch the phenomenom shall be updated
        :type epoch: Epoch

        :rtype: None
        """

        # Calculate updated phenomenom
        updated_phenomenom = Detection.check_single_phenomena(epoch, self.sat, self.phenomenom_type)

        # Copy all class variables from updated phenomenom
        self.phenomenom_type = updated_phenomenom.phenomenom_type
        self.sat = updated_phenomenom.sat
        self.distance = updated_phenomenom.distance
        self.z = updated_phenomenom.z
        self.epoch = updated_phenomenom.epoch
        self.shadow_radius_umbra = updated_phenomenom.shadow_radius_umbra
        self.shadow_radius_penumbra = updated_phenomenom.shadow_radius_penumbra
        self.moon_radius = updated_phenomenom.moon_radius
        self.phenomenom_occurs = updated_phenomenom.phenomenom_occurs
        self.shadow_type = updated_phenomenom.shadow_type

    def _calc_apparent_moon_radius(self) -> None:
        """ Calculates the apparent moon radius (Jupiter's radii) that differs because
        of differing z-Coordinates

        :rtype: None
        """

        # 1 AU in Jupiter's radii
        jr_au = 2092.512039

        # Calculate distances to observer
        # Solar phenomenom
        if self.phenomenom_type in ["EC", "OM"]:
            distance_jup = JupiterMoons.calculate_DELTA(self.epoch)[0] * jr_au
        # Earth observated phenomenom
        else:
            distance_jup = Jupiter.geometric_heliocentric_position(self.epoch)[2] * jr_au
        distance_moon = distance_jup + self.z

        # Calculate observated Angle
        alpha = atan(moons_radii[self.sat] / distance_moon)

        # Calc observated radius in plane of Jupiter with alpha as observated angle
        self.moon_radius = distance_jup * tan(alpha)

    def _check_shadow_type(self) -> None:
        """Private method that determines shadow type and whether a
        phenomenom occurs

        :rtype: None
        """

        # Calculate apparent moon radius
        self._calc_apparent_moon_radius()

        # Check z-Coordinate if Phenonomenom is possible
        if self.has_right_z_sign():
            # Set flag
            self.phenomenom_occurs = True

            # Phenomena with EXT and INT
            if self.phenomenom_type in ["PA", "OC", "OM"]:
                # Exterior contact with shadow
                if 1 - self.moon_radius < self.distance <= 1 + self.moon_radius:
                    self.shadow_type = "EXT"
                # Interior contact with shadow
                elif self.distance <= 1 - self.moon_radius:
                    self.shadow_type = "INT"
                else:
                    # Unset flag, no phenomenom detected
                    self.phenomenom_occurs = False

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
                    # Unset flag, no phenomenom detected
                    self.phenomenom_occurs = False


class Detection:
    """This class delivers methods for checking for any kind of phenomena
    of Jupiter's four gallilean moons.
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
    def check_single_phenomena(epoch: Epoch, i_sat: int, phenomenom_type: str) -> Phenomenom:
        """Calculates a Phenomenom object for a given satellite and phenomenom
        type to a given Epoch

        :param epoch: Epoch that shall be checked
        :type epoch: Epoch
        :param i_sat: Index of satellite (i.e. No. of Sat. - 1) that shall
            be checked
        :type i_sat: int
        :param phenomenom_type: Type of the phenomenom being checked
        :type phenomenom_type: str

        :returns: Phenomenom object for the given inputs
        :rtype: Phenomenom

        :raises: TypeError, if no Epoch typed input is given
        """

        # Check input type
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")

        # Calc perspective coordinates of the satellites depending on phenomenom type
        Coords = JupiterMoons.rectangular_positions_jovian_equatorial(epoch, solar=phenomenom_type in ["EC", "OM"])

        # Earth based phenomena
        if phenomenom_type in ["PA", "OC"]:
            # Calc perspective distance to Jupiter for given satellite
            distance = Detection.perspective_distance(Coords[i_sat][0], Coords[i_sat][1])

            # Return Phenomena instance
            return Phenomenom(Detection.phenomena_types_str[Detection.phenomena_types_str.index(phenomenom_type)],
                              i_sat, distance, Coords[i_sat][2], epoch)

        # Solar phenomena
        else:
            # Calc perspective distance to Jupiter for given satellite
            distance = Detection.perspective_distance(Coords[i_sat][0], Coords[i_sat][1])
            # Calc shadow cone parameter
            r_umbra, r_penumbra = Detection.cone_radius(epoch, Coords[i_sat][2])

            # Return Phenomena instance
            return Phenomenom(Detection.phenomena_types_str[Detection.phenomena_types_str.index(phenomenom_type)],
                              i_sat, distance, Coords[i_sat][2], epoch, r_umbra,
                              r_penumbra)

    @staticmethod
    def check_phenomena(epoch: Epoch) -> List[List[Phenomenom]]:
        """Calculates Phenomenom objects for all satellite and phenomena
        types for a given Epoch

        :param epoch: Epoch that shall be checked
        :type epoch: Epoch

        :returns: 2D-Array of Phenomena object, where rows correspond with
            satellites, columns with phenomena types
        :rtype: List[List[Phenomenom]]

        :raises: TypeError, if no Epoch typed input is given
        """

        # Check input type
        if not isinstance(epoch, Epoch):
            raise TypeError("Invalid input type")

        # Result matrix, where each rows is for a satellite
        # Column 0: Occultation (OC)
        # Column 1: Eclipse (EC)
        # Column 2: Transit (PA)
        # Column 3: Satellite umbra transit (OM)
        # noinspection PyTypeChecker
        result_matrix: List[List[Phenomenom]] = [[None, None, None, None],
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
            for phenomenom_type in Detection.PhenomenaTypes:
                if phenomenom_type in [Detection.PhenomenaTypes.PA, Detection.PhenomenaTypes.OC]:
                    # Calc perspective distance to Jupiter
                    distance = Detection.perspective_distance(Coords_Earth[i_sat][0], Coords_Earth[i_sat][1])

                    # Fill result matrix with Phenomena instance
                    result_matrix[i_sat][phenomenom_type] = Phenomenom(Detection.phenomena_types_str[phenomenom_type],
                                                                       i_sat, distance, Coords_Earth[i_sat][2], epoch)
                else:
                    # Calc perspective distance to Jupiter
                    distance = Detection.perspective_distance(Coords_Sun[i_sat][0], Coords_Sun[i_sat][1])
                    # Calc shadow cone parameter
                    r_umbra, r_penumbra = Detection.cone_radius(epoch, Coords_Sun[i_sat][2])

                    # Fill result matrix with Phenomena instance
                    result_matrix[i_sat][phenomenom_type] = Phenomenom(Detection.phenomena_types_str[phenomenom_type],
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
            # Distord y-coordinate to model elliptical Jupiter
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
            jupiter_radius_multiplicator = 1.071374
        else:
            jupiter_radius_multiplicator = 1.0

        jupiter_radius_au = jupiter_radius_au * 1.071374

        # Check if Epoch is given
        if epoch is not None:
            # Check types
            if isinstance(epoch, Epoch):
                # Calcuate the position of jupiter in solar-spherical coordinates
                l, b, r = Jupiter.geometric_heliocentric_position(epoch)

                # alpha is the umbra defining angle
                alpha_cone_rad = atan(r / (sun_radius_au - jupiter_radius_au))

                # beta is the penumbra defing angle
                beta_cone_rad = atan(r / (sun_radius_au + jupiter_radius_au))

                # Compute distance of the sharpest pint behind jupiter in jupiter radii
                cone_vertex_jupiter_radii = tan(alpha_cone_rad)
                cone_beta_vertex_jupiter_radii = (-1 * jupiter_radius_multiplicator * tan(beta_cone_rad))

                return alpha_cone_rad, cone_vertex_jupiter_radii, beta_cone_rad, cone_beta_vertex_jupiter_radii
            else:
                raise TypeError("Invalid input type")


if __name__ == "__main__":
    test_epoch = Epoch()
    test_epoch.set(2020, 1, 2, 12, 36, 0)

    res = Detection.check_phenomena(test_epoch)

    print()
