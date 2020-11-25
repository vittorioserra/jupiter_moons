from unittest import TestCase
from ..pymeeus.JupiterMoons import JupiterMoons
from pymeeus.base import TOL
from pymeeus.Epoch import Epoch


class TestJupiterMoons(TestCase):
    def test_rectangular_positions(self):
        """This method tests the method rectangular_positions() that calculates the rectangular geocentric
         position of Jupiter's satellites for a given epoch, using the E5-theory. """

        """ correction of light time and perspective effect are applied as described in astronomical algorithm
        chapter 44 (page 313-314)"""

        """Epoch used for the calculations in meeus chapter 44"""
        epoch = Epoch(1992, 12, 16)


        io_corr_true, europe_corr_true, ganymede_corr_true, callisto_corr_true = JupiterMoons.rectangular_positions(epoch, do_correction=True)

        assert abs(round(io_corr_true[0], 4) - 3.4502) < TOL, \
            """ERROR: 1st rectangular position (X) for Io of JupiterMoons.rectangular_position() test doesn't match"""

        assert abs(round(io_corr_true[1], 4) - 0.2137) < TOL, \
            """ERROR: 2nd rectangular position (Y) for Io of JupiterMoons.rectangular_position() test doesn't match"""

        assert abs(round(europe_corr_true[0], 4) - 7.4418) < TOL, \
            """ERROR: 1st rectangular position (X) for Europe of JupiterMoons.rectangular_position()
            test doesn't match"""

        assert abs(round(europe_corr_true[1], 4) - 0.2753) < TOL, \
            """ERROR: 2nd rectangular position for (Y) Europe of JupiterMoons.rectangular_position()
            test doesn't match"""

        assert abs(round(ganymede_corr_true[0], 4) - 1.2011) < TOL, \
            """ERROR: 1st rectangular position (X) for Ganymede of JupiterMoons.rectangular_position()
            test doesn't match"""

        assert abs(round(ganymede_corr_true[1], 4) - 0.5900) < TOL, \
            """ERROR: 1st rectangular position (X) for Ganymede of JupiterMoons.rectangular_position()
            test doesn't match"""

        assert abs(round(callisto_corr_true[0], 4) - 7.0720) < TOL, \
            """ERROR: 1st rectangular position (X) for Callisto of JupiterMoons.rectangular_position() 
            test doesn't match"""

        assert abs(round(callisto_corr_true[1], 4) - 1.0291) < TOL, \
            """ERROR: 2nd rectangular position for (Y) Callisto of JupiterMoons.rectangular_position()
            test doesn't match"""

    def test_calculate_delta(self):
        """This method tests calculate_delta() that calculates the distance between Earth and Jupiter (DELTA)
        for a given epoch by iteration. """

        """Epoch used for the calculations in meeus chapter 44"""
        epoch = Epoch(1992, 12, 16)

        delta, tau = JupiterMoons.calculate_DELTA(epoch)

        assert abs(round(delta, 4) - 5.6612) < TOL, \
            """ERROR: Distance between earth and Jupiter of JupiterMoons.calculate_DELTA()
            doesn't match"""


    def test_correct_rectangular_positions(self):
        """This method tests the method correct_rectangular_positions() that corrects the rectangular geocentric
        position of Jupiter's satellites for a given epoch as described
        in astronomical algorithm chapter 44 page: 313 -314. """

        epoch = Epoch(1992, 12, 16)

        """ calculate corrected rectangular positions for the Jupiter moons"""
        io_corr_true, europe_corr_true, ganymede_corr_true, callisto_corr_true = JupiterMoons.rectangular_positions(epoch, do_correction=True)

        """ calculate uncorrected rectangular positions for the Jupiter moons"""
        io_corr_false, europe_corr_false, ganymede_corr_true_corr_false, callisto_corr_false = JupiterMoons.rectangular_positions(epoch, do_correction=False)

        """calculate difference of corrected and uncorrected rectangular coordinates and compare  coordinate 
        for satellite I and satellite 4 with the maximum correction possible as described
        in astronomical algorithm chapter 44 page: 313."""
        x_correction_io = io_corr_true[0] - io_corr_false[0]
        x_correction_callisto = callisto_corr_true[0] - callisto_corr_false[0]
        y_correction_io = io_corr_true[1] - io_corr_false[1]
        y_correction_callisto = callisto_corr_true[1] - callisto_corr_false[1]
        z_correction_io = io_corr_true[2] - io_corr_false[2]
        z_correction_callisto = callisto_corr_true[2] - callisto_corr_false[2]

        assert abs(round(x_correction_io, 4) - 0.0003) < TOL, \
            """ERROR: correction of X coordinate of Io out of range"""
        assert abs(round(x_correction_callisto, 4) - 0.0007) < TOL, \
            """ERROR: correction of X coordinate of Callisto out of range"""
        assert abs(round(y_correction_io, 4)) < TOL, \
            """ERROR: correction of Y coordinate of Io out of range"""
        assert abs(round(y_correction_callisto, 4)) < TOL, \
            """ERROR: correction of Y coordinate of Callisto out of range"""
        assert z_correction_io == 0, \
            """ERROR: correction of Z coordinate of Io out of range"""
        assert z_correction_callisto == 0, \
            """ERROR: correction of Z coordinate of Callisto out of range"""



    def test_check_coordinates(self):
        """This method test if the method check_coordinates() returns a distance in Jupiter radi less then one
        for a phenomena."""

        """taking the ellipsoidal shape of Jupiter into account through multiplication Y coordinates
        in check_phenomena with 1.071372 """

        y_correct_ellipsoid = 1.071374

        """limits of perspective distance to Jupiter's center in Jupiter's radii for a phenomena"""
        x_start_phenomena = 1
        y_start_phenomena = 1/y_correct_ellipsoid

        """calculate perspective distance to Jupiter's center in Jupiter's radii"""
        r_distance_center = JupiterMoons.check_coordinates(x_start_phenomena, y_start_phenomena)

        assert r_distance_center == 1, \
            """ERROR:  wrong calculated perspective distance to Jupiter's center in Jupiter's radii"""



    def test_check_occulation(self):
        self.fail()

    def test_check_eclipse(self):
        self.fail()


    def test_check_phenomena(self):
        self.fail()

    def test_is_phenomena(self):
        self.fail()
