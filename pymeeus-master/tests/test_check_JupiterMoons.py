import logging
import sys

import pytest
import pytest_check as check

from pymeeus.Epoch import Epoch
from pymeeus.JupiterMoons import JupiterMoons

_logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="[%(asctime)s] %(levelname)s:%(name)s:%(message)s", datefmt="%Y-%m-%d %H:%M:%S")

EPOCH_1992_12_16_JD = Epoch(2448972.500685)
EPOCH_1992_12_16_UTC = Epoch(1992, 12, 16, utc=True)
EPOCH_1992_12_16 = Epoch(1992, 12, 16.000685)
"""Epoch used for the calculations in meeus chapter 44"""

TOL_3 = 1e-3
"""Tolerance threshold 0.001"""
TOL_4 = 1e-4
"""Tolerance threshold 0.0001"""
TOL_5 = 1e-5
"""Tolerance threshold 0.00001"""

SATELLITES = ['Io', 'Europe', 'Ganymede', 'Callisto']
COORDS_X_Y = ['X', 'Y']
COORDS_X_Y_Z = ['X', 'Y', 'Z']


def msg(idx_sat: int, idx_coord, tol_threshold: float, comp_value: float, exp_value: float):
    return f"{SATELLITES[idx_sat]}: {COORDS_X_Y[idx_coord]}: " \
           f"computed value {comp_value} NOT within tolerance threshold ({tol_threshold}) of expected value ({exp_value})"


class CheckJupiterMoons:  # class TestJupiterMoons(TestCase): - Annotation doesn't work with inheritance?!?

    @pytest.mark.parametrize('epoch', [EPOCH_1992_12_16, EPOCH_1992_12_16_UTC, EPOCH_1992_12_16_JD])
    @pytest.mark.parametrize('tol_threshold', [TOL_3, TOL_4, TOL_5])
    def test_rectangular_positions(self, epoch, tol_threshold):
        """This method tests the method rectangular_positions() that calculates the rectangular geocentric
         position of Jupiter's satellites for a given epoch, using the E5-theory. """

        """ correction of light time and perspective effect are applied as described in astronomical algorithm
        chapter 44 (page 313-314)"""

        exp_values = [[-3.4502, 0.2137],
                      [7.4418, 0.2753],
                      [1.2011, 0.5900],
                      [7.0720, 1.0291]]

        # io_corr_true, europe_corr_true, ganymede_corr_true, callisto_corr_true
        results: [[]] = JupiterMoons.rectangular_positions(epoch, do_correction=True)

        for i in range(len(SATELLITES)):
            for j in range(len(COORDS_X_Y)):
                check.less(abs(results[i][j] - exp_values[i][j]),
                           tol_threshold,
                           msg(i, j, tol_threshold, results[i][j], exp_values[i][j]))

    @pytest.mark.parametrize('epoch', [EPOCH_1992_12_16, EPOCH_1992_12_16_UTC, EPOCH_1992_12_16_JD])
    @pytest.mark.parametrize('tol_threshold', [TOL_3, TOL_4, TOL_5])
    def test_calculate_delta(self, epoch, tol_threshold):
        """This method tests calculate_delta() that calculates the distance between Earth and Jupiter (DELTA)
        for a given epoch by iteration. """
        delta, tau = JupiterMoons.calculate_DELTA(epoch)

        # exp(ected)_value: as taken from Meeus (TODO: insert reference)
        exp_value: float = 5.6611239
        comp_value = delta - exp_value
        check.less(abs(comp_value),
                   tol_threshold,
                   f"computed value of distance Earth-Jupiter ({comp_value}) "
                   f"NOT within tolerance threshold ({tol_threshold}) "
                   f"of expected value ({exp_value})")

    @pytest.mark.parametrize('epoch', [EPOCH_1992_12_16, EPOCH_1992_12_16_UTC, EPOCH_1992_12_16_JD])
    @pytest.mark.parametrize('tol_threshold', [TOL_3, TOL_4, TOL_5])
    def test_correct_rectangular_positions(self, epoch, tol_threshold):
        """This method tests the method correct_rectangular_positions() that corrects the rectangular geocentric
        position of Jupiter's satellites for a given epoch as described
        in astronomical algorithm chapter 44 page: 313 -314. """

        exp_values = [[0.0003, 0.0, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.0007, 0.0, 0.0]]

        """ calculate corrected and uncorrected rectangular positions for the Jupiter moons"""
        # io_corr_true, europe_corr_true, ganymede_corr_true, callisto_corr_true
        results_corr_true: [[]] = JupiterMoons.rectangular_positions(epoch, do_correction=True)
        results_corr_false: [[]] = JupiterMoons.rectangular_positions(epoch, do_correction=False)

        """calculate difference of corrected and uncorrected rectangular coordinates and compare  coordinate 
        for satellite I and satellite 4 with the maximum correction possible as described
        in astronomical algorithm chapter 44 page: 313."""
        for i in range(len(SATELLITES)):
            for j in range(len(COORDS_X_Y_Z)):
                delta = results_corr_true[i][j] - results_corr_false[i][j]
                check.less(abs(delta - exp_values[i][j]),
                           tol_threshold,
                           msg(i, j, tol_threshold, delta, exp_values[i][j]))
