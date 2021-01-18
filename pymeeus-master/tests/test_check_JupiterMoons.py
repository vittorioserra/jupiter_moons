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
EVENTS = ['OCC', 'ECC', 'N/A']

TOL_TEST_ALL = [TOL_3, TOL_4, TOL_5]
TOL_TEST_SIMPLE = [TOL_3, ]

EPOCH_TEST_ALL = [EPOCH_1992_12_16, EPOCH_1992_12_16_UTC, EPOCH_1992_12_16_JD]
EPOCH_TEST_SIMPLE = [EPOCH_1992_12_16, ]


def msg(idx_sat: int, idx_coord, tol_threshold: float, comp_value: float, exp_value: float):
    return f"{SATELLITES[idx_sat]}: {COORDS_X_Y[idx_coord]}: " \
           f"computed value {comp_value} NOT " \
           f"within tolerance threshold ({tol_threshold}) of expected value ({exp_value})"


@pytest.mark.parametrize('epoch', EPOCH_TEST_SIMPLE)
@pytest.mark.parametrize('tol_threshold', TOL_TEST_SIMPLE)
def test_rectangular_positions(epoch, tol_threshold):
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


@pytest.mark.parametrize('epoch', EPOCH_TEST_SIMPLE)
@pytest.mark.parametrize('tol_threshold', TOL_TEST_SIMPLE)
def test_calculate_delta(epoch, tol_threshold):
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


@pytest.mark.parametrize('epoch', EPOCH_TEST_SIMPLE)
@pytest.mark.parametrize('tol_threshold', TOL_TEST_SIMPLE)
def test_correct_rectangular_positions(epoch, tol_threshold):
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


def test_is_phenomena():
    """This method tests the method correct_rectangular_positions() that corrects the rectangular geocentric
    position of Jupiter's satellites for a given epoch as described
    in astronomical algorithm chapter 44 page: 313 -314. """

    """Matrix:
        col 0: test case identifier
        col 1: Epoch to be tested
        col 2/3: row/col-index of expected 'True'
    All other entries are expected to be 'False'!"""
    test_set = [["IO_OCC_START_2021_01_17_00_55_11", Epoch(2021, 1, 17.0383217592593), 0, 0],
                ["IO_ECC_START_2021_02_12_14_19_14", Epoch(2021, 2, 12.5966898148148), 0, 1],
                ["EU_OCC_START_2021_10_02_13_04_19", Epoch(2021, 10, 2.54466435185185), 1, 0],
                ["EU_ECC_START_2021_02_13_14_26_18", Epoch(2021, 2, 13.6015972222222), 1, 1],
                ["GA_OCC_START_2021_03_29_00_43_46", Epoch(2021, 3, 29.0303935185185), 2, 0],
                ["GA_ECC_START_2021_02_06_16_53_50", Epoch(2021, 2, 6.70405092592593), 2, 1],
                ["CA_OCC_START_2021_03_09_13_26_47", Epoch(2021, 3, 9.5602662037037), 3, 0],
                ["CA_ECC_START_2021_04_28_13_32_36", Epoch(2021, 4, 28.5643055555556), 3, 1],
                ["IO_OCC_END_2021_02_14_11_20_01", Epoch(2021, 2, 14.4722337962963), 0, 0],
                ["IO_ECC_END_2021_10_09_14_44_31", Epoch(2021, 10, 9.61424768518519), 0, 1],
                ["EU_OCC_END_2021_03_07_02_14_12", Epoch(2021, 3, 7.09319444444444), 1, 0],
                ["EU_ECC_END_2021_10_06_07_12_45", Epoch(2021, 10, 6.30052083333333), 1, 1],
                ["GA_OCC_END_2021_12_18_23_58_55", Epoch(2021, 12, 18.9992476851852), 2, 0],
                ["GA_ECC_END_2021_12_04_20_34_04", Epoch(2021, 12, 4.85699074074074), 2, 1],
                ["CA_OCC_END_2021_11_15_07_29_43", Epoch(2021, 11, 15.3123032407407), 3, 0],
                ["CA_ECC_END_2021_08_24_00_58_49", Epoch(2021, 8, 24.0408449074074), 3, 1]
                ]
    for test_label, epoch, row_idx, col_idx in test_set:
        result = JupiterMoons.is_phenomena(epoch)
        exp_event = False
        other_events: [[]] = []
        for i in range(len(result)):
            events = result[i]
            for j in range(len(events) - 1):
                if result[i][j]:
                    if i == row_idx and j == col_idx:
                        exp_event = True
                    else:
                        other_events.append([SATELLITES[i], EVENTS[j]])

        # assert success, "Test failed!"
        # check.is_true(success, f"Test failed for {epoch_label}")
        check.is_true(exp_event and len(other_events) == 0, f"Test={test_label}: "
                                                            f"'expected event'={exp_event}, "
                                                            f"'other event(s)'={other_events}")


if __name__ == "__main__":
    test_is_phenomena()
