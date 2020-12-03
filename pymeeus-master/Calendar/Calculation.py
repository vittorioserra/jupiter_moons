# -*- coding: utf-8 -*-

import functools
import multiprocessing
from math import floor

from typing import List, Tuple

from Calendar.PhenomenaDetection import Detection, Phenomenom, Result

from pymeeus.Epoch import Epoch
from pymeeus.JupiterMoons import JupiterMoons

# 1 second in jd = 1.157401129603386e-05
s_jd = 1.157401129603386e-05


class CalculationResult(Result):
    """ This class handles results calculated by the pymeeus
    JupiterMoons.rectangular_positions() method.
    It supports sorting the results by date of the calculation
    """

    def __init__(self, epoch: Epoch, dist_matrix: list):
        """ Constructor for Calculation results class

        :param epoch: Epoch the Calculation was made for
        :type epoch: Epoch
        :param dist_matrix: perspective distance to phenomena (column)
            for all satellites (rows) (if < 1: phenomena taking place)
        :type dist_matrix: list
        """

        # Call Result constructor
        super().__init__(epoch)

        # Initiate variables
        self.dist_matrix = dist_matrix


@functools.total_ordering
class TimingResult(Result):
    """ This class handles timing results. It supports sorting the
    results by date of the calculation.
    """

    def __init__(self, epoch: Epoch, appereance_type: str, row: int, col: int, rough_time_step: int,
                 phenomenom: Phenomenom):
        """ Constructor for TimingResult class

        :param epoch: Exact Epoch of phenomenom, could differ from
            phenomenom.epoch if time is corrected
        :type epoch: Epoch
        :param appereance_type: Whether it is start or end of phenomenom
        :type: str
        :param row: Corresponding row in the timing_lists array (i.e.
            row = Number of satellite - 1)
        :type row: int
        :param col: Corresponding column in the timing_lists array (i.e.
            col = 0: Ocultation, col = 1: Eclipse, col = 2: Penumbra)
        :type col: int
        :param phenomenom: Phanomenom that occurs
        :type phenomenom: Phenomenom
        """

        # Call Result constructor
        super().__init__(epoch)

        # Initiate variables
        self.appereance_type = appereance_type
        self.row = row
        self.col = col
        self.rough_time_step = rough_time_step
        self.phenomenom = phenomenom

        self.is_time_corrected = False

    def find_other_shadow_phenomena(self) -> list:
        """Method that searches for the remaining shadow phenomena for a
        given start or end of a phenomenom extreme shadow types.
        For Eclipse: PEN
        For other phenomena types: EXT

        :returns: List of TimingResults with starts or end of the remaining
            shadow phenomena
        :rtype: List[TimingResult]
        """

        # List for other timings
        other_timings = []

        # Switch appearing shadow phenomena
        if self.phenomenom.phenomenom_type == "EC":
            phenomenom_flag = ["PEN", "EXT", "INT"]
        else:
            phenomenom_flag = ["EXT", "INT"]

        # Time step for search
        time_step = s_jd * 2

        # If timing is start of phenomenom
        if self.appereance_type == "start":
            # Add self to list as start of extrem shadow phenomenom
            other_timings.append(self)

            # Set starting Phenomenom(Epoch) for iteration
            cur_phenomenom = self.phenomenom

            # Iterate through shadow types start of first shadow phenomenom already found -> start at index 1
            for i in range(1, len(phenomenom_flag)):
                # While shadow type is not reached -> calculate for next time step
                while not (cur_phenomenom.shadow_type == phenomenom_flag[i] and cur_phenomenom.phenomenom_occurs):
                    cur_phenomenom = Detection.check_single_phenomena(cur_phenomenom.epoch + time_step,
                                                                      cur_phenomenom.sat,
                                                                      cur_phenomenom.phenomenom_type)
                # Add Phenomenom to list
                other_timings.append(
                    TimingResult(cur_phenomenom.epoch, self.appereance_type, self.row, self.col,
                                 self.rough_time_step, cur_phenomenom))

        # If timing is start of phenomenom
        elif self.appereance_type == "end":

            # Set starting Phenomenom(Epoch) for iteration
            cur_phenomenom = self.phenomenom

            # Iterate through shadow types start of first shadow phenomenom already found -> start at index 1
            for i in range(1, len(phenomenom_flag)):
                # While shadow type is not reached -> calculate for previous time step
                while not (cur_phenomenom.shadow_type == phenomenom_flag[i] and cur_phenomenom.phenomenom_occurs):
                    cur_phenomenom = Detection.check_single_phenomena(cur_phenomenom.epoch - time_step,
                                                                      cur_phenomenom.sat,
                                                                      cur_phenomenom.phenomenom_type)
                # Add Phenomenom to list
                other_timings.append(
                    TimingResult(cur_phenomenom.epoch, self.appereance_type, self.row, self.col,
                                 self.rough_time_step,
                                 Detection.check_single_phenomena(cur_phenomenom.epoch,
                                                                  cur_phenomenom.sat,
                                                                  cur_phenomenom.phenomenom_type)))
            # Add self TimingResult - 1 timestep as last appereance of extrem shadow phenomenom to list
            other_timings.append(TimingResult(self.epoch - time_step, self.appereance_type, self.row, self.col,
                                              self.rough_time_step,
                                              Detection.check_single_phenomena(self.epoch - time_step,
                                                                               self.phenomenom.sat,
                                                                               self.phenomenom.phenomenom_type)))

        return other_timings

    def time_correction(self) -> None:
        """Corrects a bug in pymeeus for the differential Jupiter-Earth
        light travel time when computing eclipses

        :rtype: None
        """
        jr_au = 2092.512039

        if self.phenomenom.phenomenom_type in ["EC", "OM"] and self.is_time_corrected is False:
            # Compute Distance Jupiter - Earth and Coordinates of Jupiter
            dist = JupiterMoons.calculate_DELTA(self.epoch)

            # DELTA - r_jupiter
            correction_term_au = dist[0] + self.phenomenom.z / jr_au - dist[4]

            # Light time delay for the distance
            time_correction = correction_term_au * 149597870700 / 299792458  # (150*(10**9))/(3*(10**8))

            time_correction_jd = s_jd * time_correction

            # Correct Epoch
            self.epoch = self.epoch + time_correction_jd
            self.is_time_corrected = True
        elif self.is_time_corrected is False:
            correction_term_au = self.phenomenom.z / jr_au

            # Light time delay for the distance
            time_correction = correction_term_au * 149597870700 / 299792458  # (150*(10**9))/(3*(10**8))

            time_correction_jd = s_jd * time_correction

            # Correct Epoch
            self.epoch = self.epoch + time_correction_jd
            self.is_time_corrected = True

    def __lt__(self, other):
        """Compares Timing result to an other Result depending
        on the Epoch (if set) the calculation was made for.
        Otherwise it will be sorted by rough time step.

        :param other: Other result
        :type other: Any

        :returns: Whether other TimingResult has higher Epoch or not
        :rtype: bool

        :raises: TypeError, if other is no instance of TimingResult
        """

        # Check type
        if not isinstance(other, TimingResult):
            raise TypeError("Can't compare different types")

        # If Epoch is not set
        if self.epoch is None:
            return self.rough_time_step < other.rough_time_step

        # If Epoch is set, call overwritten method
        return super().__lt__(other)

    def __eq__(self, other):
        """ Checks if TimingResult is equal to an other TimingResult
        depending on the Epoch of the timing. If not set, the
        TimingResult will be compared by rough time step.

        :param other: Other Result
        :type other: Any

        :returns: Whether other TimingResult has same Epoch or not
        :rtype: bool

        :raises: TypeError, if other is no instance of Result
        """

        # Check type
        if not isinstance(other, TimingResult):
            raise TypeError("Can't compare different types")

        # If Epoch is not set
        if self.epoch is None:
            return self.rough_time_step == other.rough_time_step

        # If Epoch is set, call overwritten method
        return super().__eq__(other)


class Calculation:
    """ Class that handles Calculation of phenomena for a given timespan.
    First calculates rough timings of begin and end of an phenomena for
    the given time step, then increases accuracy by iterating around the
    rough timings using binary search.
    """

    # Core count of calculating CPU
    cpu_core_count = 8

    def __init__(self, start_epoch: Epoch, end_epoch: Epoch, time_step: float, tol: float = 0.0):
        """ Constructor for Calculation class

        :param start_epoch: Epoch to start calculating
        :type start_epoch: Epoch
        :param end_epoch: Epoch to end calculation
        :type end_epoch: Epoch
        :param time_step: Rough time step in julian date
        :type time_step: float
        :param tol: Tolerance for detecting phenomena (0: rough timing
            triggers if distance to phenomena < 1 [Jupiter's radii])
        :type tol: float

        :raises: ValueError, if start is later than end of calculation
        """

        # Check start_epoch and end_epoch
        if start_epoch > end_epoch:
            raise ValueError("Start is later than End of Calculation")

        # Initialize variables
        self.start_epoch = start_epoch
        self.end_epoch = end_epoch
        self.time_step = time_step
        self.tol = tol
        self.number_of_timings = 0
        self.all_timings_sorted: List[TimingResult] = []

        # List for all CalculationResults during the timespan
        self.calc_res: List[CalculationResult] = []

        # Array full of lists for the distances over time
        # Rows: satellites
        # Columns: Phenomena (0: Ocultation, 1: Eclipse, 2: Transit, 3: Shadow transit)
        self.distances = [[[], [], [], []],
                          [[], [], [], []],
                          [[], [], [], []],
                          [[], [], [], []]]

        # Array full of lists for the timings of the phenomena
        self.timing_lists = [[[], [], [], []],
                             [[], [], [], []],
                             [[], [], [], []],
                             [[], [], [], []]]

        # Start calculation
        self.calculate(self.start_epoch, self.end_epoch, self.time_step, self.tol)

    def calculate(self, start_epoch: Epoch, end_epoch: Epoch, time_step: float, tol: float = 0.0) -> None:
        """ This method first calculates the perspective coordinates
        for the given timespan with the given time step, then finds
        the rough timings depending on the course of the perspective
        coordinates. Later on, the rough timings are used to calculate
        the exact timings that will be stored in timing_lists.

        :param start_epoch: Epoch the calculation starts
        :type start_epoch: Epoch
        :param end_epoch: Epoch the calculation stops
        :type end_epoch: Epoch
        :param time_step: Time step for rough calculation in julian date
        :type time_step: float
        :param tol: Tolerance for detecting phenomena (0: rough timing
            triggers if distance to phenomena < 1 [Jupiter's radii])
        :type time_step: float

        :rtype: None
        """

        # Calc time span each process has to calculate
        time_span_per_core = (end_epoch - start_epoch) / self.cpu_core_count

        # Set up Queue
        procs = []
        queue_cal = multiprocessing.Queue()

        # Make sure all calculations have same time step
        desired_epochs = []
        actual_epochs = []

        # Calculate desired time spans for all processes
        for i in range(self.cpu_core_count):
            desired_epochs.append((start_epoch + i * time_span_per_core, start_epoch + (i + 1) * time_span_per_core))

        last_end_epoch = desired_epochs[0][0]

        # Set up actual time span for all processes
        for i in range(len(desired_epochs)):
            actual_start = last_end_epoch
            actual_end = last_end_epoch

            # Move end until desired end is reached
            while actual_end < desired_epochs[i][1]:
                actual_end += time_step

            last_end_epoch = actual_end

            actual_epochs.append((actual_start, actual_end))

        # Start process for each core
        for i in range(len(actual_epochs)):
            proc = multiprocessing.Process(target=self.mp_calculate_position, args=(
                queue_cal, actual_epochs[i][0], actual_epochs[i][1], time_step))
            procs.append(proc)
            proc.start()

        results = []

        # Number of already finished processes
        processes_finished = 0

        # Read Queue until all processes have finished
        while processes_finished < self.cpu_core_count:
            res = queue_cal.get()
            results.append(res)
            processes_finished += 1

        # Wait for remaining processes to finish
        for i in procs:
            i.join()

        # Add results to calculation
        self.calc_res.clear()
        for lst in results:
            self.calc_res.extend(lst)

        # Calculate coordinates (distances) in respect to phenomena
        # Fill distances array
        self.list_distances()

        # Clear timing_list
        self.timing_lists = self.apply_to_2d_array(lambda x: [], self.timing_lists)

        # Getting rough timings for the phenomena
        # Iterate through distances array
        for row in range(len(self.distances)):
            for col in range(len(self.distances[row])):
                # Find phenomena for current list in distances[row][col]
                self.timing_lists[row][col].extend(self.find_phenomena(self.distances[row][col], row, col, tol))
        print("Got rough timings")

        # Set up Queue for comunication between Processes
        queue_bs = multiprocessing.Queue()
        # List for processes
        procs = []

        # Get exact timings
        # Iteration through timing_lists
        for row in range(len(self.timing_lists)):
            for col in range(len(self.timing_lists[row])):
                # Set up process
                proc = multiprocessing.Process(target=self.mp_binarysearch,
                                               args=(queue_bs, self.timing_lists[row][col], row, col))
                # Start process
                procs.append(proc)
                proc.start()

        results = []

        # Number of finished processes
        processes_finished = 0

        # Read Queue while not all processes have finished
        while processes_finished < len(procs):
            res = queue_bs.get()
            results.append(res)
            processes_finished += 1

        # Wait for all processes to finish
        for i in procs:
            i.join()

        # Clear timing_lists
        self.timing_lists = Calculation.apply_to_2d_array(lambda x: [], self.timing_lists)

        # Rearrange results in timing_lists
        for lst in results:
            for timing in lst:
                self.timing_lists[timing.row][timing.col].append(timing)

        print("Got exact timings")

        # Set up sorted list for timings
        extreme_timings: List[TimingResult] = []

        for row in self.timing_lists:
            for lst in row:
                extreme_timings.extend(lst)

        # Sort timings by epoch
        extreme_timings = sorted(extreme_timings)

        for i in range(len(extreme_timings)):
            # Update Phenomena for exact epochs
            extreme_timings[i].phenomenom.update(extreme_timings[i].epoch)

        # Set up multiprocessing for finding other shadow phenomena
        # Count of processes double the cpu core count to better distribute work
        # when processes with less work already have finished
        proc_count = self.cpu_core_count * 2
        calcs_per_core = floor(len(extreme_timings) / proc_count)

        # Set up Queue for comunication between Processes
        queue_sst = multiprocessing.Queue()
        # List for processes
        procs = []

        # If there are too less timings for all processes
        if calcs_per_core <= 1:
            indices = [(0, len(extreme_timings) - 1)]
        else:
            # List for start and end index any process has to calculate
            indices = [(0, calcs_per_core - 1)]

            for i in range(1, proc_count - 1):
                indices.append((indices[i - 1][1] + 1, (indices[i - 1][1]) + calcs_per_core))

            indices.append((indices[proc_count - 2][1] + 1, len(extreme_timings) - 1))

        print("Search shadow phenomena")
        # Start process for given list sector
        for inds in indices:
            timing_list = extreme_timings[inds[0]:inds[1] + 1]
            proc = multiprocessing.Process(target=self.mp_search_shadow_timings,
                                           args=(queue_sst, timing_list))
            procs.append(proc)
            proc.start()

        # Number of finished processes
        processes_finished = 0

        results = []
        # Read Queue while not all processes have finished
        while processes_finished < len(procs):
            res = queue_sst.get()
            results.extend(res)
            processes_finished += 1
            print(processes_finished, "of", len(procs), "Shadow search process finished")

        # Wait for all processes to finish
        for i in procs:
            i.join()

        # Do time correction
        for res in results:
            res.time_correction()

        # Trim and sort results
        self.all_timings_sorted = sorted(self._cut_to_timespan(results))

        self.number_of_timings = len(self.all_timings_sorted)
        print("Found ", self.number_of_timings, " timings")
        print()

        return

    @staticmethod
    def mp_search_shadow_timings(queue: multiprocessing.Queue, timing_list: List[TimingResult]):
        """Method for handling multiprocessing for calculating other shadow
        phenomena for a given list of start or end timings of extreme shadow
        phenomena.

        :param queue: Queue for communication with other processes
        :type queue: multiprocessing.Queue
        :param timing_list: List of start or end timings of extreme shadow
            phenomena
        :type timing_list: List[TimingResult]

        :rtype: None
        """

        res = []

        # Calculate other shadow phenomena
        for timing in timing_list:
            res.extend(timing.find_other_shadow_phenomena())

        # Write result to queue
        queue.put(res)

    @staticmethod
    def mp_calculate_position(queue: multiprocessing.Queue, start_epoch: Epoch, end_epoch: Epoch,
                              time_step: float):
        """Method that handels a multiprocessing task for calculating coordinates
        for Epochs in the given timespan with given time step. The results will
        be written to the queue.

        :param queue: Queue for communication with other processes
        :type queue: multiprocessing.Queue
        :param start_epoch: Start Epoch for calculation
        :type start_epoch: Epoch
        :param end_epoch: Start Epoch for calculation
        :type end_epoch: Epoch
        :param time_step: Time step for the calculation
        :type time_step: float

        :rtype: None
        """

        # Write calculation result to Queue for communication with other processes
        queue.put(Calculation.make_calculation(start_epoch, end_epoch, time_step))

    def mp_binarysearch(self, queue: multiprocessing.Queue, timing_list: list, row: int, col: int) -> None:
        """Method that handels task (binary search to find exact timings)
        of one process in order to implement multiprocessing to use multiple
        CPU cores for the calculation.

        :param queue: Queue for communication with other processes
        :type queue: multiprocessing.Queue
        :param timing_list: List calculation of exact timings will be made
        :type timing_list: list
        :param row: Row of the timing_list in timing_lists
        :type row: int
        :param col: Column of the timing_list in timing_lists
        :type col: int

        :rtype: None
        """

        result: List[TimingResult] = []

        # Iteration through given timing_list
        for i in range(len(timing_list)):
            # Get current content of the list
            timing: TimingResult = timing_list[i]
            # Calculate and write exact timing together with former content and array coordinates
            timing.epoch = self.find_start_end_cons_tol(self.calc_res[timing.rough_time_step].epoch, self.time_step,
                                                        row, col, timing.appereance_type, self.tol)
            result.append(timing)

        # Add result to Queue
        queue.put(result)

        print("Exact timings for Satellite ", row + 1, ", phenomenom: ", col, " calculated")

    @staticmethod
    def make_calculation(start_epoch: Epoch, end_epoch: Epoch, time_step: float) -> List[CalculationResult]:
        """ Calculates distances to phenomena for all satellites in the
        given timespan with the given time step

        :param start_epoch: Epoch the calculation starts
        :type start_epoch: Epoch
        :param end_epoch: Epoch the calculation stops
        :type end_epoch: Epoch
        :param time_step: Time step for rough calculation in julian date
        :type time_step: float

        :returns: List of Calculation Results
        :rtype: List[CalculationResult]
        """

        # Set up return variable
        dist_list: List[CalculationResult] = []

        # Set current Epoch for iteration
        curr_epoch = start_epoch

        # Loop as long as end_epoch is not reached
        while curr_epoch.jde() <= end_epoch.jde():
            # Add CalculationResult for the calculation result and current Epoch to result
            dist_list.append(CalculationResult(curr_epoch, Detection.check_phenomena(curr_epoch)))
            # Increase current Epoch by time step
            curr_epoch += time_step
        print("Coordinate Calculation finished for: ", start_epoch, " - ", end_epoch)

        return dist_list

    def list_distances(self) -> None:
        """Fills distances list with respect to the CalculationResults
        stored in calc_res. This method rearranges the data to be able
        to see the course of the coordinates over time in a seperate
        list in the distances array.

        :rtype: None
        """

        # Sort calc_res by date
        self.calc_res = sorted(self.calc_res)

        # Clear lists in distances array
        self.distances = self.apply_to_2d_array(lambda x: [], self.distances)

        # Iterate through CalculationResults in sorted calc_res list
        for ele in self.calc_res:
            # Iterate through distance array
            for row in range(len(ele.dist_matrix)):
                for col in range(len(ele.dist_matrix[row])):
                    # Fill distance list
                    self.distances[row][col].append(ele.dist_matrix[row][col])

    def find_phenomena(self, distance_list: List[Phenomenom], row: int, col: int, tol: float = 0.0) -> List[
        TimingResult]:
        """This method finds rough timings of start and end of phenomena
        for a given course of distances.

        :param distance_list: List of distances for concurrent time steps
        :type distance_list: list
        :param row: Row of distance_list in distances array
        :type row: int
        :param col: Column of distance_list in distances array
        :type col: int
        :param tol: Tolerance for detecting phenomena (Jupiter's radii)
        :type tol: float

        :returns: List of Tuples with time step of first appereance and
            type of timing ("start" or "end") of a phenomena
        :rtype: List[Tuple[int, str]]
        """

        first_entry = distance_list[0]

        # Set flag
        if first_entry.distance > first_entry.get_critical_radius() + tol or not first_entry.has_right_z_sign():
            is_phenomena = False
        else:
            result_timestep_earlier = Detection.check_phenomena(self.start_epoch - self.time_step)[row][col]
            is_phenomena = result_timestep_earlier.distance <= result_timestep_earlier.get_critical_radius() + tol and result_timestep_earlier.has_right_z_sign()

        # Set up return variable
        timing_list: List[TimingResult] = []

        # Iterate through distance_list
        for i in range(len(distance_list)):
            # Check for start of a phenomena
            if is_phenomena is False and distance_list[i].has_right_z_sign() and distance_list[i].distance <= \
                    distance_list[i].get_critical_radius() + tol:
                timing_list.append(TimingResult(None, "start", row, col, i, distance_list[i]))
                is_phenomena = True
            # Check for end of phenomena
            elif is_phenomena is True and not (
                    distance_list[i].has_right_z_sign() and distance_list[i].distance <= distance_list[
                i].get_critical_radius() + tol):
                timing_list.append(TimingResult(None, "end", row, col, i, distance_list[i]))
                is_phenomena = False

        return timing_list

    @staticmethod
    def apply_to_2d_array(func, array) -> tuple:
        """ This method applies a function to all elements of an two-
        dimensional array

        :param func: Function that will be applied
        :type func: lambda
        :param array: 2d-Array the function will be applied to
        :type array: tuple, list

        :returns: Manipulated array
        :rtype: tuple, list
        """

        # Copy input array
        res = array.copy()

        # Iterate through rows
        for row in range(len(array)):
            # Apply function on row and store in result variable
            res[row] = tuple(map(func, array[row]))

        return res

    @staticmethod
    def find_start_end_cons_tol(first_appereance: Epoch, time_step: float, row: int, col: int, appereance_type: str,
                                tol: float = 0.0) -> Epoch:
        """This method implements binary search for finding
        exact Epochs for given rough timings of phenomena. It changes
        behavior considering whether there was a tolerance for finding rough
        timings.

        :param first_appereance: Rough Epoch the phenomenom (start or
            end) first appeared
        :type first_appereance: float
        :param time_step: Time step the calculation was done
        :type time_step: float
        :param row: Corresponding row in the timing_lists array (i.e.
            row = Number of satellite - 1)
        :type row: int
        :param col: Corresponding column in the timing_lists array (i.e.
            col = 0: Ocultation, col = 1: Eclipse, col = 2: Penumbra)
        :type col: int
        :param appereance_type: Which type of appereance was detected
            ("start" or "end")
        :type appereance_type: str
        :param tol: Tolerance for finding rough timings
        :type tol: float

        :returns: Exact Epoch for begin or end of the phenomenom
        :rtype: Epoch

        :raises: ValueError, if appereance_type is invalid
        """

        if tol == 0.0:
            # Regular behavior if no tolerance
            return Calculation.find_start_end(first_appereance, time_step, row, col, appereance_type)
        else:
            # Safe search time step to find when the phenomenom really occurs
            search_time_step = 60 * 60 * s_jd

            # Condition for phenomenom start
            if appereance_type == "start":
                # If there is a phenomenom in the given time step
                if Detection.check_phenomena(first_appereance)[row][col].phenomenom_occurs and not \
                        Detection.check_phenomena(first_appereance - time_step)[row][col].phenomenom_occurs:
                    return Calculation.find_start_end(first_appereance, time_step, row, col, appereance_type)

                cur_epoch = first_appereance

                # Whether there is an phenomenom at the end of the first appereance time step
                is_phenomena_end = Detection.check_phenomena(cur_epoch)[row][col].phenomenom_occurs

                # Move time step until there is a phenomenom
                while not is_phenomena_end:
                    cur_epoch += search_time_step
                    is_phenomena_end = Detection.check_phenomena(cur_epoch)[row][col].phenomenom_occurs

                return Calculation.find_start_end(cur_epoch, search_time_step, row, col, appereance_type)

            # Condition for phenomenom end
            elif appereance_type == "end":
                # If there is a phenomenom in the given time step
                if not Detection.check_phenomena(first_appereance)[row][col].phenomenom_occurs and \
                        Detection.check_phenomena(first_appereance - time_step)[row][col].phenomenom_occurs:
                    return Calculation.find_start_end(first_appereance, time_step, row, col, appereance_type)

                cur_epoch = first_appereance - time_step

                # Whether there is an phenomenom at the start of the first appereance time step
                is_phenomena_start = Detection.check_phenomena(cur_epoch)[row][col].phenomenom_occurs

                # Move time step until there is a phenomenom
                while not is_phenomena_start:
                    cur_epoch -= search_time_step
                    is_phenomena_start = Detection.check_phenomena(cur_epoch)[row][col].phenomenom_occurs

                return Calculation.find_start_end(cur_epoch + search_time_step, search_time_step, row, col,
                                                  appereance_type)

            # Input string invalid
            else:
                raise ValueError("Input string invalid")

    @staticmethod
    def find_start_end(first_appereance: Epoch, time_step: float, row: int, col: int, appereance_type: str) -> Epoch:
        """ This recursive method implements binary search for finding
        exact Epochs for given rough timings of phenomena

        :param first_appereance: Rough Epoch the phenomenom (start or
            end) first appeared
        :type first_appereance: float
        :param time_step: Time step the calculation was done
        :type time_step: float
        :param row: Corresponding row in the timing_lists array (i.e.
            row = Number of satellite - 1)
        :type row: int
        :param col: Corresponding column in the timing_lists array (i.e.
            col = 0: Ocultation, col = 1: Eclipse, col = 2: Penumbra)
        :type col: int
        :param appereance_type: Which type of appereance was detected
            ("start" or "end")
        :type appereance_type: str

        :returns: Exact Epoch for begin or end of the phenomenom
        :rtype: Epoch

        :raises ValueError, if appeareanceType is neither "start"
            nor "end", or no phenomenom was detected for the given
            inputs
        """

        # Termination condition, if time_step (accuracy) is below 1 second
        if time_step >= s_jd:
            # Calculate half time step before
            iteration_between = first_appereance - time_step / 2

            # Condition for start of phenomenom
            if appereance_type == "start":
                # Whether phenomenom was already True a half time step ago
                if Detection.check_phenomena(iteration_between)[row][col].phenomenom_occurs is True:
                    # Find exact Epoch between iteration_between and iteration_between - time_step / 2
                    return Calculation.find_start_end(iteration_between, time_step / 2, row, col, appereance_type)
                else:
                    # Find exact Epoch between first_appereance and first_appereance - time_step / 2
                    return Calculation.find_start_end(first_appereance, time_step / 2, row, col, appereance_type)
            # Condition for end of phenomenom
            elif appereance_type == "end":
                # Whether phenomenom was False a half time step ago
                if Detection.check_phenomena(iteration_between)[row][col].phenomenom_occurs is False:
                    # Find exact Epoch between iteration_between and iteration_between - time_step / 2
                    return Calculation.find_start_end(iteration_between, time_step / 2, row, col, appereance_type)
                else:
                    # Find exact Epoch between first_appereance and first_appereance - time_step / 2
                    return Calculation.find_start_end(first_appereance, time_step / 2, row, col, appereance_type)
            # Input string invalid
            else:
                raise ValueError("Input string invalid")
        # Termination condition met
        else:
            # Check if there really is a phenomenom
            if Detection.check_phenomena(first_appereance - time_step)[row][col].phenomenom_occurs is \
                    Detection.check_phenomena(first_appereance + time_step)[row][col].phenomenom_occurs:
                raise ValueError("No phenomenom detected")
            return first_appereance

    def _cut_to_timespan(self, timing_list: List[TimingResult]) -> List[TimingResult]:
        """Private method to remove content from a given list that is beyond
        the given calculation timespan

        :param timing_list: List that shall be trimmed
        :type timing_list: List[TimingResult]

        :returns: Trimmed list
        :rtype: List[TimingResult]
        """

        # List of indices that will be removed from the list
        to_pop: List[int] = []

        for i in range(len(timing_list)):
            # Find timings later / before calculation timespan
            if timing_list[i].epoch > self.end_epoch or timing_list[i].epoch < self.start_epoch:
                to_pop.append(i)

        # Remove those timings
        for i in sorted(to_pop, reverse=True):
            timing_list.pop(i)

        return timing_list


if __name__ == "__main__":
    epoch_start = Epoch()
    epoch_start.set(2020, 1, 1, 0)

    epoch_stop = Epoch()
    epoch_stop.set(2021, 1, 1, 0)

    calc_time_step = 60 * 120 * 1.157401129603386e-05

    calc = Calculation(epoch_start, epoch_stop, calc_time_step, 0)

    print()

    for tmng in calc.all_timings_sorted:
        print(tmng.epoch.get_full_date(), tmng.row, tmng.phenomenom.phenomenom_type, tmng.appereance_type,
              tmng.phenomenom.shadow_type)

    print()
