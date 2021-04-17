from typing import TextIO, List, Any

from Calendar.Calculation import Calculation
from pymeeus_optimized.Epoch import Epoch


class CreateList(object):
    def __init__(self, year: int):
        self.print_list(year)

    @staticmethod
    def print_file_header(f: TextIO) -> None:
        """ includes all required packages and opens tex file
            :rtype: None
            :param f: tex file into which is written
            :type f: TextIO
            """
        f.write("\\documentclass[12pt, a4paper]{article}\n")
        f.write(
            "\\usepackage[lmargin={2cm},rmargin={2cm},tmargin={0.5cm},bmargin={3.05cm},head=21.75pt]{geometry}\n")
        f.write("\\usepackage{setspace, multicol, graphicx}\n")
        f.write("\\usepackage[footsepline]{scrlayer-scrpage}\n")
        f.write("\\usepackage[english]{babel}\n")
        # f.write("\\usepackage[utf8]{inputenc}\n")
        f.write("\\linespread{0.93}\n")
        # f.write("\\pagestyle{scrheadings}\n")

        f.write("\\begin{document}\n")
        f.write("\\centering\n")

    @staticmethod
    def print_page_header(f: TextIO, epoch: Epoch) -> None:
        """ prints the title on top of the page
        :param f: tex file into which is written
        :type f: TextIO

        :param epoch: epoch in which the events on the page are taking place
        :type epoch: Epoch
        """

        year, month, day = epoch.get_date()

        f.write("\\textbf{" + str(year) + " - PHENOMENA OF JUPITER'S GALILEAN SATELLITES}\\\\(Terrestrial Time) \n")
        f.write("\\vspace{0.1cm} \\hrule \\vspace{0.1cm}\n")

        switcher = {
            1: "JANUARY",
            2: "FEBRUARY",
            3: "MARCH",
            4: "APRIL",
            5: "MAY",
            6: "JUNE",
            7: "JULY",
            8: "AUGUST",
            9: "SEPTEMBER",
            10: "OCTOBER",
            11: "NOVEMBER",
            12: "DECEMBER"
        }
        f.write(switcher.get(month))
        f.write("\\vspace{0.1cm}\n")
        f.write("\\hrule\n\\vspace{-0.2cm}\n")

    @classmethod
    def print_event(cls, f: TextIO, epoch1: Epoch, satellite: int, type1: str, type2: str, type3: str,
                    row_counter: int, column_counter: int, previous_day: int = 0) -> tuple:
        """ prints a single event to a tex file, tex file has to be opened and a tabular has to be already started
                :param f: tex file into which is written
                :type f: TextIO

                :param epoch1: epoch in which the event takes place
                :type epoch1: Epoch

                :param satellite:
                    0: Io
                    1: Europa
                    2: Ganymed
                    3: Callisto
                :type satellite: int

                :param type1:
                    EC: eclipse
                    OC: occultation
                    PA: passgge
                    OM: shadow passage
                :type type1: str

                :param type2:
                    D: start
                    F: finish
                :type type2: str

                :param type3:
                    PEN:
                    EXT:
                    INT:
                :type type3: str

                :param row_counter: number of rows from the beginning of the tabular
                :type row_counter: int

                :param column_counter: number of the tabular on the page
                :type column_counter: int

                :param previous_day: calendar day of the last event
                :type previous_day: int
                """

        year, month, day, hour, minute, second = epoch1.get_full_date()
        second = int(second)

        if row_counter == 74 or row_counter == 75:
            f.write("\t \t \\end{tabular}\n \t}\n")
            row_counter = 0

        if column_counter == 3 and row_counter == 0:
            f.write("\\end{multicols}\n")
            f.write("\\pagebreak\n")
            cls.print_page_header(f, epoch1)
            f.write("\\begin{multicols}{3}\n")
            column_counter = 0

        if row_counter == 0:  # print header
            f.write("\t \\resizebox{0.325\\textwidth}{!}{\n")
            f.write("\t \t \\begin{tabular}{c c c c c c}\n")
            f.write(
                "\t \t \t \\textbf{day} & \\textbf{h} & \\textbf{m} & \\textbf{s} & \\textbf{SAT.} & \\textbf{TYPE}\\\\\n")
            column_counter += 1

        row_counter += 1

        if previous_day != day:
            if row_counter != 1:
                f.write("\t \t \t \t & & & & & \\\\")
                f.write("%row_counter =" + str(row_counter) + " , column_counter =" + str(column_counter) + "\n")
                row_counter += 1
            f.write("\t \t \t \t" + str(day) + " & ")
        else:
            f.write("\t \t \t \t & ")

        f.write(str(hour) + " & " + str(minute) + " & " + str(second) + " & ")
        switcher = {
            0: "I",
            1: "II",
            2: "III",
            3: "IV",
        }
        f.write(switcher.get(satellite) + " & ")
        f.write(type1 + "." + type2 + "." + type3 + "\\\\")
        f.write("%row_counter =" + str(row_counter) + " , column_counter =" + str(column_counter) + "\n")
        previous_day = day
        return row_counter, column_counter, previous_day

    def print_list(self, year) -> None:
        """
        :param year: Year for which the list is created
        :type year: int

        :rtype: None
        """
        f: TextIO = open("../Calendar/"+str(year) + ".tex", "w")
        self.print_file_header(f)
        epoch_start = Epoch()
        epoch_stop = Epoch()
        # 1 s in jd = 1.157401129603386e-05
        calc_time_step = 60 * 120 * 1.157401129603386e-05  # 2h

        row_counter = 0
        column_counter = 0
        previous_day: int = 0
        previous_month: int = 0

        epoch_start.set(year, 1, 1, 0)
        epoch_stop.set(year + 1, 1, 1, 0)

        calc = Calculation(epoch_start, epoch_stop, calc_time_step, 0.0)

        x = len(calc.all_timings_sorted)
        for pos in range(0, x):
            epoch = calc.all_timings_sorted[pos].epoch
            year, month, day = epoch.get_date()
            row = calc.all_timings_sorted[pos].row
            type1 = calc.all_timings_sorted[pos].phenomenon.phenomenon_type
            if calc.all_timings_sorted[pos].appearance_type == 'start':
                type2 = "D"
            else:
                type2 = "F"
            type3 = calc.all_timings_sorted[pos].phenomenon.shadow_type
            if previous_month != month:
                previous_month = month
                if row_counter < 76 and month != 1:
                    f.write("\t \t \\end{tabular}\n \t}\n")
                    f.write("\\end{multicols}\n")
                    row_counter = 0
                    column_counter = 0
                self.print_page_header(f, epoch)
                f.write("\\begin{multicols}{3}\n")

            row_counter, column_counter, previous_day = self.print_event(f, epoch, row, type1, type2, type3,
                                                                         row_counter, column_counter, previous_day)
        if row_counter < 76:
            f.write("\t \t \\end{tabular}\n \t}\n")
            f.write("\\end{multicols}\n")
        f.write("\\end{document}\n")


if __name__ == "__main__":
    CreateList(2021)
