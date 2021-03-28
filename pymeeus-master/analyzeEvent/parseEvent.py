import sqlite3
from datetime import datetime
import copy


class ParseEvent(object):

    def __init__(self, file: str, origin: str):

        self.file = file
        self.origin = origin

        # parser state
        self.state = 'searchYear'

        # define 300x18 list
        self.cleanTable = [[0, 0, 0, 0, 0, '', 0, 0, 0, 0, 0, '', 0, 0, 0, 0, 0, ''] for i in range(300)]
        self.table = copy.deepcopy(self.cleanTable)

        # current Date and Time
        self.dateTimeAdded = format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

        # SQL database
        self.database = "Event.sqlite"

        # open SQL database
        self.conn = sqlite3.connect(self.database)
        self.cursor = self.conn.cursor()

        # call parser (i.e. start processing)
        self.parse()

        # close SQL database
        self.conn.commit()
        self.conn.close()

    def writeTable(self, rowTotal: int):

        for column in range(3):
            for row in range(rowTotal):
                elements = self.table[row][6 * column:6 * (column + 1)]

                if elements[4] > 0:
                    # row contains event data for a satellite

                    date = datetime(self.year, self.month, elements[0], elements[1], elements[2], elements[3]).strftime(
                        '%Y-%m-%d %H:%M:%S')
                    types = elements[5].split('.')

                    add = ("INSERT INTO Event "
                           "(DateTimeEventTT, Satellite, Type1, Type2, Type3, Origin, DateTimeAdded) "
                           "VALUES (?, ?, ?, ?, ?, ?, ?)")
                    data = (date, elements[4], types[0], types[1], types[2], self.origin, self.dateTimeAdded)
                    self.cursor.execute(add, data)

    def monthNo(self, monthText: str):
        monthList = {
            "JANVIER": 1,  # French
            "FEVRIER": 2,
            "MARS": 3,
            "AVRIL": 4,
            "MAI": 5,
            "JUIN": 6,
            "JUILLET": 7,
            "AOUT": 8,
            "SEPTEMBRE": 9,
            "OCTOBRE": 10,
            "NOVEMBRE": 11,
            "DECEMBRE": 12,

            "JANUARY": 1,  # English
            "FEBRUARY": 2,
            "MARCH": 3,
            "APRIL": 4,
            "MAY": 5,
            "JUNE": 6,
            "JULY": 7,
            "AUGUST": 8,
            "SEPTEMBER": 9,
            "OCTOBER": 10,
            "NOVEMBER": 11,
            "DECEMBER": 12
        }
        month = monthList.get(monthText, -1)
        # try to fix not-recognized months due to special french characters
        if month == -1:
            if monthText[:1] == 'F' and monthText[-5:] == 'VRIER':
                month = 2
            if monthText[:2] == 'AO' and monthText[-1:] == 'T':
                month = 8
            if monthText[:1] == 'D' and monthText[-6:] == 'CEMBRE':
                month = 12
        return month

    def satNo(self, satText: str):
        satList = {
            "I": 1,
            "II": 2,
            "III": 3,
            "IV": 4
        }
        return satList.get(satText, -1)

    def fixTable(self, rowTotal: int):

        # detect empty rows in individual columns, if the file is of multi column type
        # do from bottom to top to work with fixed rows, as always row above is fixed based on row below
        # start with row that has three columns, i.e. nothing has to be moved, i.e. infos are in correct column
        if self.multicolumn:
            for row in range(rowTotal - 1, 0, -1):
                # find bottom row that has three columns, i.e. nothing has to be moved, i.e. infos are in correct column
                if self.table[row][16] > 0:
                    bottomThreeColumnRow = row
                    continue

            column = 1
            for row in range(bottomThreeColumnRow, 0, -1):
                if self.table[row][6 * (column - 1)] != 0:
                    # detected new day, i.e. empty column in row above
                    # i.e. move contents to the right
                    self.table[row - 1][2 * 6:3 * 6] = self.table[row - 1][1 * 6:2 * 6]
                    self.table[row - 1][1 * 6:2 * 6] = self.table[row - 1][0 * 6:1 * 6]
                    self.table[row - 1][0:6] = [0, 0, 0, 0, 0, '']

            column = 2
            for row in range(bottomThreeColumnRow, 0, -1):
                if self.table[row][6 * (column - 1)] != 0:
                    # detected new day, i.e. empty column in row above
                    # i.e. move contents to the right
                    self.table[row - 1][2 * 6:3 * 6] = self.table[row - 1][1 * 6:2 * 6]
                    self.table[row - 1][1 * 6:2 * 6] = [0, 0, 0, 0, 0, '']

        # detect empty rows in individual columns, last row, if the file is of multi column type
        if self.multicolumn:
            column = 1
            if self.table[0][6 * column] != 0:
                # detect new day top second column means empty last row first column
                # i.e. move contents to the right
                self.table[rowTotal - 1][2 * 6:3 * 6] = self.table[rowTotal - 1][1 * 6:2 * 6]
                self.table[rowTotal - 1][1 * 6:2 * 6] = self.table[rowTotal - 1][0 * 6:1 * 6]
                self.table[rowTotal - 1][0:6] = [0, 0, 0, 0, 0, '']

            column = 2
            if self.table[0][6 * column] != 0:
                # detect new day top third column means empty last row second column
                # i.e. move contents to the right
                self.table[rowTotal - 1][2 * 6:3 * 6] = self.table[rowTotal - 1][1 * 6:2 * 6]
                self.table[rowTotal - 1][1 * 6:2 * 6] = [0, 0, 0, 0, 0, '']

        # detect empty rows in individual columns, if the file is of multi column type
        # do from bottom to top to work with fixed rows, as always row above is fixed based on row below
        # do "rest", i.e. from bottom up to row with three columns where we started before
        if self.multicolumn:
            column = 1
            for row in range(rowTotal - 1, bottomThreeColumnRow, -1):
                if self.table[row][6 * (column - 1)] != 0:
                    # detected new day, i.e. empty column in row above
                    # i.e. move contents to the right
                    self.table[row - 1][2 * 6:3 * 6] = self.table[row - 1][1 * 6:2 * 6]
                    self.table[row - 1][1 * 6:2 * 6] = self.table[row - 1][0 * 6:1 * 6]
                    self.table[row - 1][0:6] = [0, 0, 0, 0, 0, '']

            column = 2
            for row in range(rowTotal - 1, bottomThreeColumnRow, -1):
                if self.table[row][6 * (column - 1)] != 0:
                    # detected new day, i.e. empty column in row above
                    # i.e. move contents to the right
                    self.table[row - 1][2 * 6:3 * 6] = self.table[row - 1][1 * 6:2 * 6]
                    self.table[row - 1][1 * 6:2 * 6] = [0, 0, 0, 0, 0, '']

        # transfer day from previous page if needed
        if self.table[0][0] == 0:
            self.table[0][0] = self.dayFoundLastPage
        # save self.dayFound from this page for usage in next page
        self.dayFoundLastPage = self.dayFound

        # in first column: transfer day from row above if zero
        column = 1
        for row in range(1, rowTotal):
            if self.table[row][6 * (column - 1)] == 0:
                self.table[row][6 * (column - 1)] = self.table[row - 1][6 * (column - 1)]

        if self.multicolumn:
            # transfer day from first bottom to second top
            if self.table[0][6] == 0:
                self.table[0][6] = self.table[rowTotal - 1][0]

            # in second column: transfer day from row above
            column = 2
            for row in range(1, rowTotal):
                if self.table[row][6 * (column - 1)] == 0:
                    self.table[row][6 * (column - 1)] = self.table[row - 1][6 * (column - 1)]

            # transfer day from second bottom to third top
            if self.table[0][12] == 0:
                self.table[0][12] = self.table[rowTotal - 1][6]

            # in third column: transfer day from row above
            column = 3
            for row in range(1, rowTotal):
                if self.table[row][6 * (column - 1)] == 0:
                    self.table[row][6 * (column - 1)] = self.table[row - 1][6 * (column - 1)]

    def parse(self):

        # initialize row
        row = 0

        with open(self.file) as f:

            for line in f:

                elements = line.split()

                if self.state == 'searchYear' or self.state == 'searchMonth':
                    # could be encountered when searching for the year at the beginning of the file
                    # or when searching for the next month, will evaluate year at each encounter
                    # in case a file contains multiple years

                    if len(elements) > 4:
                        if elements[4] == 'SATELLITES' or elements[2] == 'PHENOMENA':
                            self.year = int(elements[0])
                            print('Found year: ', self.year)
                            self.state = 'searchMonth'
                            continue

                if self.state == 'searchMonth':

                    if len(elements) > 0:
                        if elements[0] == '(Temps' or elements[0] == '(Terrestrial':
                            continue
                        else:
                            self.month = self.monthNo(elements[0])
                            print('Found month:', elements[0])
                            self.state = 'searchHeading'
                            # initialize new table
                            self.table = copy.deepcopy(self.cleanTable)
                            self.multicolumn = False  # initialize as single column format, set to multicolumn once
                            # a multicolumn is found
                            row = 0
                            continue

                if self.state == 'searchHeading' or self.state == 'searchEvents':
                    # heading could be encountered while searching for heading (Thuillot)
                    # or while searching for events if heading is repeated (PyMeeus)

                    if len(elements) > 0:
                        if elements[0] == 'jour' or elements[0] == 'day':
                            self.state = 'searchEvents'
                            continue

                if self.state == 'searchEvents':

                    if len(elements) > 4:
                        if elements[4] == 'SATELLITES' or elements[2] == 'PHENOMENA':
                            # reached start of next page
                            self.state = 'searchMonth'

                            # fix and write table
                            self.fixTable(row)
                            self.writeTable(row)

                            continue

                    if len(elements) >= 15:
                        columns = 3
                        self.multicolumn = True
                    elif len(elements) >= 10:
                        columns = 2
                        self.multicolumn = True
                    elif len(elements) >= 5:
                        columns = 1
                    else:
                        # reached end of table
                        continue

                    for column in range(0, columns):

                        day = 0
                        if elements[4][0] == 'I':
                            # examines first char of sat string, catches all four cases 'I', 'II', 'III', 'IV'
                            # to determine that row contains day-element
                            self.dayFound = int(elements[0])  # store day found in object variable, needed during
                            # fixing table in case of single column table to transfer to next page
                            day = self.dayFound
                            print('Found day: ', day)
                            del elements[:1]  # remove day-element

                        self.table[row][6 * column:6 * (column + 1)] = \
                            (day, int(elements[0]), int(elements[1]), int(elements[2]),
                             self.satNo(elements[3]), elements[4])

                        del elements[:5]  # remove h m s sat type

                    row = row + 1  # gets only executed if a row is processed and populated in table

        if row > 0:
            # fix and write last table
            self.fixTable(row)
            self.writeTable(row)


if __name__ == "__main__":
    # first argument is filename of data to be imported (to be parsed) into SQLite database named Event.sqlite
    # set second argument to "Thuillot" if Thiullot data is imported, set to "PyMeeus" if PyMeeus data is imported
    ParseEvent("jupiterConfigPyMeeus-18.txt", "PyMeeus")
