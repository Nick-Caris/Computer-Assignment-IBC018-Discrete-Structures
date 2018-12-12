import csv
import sys
import numpy as np

'''
   importFunction is a function where you can put in the filename and you get the function and the intialCondtions as output.

   Input:
   Filename     The location and name of the file you want to convert to an np array

   Output:
   function             The function i.e. the line where it say's s(n)
   initialConditions    The initialConditions i.e. the line where it say's s(int)
 '''


def importFunction(filename):
    function = np.empty(1, dtype=str)
    initialConditions = []

    with open(filename, 'r') as f:
        reader = csv.reader(f)
        try:
            for row in reader:
                if row[0][0] == 's':
                    if row[0][2] == 'n':
                        function = row[0]
                    else:
                        initialConditions.append(row)
        except csv.Error as e:
            sys.exit('file %s, line %d: %s' % (filename, reader.line_num, e))
        initialConditions = np.array(initialConditions)
        return function, initialConditions


# Hier een teste om te checken of die werkt
# een, twee = importFunction('../testData/comass03.txt')
#
# print('--------------')
# print(een)
# print('---------------')
# print(twee)
