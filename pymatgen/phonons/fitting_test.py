__author__ = 'maarten'

import numpy as np

# try fitting if there are NaN's appearing

arr1 = [-0.0099500000000000144, 0.0050124999999998643, -0.0049874999999999781, 0.010050000000000003]
arr2 = [26.144552839999999, -12.63525993, 12.811005740000001, -25.07676682]

p1 = np.polyfit(arr1, arr2, 1)
#print p1

#print np.sort(arr1)


# Imagine we don't have sigma corresponding to -0.004987...... strain.
# Let's set it to a very high value: 999999

#print arr1
arr2[1] = 99999

print arr1
print arr2
