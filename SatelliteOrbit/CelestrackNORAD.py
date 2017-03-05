'''
TURKSAT 4A              
1 39522U 14007A   15301.78105273  .00000128  00000-0  00000+0 0  9996
2 39522   0.0299 272.9737 0004735 326.3457 120.6614  1.00271335  6265
'''

import math
import sys

import astropy.units as u


name = sys.stdin.readline().strip()
line1 = sys.stdin.readline()
line2 = sys.stdin.readline()
number = int(line1[2:7]) # + line1[7]
year = int(line1[9:11]) + 1000
if year < 1057: # 2k
    year = year + 1000
launch = line1[11:14]
piece = line1[14]
epoch = line1[18:32]

i = float(line2[8:16])
raan = float(line2[17:25])
e = float("0." + line2[26:33].strip())
ap = float(line2[34:42])
ma = float(line2[43:51])
f = float(line2[52:63])
revs = int(line2[63:68])

t = 1.0 / f * u.day

print ap
print t
print t.to(u.second)
print 1.0 / t.to(u.second)
