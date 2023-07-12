import math



location = [1, 4, 12]
recipBoxSize = [0.5, 0.5, 0.5]
periodicBoxVectors = [
                [2, 0, 0],
                [3, 2, 0],
                [0, 1, 2]
]

scale2 = math.floor(location[2]*recipBoxSize[2])
yperiodic = location[1]-periodicBoxVectors[2][1]*scale2
print(scale2)
zperiodic = location[2]-periodicBoxVectors[2][2]*scale2
scale1 = math.floor(yperiodic*recipBoxSize[1])
yperiodic -= periodicBoxVectors[1][1]*scale1
print(scale1, periodicBoxVectors[1][0])

print(yperiodic, zperiodic)