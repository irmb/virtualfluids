from pyfluids.geometry import *

point1 = GbPoint3D(1, 2, 3)
point2 = GbPoint3D(4, 5, 6)
line = GbLine3D()
cube: GbCuboid3D = GbCuboid3D()

line.point1 = point1
line.point2 = point2

print(point1)
print(line)
print("Distance", point1.get_distance(point2))
