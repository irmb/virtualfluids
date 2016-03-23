print 30*27
print 675/27
print 3*1*9

bMaxX1 = 10
bMaxX2 = 10
bMaxX3 = 10

minX1 = 0
minX2 = 0
minX3 = 0
maxX1 = bMaxX1 - 1
maxX2 = bMaxX2 - 1
maxX3 = bMaxX3 - 1

lMin1X1 = maxX1 - 6
lMax1X1 = lMin1X1 + 4
lMin1X3 = maxX3 - 6
lMax1X3 = lMin1X3

lMin2X1 = maxX1 - 6
lMax2X1 = lMin2X1
lMin2X3 = maxX3 - 6
lMax2X3 = lMin2X3 + 4

lMin1X2 = minX2;
lMax1X2 = maxX2 + - 1;

print 'FC send: ' + str(((lMax1X1-lMin1X1+1)*(lMax1X2-lMin1X2+1)*(lMax1X3-lMin1X3+1) + (lMax2X1-lMin2X1+1)*(lMax1X2-lMin1X2+1)*(lMax2X3-lMin2X3+1)) / 2)

lMin1X1 = maxX1 - 3
lMax1X1 = lMin1X1 + 2
lMin1X3 = maxX3 - 3
lMax1X3 = lMin1X3

lMin2X1 = maxX1 - 3
lMax2X1 = lMin2X1
lMin2X3 = maxX3 - 3
lMax2X3 = lMin2X3 + 2

lMin1X2 = minX2
lMax1X2 = 4

print 'CF receive: ' + str(((lMax1X1-lMin1X1+1)*(lMax1X2-lMin1X2+1)*(lMax1X3-lMin1X3+1) + (lMax2X1-lMin2X1+1)*(lMax1X2-lMin1X2+1)*(lMax2X3-lMin2X3+1)))
