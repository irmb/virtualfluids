#script create data file for surface averaged velocity plot
 
from paraview.simple import *

file='/hpc3lustre/work/sonjaOutputs/BKanaltest0Ref2nue0312/steps/stepAV_150000.pvtu'
print file
#read data set file
reader = XMLPartitionedUnstructuredGridReader(FileName=file)

#create slice filter for z
sliceFilter = Slice(reader)
sliceFilter.SliceType.Normal = [0, 0, 1]

#create calculator filter for calculate surface of elements
calc = Calculator(sliceFilter)
calc.ResultArrayName = "one"
calc.Function = "1"

#create integrate variables filter for integration of velocity (integration = summation * surface of elements)
intVal = IntegrateVariables(calc)

f = open('/work/koskuche/SFB880/BKanalPlotData/BKanalAvVx400s2.dat', 'w')   
f.write("z")
f.write(" ")
f.write("avVx")
f.write(" ")
f.write("avVxx")
f.write(" ")
f.write("avVyy")
f.write(" ")
f.write("avVzz")
f.write(" ")
f.write("avVxz")
f.write("\n")

#for i in range (1, 400, 1):
z = 0.995
while z <= 396:
    #shift a slice
    sliceFilter.SliceType.Origin = [298.668884754181, 198.835281848907, z]
    sliceFilter.UpdatePipeline()
    print z
    numPoints=sliceFilter.GetDataInformation().GetNumberOfPoints()
    print numPoints
    intVal.UpdatePipeline()
    pointData=intVal.GetPointDataInformation()
    vxArray=pointData.GetArray("vx")
    vxxArray=pointData.GetArray("vxx")
    vyyArray=pointData.GetArray("vyy")
    vzzArray=pointData.GetArray("vzz")
    vxzArray=pointData.GetArray("vxz")
    sArray  =pointData.GetArray("one")
    surface =sArray.GetRange()[0]
    if numPoints > 0:
      avVx=vxArray.GetRange()[0]/surface
      avVxx=vxxArray.GetRange()[0]/surface
      avVyy=vyyArray.GetRange()[0]/surface
      avVzz=vzzArray.GetRange()[0]/surface
      avVxz=vxzArray.GetRange()[0]/surface

      f.write(str(z))
      f.write(" ")
      f.write(str(avVx))
      f.write(" ")
      f.write(str(avVxx))
      f.write(" ")
      f.write(str(avVyy))
      f.write(" ")
      f.write(str(avVzz))
      f.write(" ")
      f.write(str(avVxz))
      f.write("\n")

    if z < 280:
      z=z+1.99005
    else:
      z=z+3.9801

f.close()

print "!!!ready!!!"