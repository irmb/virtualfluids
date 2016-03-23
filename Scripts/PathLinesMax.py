#script create data file for maximum strain rate in path lines
 
from paraview.simple import *

file='d:/temp/pathLine/pathLine.pvtu'
print file
#read data set file
reader = XMLPartitionedUnstructuredGridReader(FileName=file)

threshold = Threshold(reader)
threshold.Scalars = 'ID'
r=1.0
threshold.ThresholdRange = [r, r]
print threshold.ThresholdRange
#help(threshold)

pcalc = PythonCalculator(threshold)
pcalc.Expression = 'max(tauxy)'
pcalc.UpdatePipeline()
#pcalc.GetOutput()
#help(pcalc)
#pcalc.inputs[0].Points[:,0]
print pcalc.PointData.Points
points = pcalc.Points
#help(points)

print pcalc.GetPointDataInformation().GetArray("result").GetRange()[0]