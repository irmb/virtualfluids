import os
os.system("date")

pathin = "/gfs1/work/niivfcpu/scratch/kucher/BKanaltest/steps/step*"
pathout = "/gfs1/work/niivfcpu/scratch/kucher/BKanaltest/tar/step"
tar = "tar -cvf "

start = 281750
step = 2250
stop = start + step * 10

print stop

for i in range(start, stop, step):
  command = tar + pathout + str(i) + ".tar " + pathin + str(i) + "*"
  #print command
  os.system(command)