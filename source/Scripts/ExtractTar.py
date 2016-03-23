import os
os.system("date")

pathin = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/tar/step"
pathout = "/hpc3lustre/work/koskuche/SFB880/BKanalTestBerlin/steps"
tar = "tar -xvf "

start = 491000
step = 2250
stop = start + step * 2

print stop

for i in range(start, stop, step):
  command = tar + pathin + str(i) + ".tar -C " + pathout + " --strip-components=7"
  #print command
  os.system(command)