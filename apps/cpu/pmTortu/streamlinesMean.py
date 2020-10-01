import os
import scipy
import numpy
import math

#li=[[5,0.005]]

n=409;
i=0;
j=0;
k=0;
length=[0];
lengthx=[0];
lengthT=[0,0,0,0,0];
lengthxT=[0,0,0,0,0];
ItortuList=[0];
outKOx=78.0;#0.105;
try:
 for i in range(1,n):
  length.append(0);
  lengthx.append(0);
  #if (i<10):
  # dateiIn="C:/Users/Sonja/Documents/pathlines3/pathline30"+str(i)+".dat.ascii.vtu"
  #else:
  dateiIn="C:/Users/Sonja/Documents/blech2c/pathline"+str(i)+".dat.ascii.vtu"
  print dateiIn
  datei = open(dateiIn,"r")
  j=0; k=0;
  for line in datei:
     j=j+1
     #print line
 ##   if ((i>6568) and (i<17261)  ):
     #zuordnung = line.split("  ")
     #if (j==1): print line
     if (k==1): ##marker der anfang der koordinaten bezeichnet
         zuordnung = line.split("  ")
         #print zuordnung
         #pointsxyz=zuordnung[7].split(" ")
         #print pointsxyz
         t=0;
         for entry in zuordnung:
             pointsxyz=entry.split(" ")
             #print pointsxyz
             t=t+1;
             #if (i==1 | i==2):
             lengthT.append(0);
             lengthxT.append(0);
             if (t>7):
              if (pointsxyz[1]!="\n"):
               if(float(pointsxyz[1])<outKOx): ##ende messbereich
                 #print pointsxyz
                 if (t==8):
                     xalt=float(pointsxyz[1]);
                     yalt=float(pointsxyz[2]);
                     zalt=float(pointsxyz[3]);
                 xneu=float(pointsxyz[1]);
                 yneu=float(pointsxyz[2]);
                 zneu=float(pointsxyz[3]);
                 if (xalt>20.0):              ##beginn messbereicht
                  length[i]=length[i]+math.sqrt((xneu-xalt)*(xneu-xalt)+(yneu-yalt)*(yneu-yalt)+(zneu-zalt)*(zneu-zalt));
                  lengthx[i]=lengthx[i]+(xneu-xalt);
                  lengthT[t]=lengthT[t]+length[i];
                  lengthxT[t]=lengthxT[t]+lengthx[i];
                  #print lengthT[t]
                  #print xneu
                  #print lengthx[i]
                 xalt=xneu; yalt=yneu; zalt=zneu;
                 
 
         k=2;
     #if (str(line)=="""            <DataArray type="Float64" NumberOfComponents="3" format="ascii">"""):
 
     if(j==5):
       print line
       k=1;
     #print zuordnung
     #print zuordnung[0]
     #print zuordnung[1]
     #test0=float(zuordnung[0])
     #test1=float(zuordnung[1])
     #li.append([test0,test1])
 ##print float(li[10])/20
 #print li
  datei.close();
 i=0;
 j=0;
 length.pop(0);
 lengthx.pop(0);
 #print length
 #print lengthx
 tortuGes=0;
 LGes=0.0;
 LxGes=0.0;
 fFile = open("f:/temp/pathlinesb2cLength.dat", "w")
 for entry in length:
     #print entry;
     #print lengthx[i];
     LGes=LGes+length[i];
     LxGes=LxGes+lengthx[i];
     ItortuList.append(entry/max(lengthx[i],0.00000001));
     if (length[i]>2.0):
      Itortu=entry/lengthx[i]
      print Itortu
      j=j+1;
      tortuGes=tortuGes+Itortu;
     i=i+1
     fFile.write(str(i))
     fFile.write(" ")
     fFile.write(str(entry))
     fFile.write(" ")
     #fFile.write(str(lengthx[i]))
     #fFile.write(" ")
     #fFile.write(str(entry/max(lengthx[i],0.00000001)))   
     fFile.write("\n")
 tortuGes=tortuGes/j;
 print "berücksichtigte Stromlinien:" 
 print j
 fFile.close();
 ItortuList.pop(0);
 print "TortuGes:"
 print tortuGes;
 print "Lges:"
 print LGes;
 print "Lxges:"
 print LxGes;
 print "Lges/LxGes:"
 print LGes/LxGes;
 erg=[lengthx,length,ItortuList];
 #print erg
 erg=numpy.asarray(erg).T.tolist() #does a list-transpose
 #print erg
 erg=sorted(erg);
 fFile = open("f:/temp/pathlinesb2cLengthSortt1000.dat", "w")
 i=0;
 #print erg[1][1]
 #print erg[0][1]
 for entry in erg:
     i=i+1;
     #print i
     fFile.write(str(entry[0]))
     fFile.write(" ")
     fFile.write(str(entry[1]))
     fFile.write(" ")
     fFile.write(" ")
     fFile.write(str(entry[2]))
     fFile.write("\n")
 fFile.close();
 fFile = open("f:/temp/pathlinesbcbwithTime.dat", "w")
 i=0;
 for entry in lengthxT:
     i=i+1;
     #print i
     fFile.write(str(entry))
     fFile.write(" ")
     fFile.write(str(lengthT[i]))
     fFile.write(" ")
     fFile.write("\n")
 fFile.close();
except IOError:
 datei.close()
 print "caught error couldnt process datafile"
 print i
