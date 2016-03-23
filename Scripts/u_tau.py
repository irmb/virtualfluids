import math

k = 0.41
u1=1.2
u2=0
y1=6
y2=0.05

u_tau = k/2.3 * (u2 - u1) / math.log10(y1/y2)

print 'u_tau = ', u_tau 
