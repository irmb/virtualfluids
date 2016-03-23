import math

u = 0.0154946
nu = 7.05882e-05 * pow(2,4) 
y0 = 0.00595276
deltax = 0.078125

q = y0/deltax

utau = math.sqrt(nu*(u/q))

yplus = (utau*q)/nu

print 'u_tau = ', utau
print 'y_plus = ', yplus

