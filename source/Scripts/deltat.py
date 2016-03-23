u_real = 17 
u_LB = 0.1 
deltax =  0.078125  
number_of_timesteps = 84000

deltat=(u_LB/u_real)*deltax

print 'delta t =', deltat, 's'
print 'simulated time =', deltat*number_of_timesteps, 'sec'
print 'simulated time =', deltat*number_of_timesteps/60.0, 'min'
