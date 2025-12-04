import numpy as np
import matplotlib.pyplot as plt

t_area=5 #throat area inv m^2
e_area=20 #exit area in m^2
eps=e_area/t_area #expanision area ratio
i_theta=np.deg2rad(23)#typical solid rocket inflation angle as per Sutton
e_theta=np.deg2rad(12)#typical solid rocket exit angle as per Sutto
r_t=np.sqrt(t_area/np.pi) #throat radius in m
xlim=0.175348*r_t #limit for the inflation contour
x=np.linspace(0,xlim,100)
y=1.4*r_t-np.sqrt((0.4*r_t)**2-x**2)
plt.plot(x,y)
plt.gca().set_aspect('equal', 'box')
plt.show()
print(y)
print("Half circle has been plotted")