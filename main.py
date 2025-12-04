import numpy as np
import matplotlib.pyplot as plt

t_area=5 #throat area inv m^2
e_area=20 #exit area in m^2
eps=e_area/t_area #expanision area ratio
i_theta=np.deg2rad(23)#typical solid rocket inflation angle as per Sutton
e_theta=np.deg2rad(12)#typical solid rocket exit angle as per Sutton
r_t=np.sqrt(t_area/np.pi) #throat radius in m
i_xlim=0.175348*r_t #limit for the inflation contour
i_rad=0.4*r_t #radius of the inflation contour as per Sutton/Rao
i_center_y=1.4*r_t #center of the inflation contour
cc_rad=12.5 #combustion chamber radius in m
x_c=np.linspace(-1.5*r_t,0,100)
y_c=2.5*r_t-np.sqrt((1.5*r_t)**2-x_c**2)
x_i=np.linspace(0,i_xlim,100) #x values for inflation contour
y_i=i_center_y-np.sqrt(i_rad**2-x_i**2) #inflation contour
plt.plot(x_i,y_i)
plt.plot(x_c,y_c)
plt.gca().set_aspect('equal', 'box')
plt.show()
print("Half circle has been plotted")