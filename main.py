import numpy as np
import matplotlib.pyplot as plt
import csv

r_t=5 #throat radius in mm
r_e=23 #exit plane radius in mm
i_theta=np.deg2rad(26)#typical solid rocket inflation angle as per Sutton
e_theta=np.deg2rad(16)#typical solid rocket exit angle as per Sutton

#r_t=np.sqrt(t_area/np.pi) #throat radius in m

i_rad=0.4*r_t #radius of the inflation contour as per Sutton/Rao
i_center_y=r_t+i_rad #center of the inflation contour
i_xlim=np.sqrt(((i_rad)**2*np.tan(i_theta)**2/(1+np.tan(i_theta)**2))) #limit for the inflation contour
i_y_lim=i_center_y-np.sqrt(i_rad**2-i_xlim**2)
x_i=np.linspace(0,i_xlim,100) #x values for inflation contour
y_i=i_center_y-np.sqrt(i_rad**2-x_i**2) #inflation contour

cc_rad=12.5 #combustion chamber radius in m
c_deg=np.deg2rad(-20) #convergence half-angle

c_rad=1.5*r_t #convergence contour radius
c_cir_cent=r_t+1.5*r_t #convergence contour circle centre
c_xlim=-np.sqrt(((1.5*r_t)**2*np.tan(c_deg)**2/(1+np.tan(c_deg)**2)))
x_c=np.linspace(c_xlim,0,150)
y_c=c_cir_cent-np.sqrt((c_rad)**2-x_c**2)

r_mat=np.array([[2*i_y_lim,1,0],[2*r_e,1,0],[i_y_lim**2,i_y_lim,1]]) #radius matrix as per 'Design and analysis of contour bell nozzle and comparison with dual bell nozzle'
b_vec=np.array([[1/np.tan(i_theta)],[1/np.tan(e_theta)],[i_xlim]]) #vector of tangens/x_n as per 'Design and analysis of contour bell nozzle and comparison with dual bell nozzle'
r_mat_inv=np.linalg.inv(r_mat)

const_vec=np.dot(r_mat_inv,b_vec)
y_parabola=np.linspace(i_y_lim,r_e,100)
x_parabola=const_vec[0]*y_parabola**2+const_vec[1]*y_parabola+const_vec[2]

x_full=np.hstack((x_c,x_i,x_parabola)).ravel()
x_unique,idx=np.unique(x_full,return_index=True) #deleting duplicates from vector
y_full=np.hstack((y_c,y_i,y_parabola)).ravel()
y_unique=y_full[idx] #deleting duplicates from vector

with open('Nozzle.csv','w',newline='') as csvfile:
    coordinates_writer=csv.writer(csvfile,delimiter=';')
    for i in range(len(x_unique)):
        coordinates_writer.writerow([x_unique[i],y_unique[i],0])
        


plt.plot(x_unique,y_unique)
plt.gca().set_aspect('equal', 'box')
plt.xlabel('Współrzędna x [mm]')
plt.ylabel('Współrzędna y [mm]')
plt.grid(True,which='major',alpha=0.6)
plt.grid(True,which='minor',alpha=0.3)
plt.minorticks_on()
plt.savefig('Nozzle.png',bbox_inches='tight')
plt.show()