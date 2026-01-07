import numpy as np
import matplotlib.pyplot as plt

def centerline_point_finder(theta1,theta2,mu1,mu2,x1,y1): #arg: theta1 value at a point through which the right runnin characteristics goes, theta2 - value at a point on a centerline
    neg_char_slope=np.tan(0.5*(theta1+theta2)-0.5*(mu1+mu2))
    neg_char_constant=y1-neg_char_slope*x1
    x=-neg_char_constant/neg_char_slope
    y=0
    return x,y

def flow_point_finder(theta1,theta2,theta3,mu1,mu2,mu3,x1,y1,x2,y2): #physical location point finder, arg: theta value at point 1, theta value at point 2, theta value at point 3, mu value at point 1, mu value at point 2, mu value at point 3, location of point 1, location of point 2
    neg_char_slope=np.tan(0.5*(theta1+theta3)-0.5*(mu1+mu3))
    pos_char_slope=np.tan(0.5*(theta2+theta3)+0.5*(mu2+mu3))
    neg_char_constant=y1-neg_char_slope*x1
    pos_char_constant=y2-pos_char_slope*x2
    x=(pos_char_constant-neg_char_constant)/(neg_char_slope-pos_char_slope)
    y=pos_char_slope*x+pos_char_constant
    return x,y

def node_volume_finder(n_characteristics):
    n_nodes_buffer=int(0)
    for i in range(1,n_characteristics):
        n_nodes_buffer+=int(n_characteristics-i)
    output=int(n_nodes_buffer+n_characteristics)
    return output

def Mach_angle_calculator(M):
    return np.arcsin(1/M)

def PM_fun(M,gamma):
    return np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M**2-1)))-np.arctan(np.sqrt(M**2-1))

def Hall_IPM_fun(nu,gamma):
    half_life=np.sqrt((gamma-1)/(gamma+1))
    nu_inf=np.pi/2*(1/half_life-1)
    y=(nu/nu_inf)**(2/3)
    K=(4/(3*np.pi))*(1+1/half_life)
    eta_inf=((3*np.pi)/(2*half_life*(1+half_life)))**(2/3)
    a1=0.5*eta_inf
    a2=((3+8*half_life**2)/40)*eta_inf**2
    a3=((-1+328*half_life**2+104*half_life**4)/2800)*eta_inf**3
    d1=a1-1-(a3-K)/(a2-K)
    d2=a2-a1-((a1-1)*(a3-K))/(a2-K)
    d3=((a3-K)*(a1-K))/(a2-K)-a2+K
    e1=-1-(a3-K)/(a2-K)
    e2=(a3-K)/(a2-K)
    M=(1+d1*y+d2*y**2+d3*y**3)/(1+e1*y+e2*y**2)
    return M

def kernel_node_value_calculator(n_char,theta_start,gamma,Mach_exit): #function for calculationg the characteristic values in the expansion fan of the nozzle
    #initial values calculation
    n_nodes=node_volume_finder(n_char) #number of nodes in the characteristic mesh
    theta_end=PM_fun(Mach_exit,gamma)/2 #throat exit flow halfangle in radians
    theta_step=(theta_end-theta_start)/(n_char-1) #step size for theta values
    #memory allocation for nodes
    K_plus_values=np.zeros(n_nodes)  #Right characteristic values at each node
    K_minus_values=np.zeros(n_nodes) #Left characteristic values at each node
    theta_values=np.zeros(n_nodes) #flow angle values at each node
    nu_values=np.zeros(n_nodes) #Prandtl-Meyer function values at each node
    #initial values on the characteristic lines
    theta_values[0:n_char]=nu_values[0:n_char]=np.arange(theta_start,theta_end+theta_step/2,theta_step) #initial theta and nu values along the characteristic lines
    K_plus_values[0:n_char]=theta_values[0:n_char]-nu_values[0:n_char] #initial K+ values along the characteristic lines
    K_minus_values[0:n_char]=theta_values[0:n_char]+nu_values[0:n_char] #initial K- values along the characteristic lines

    #values for nodes on the centerline
    for i in range(0,int(n_char-1)):
        diff=int(i*(i+3)/2+1)
        theta_values[n_nodes-diff]=0
        K_minus_values[n_nodes-diff]=K_minus_values[n_char-1-i]
        K_plus_values[n_nodes-diff]=2*theta_values[n_nodes-diff]-K_minus_values[n_nodes-diff]
        nu_values[n_nodes-diff]=0.5*(K_minus_values[n_nodes-diff]-K_plus_values[n_nodes-diff])
    #values for internal nodes
    for i in range (3,n_char+1):
        j=int(1)
        k_idx_vec=[]
        k_idx_vec.append(j)
        while j<i-1:
            i_idx=i-1
            for k in k_idx_vec:
                buf_test=n_char-k
                w_idx=i_idx+buf_test
                K_minus_values[w_idx]=K_minus_values[i_idx]
                i_idx=w_idx
            j+=1
            k_idx_vec.append(j)
        k_idx_vec.clear()
    i_idx=n_nodes-3
    for i in range (2,n_char):
        j=int(1)
        k_idx_vec=[]
        k_idx_vec.append(j)
        while j<i:
            for k in k_idx_vec:
                K_plus_values[i_idx+k]=K_plus_values[i_idx]
            j+=1
            k_idx_vec.append(j)
        i_idx=i_idx-i-1
    k_idx_vec.clear()
    for i in range(0,len(nu_values)-1):
        if nu_values[i]==0.:
            theta_values[i]=0.5*(K_minus_values[i]+K_plus_values[i])
            nu_values[i]=0.5*(K_minus_values[i]-K_plus_values[i])
    #Mach number and angle calculation
    M_values=Hall_IPM_fun(nu_values,gamma) #Mach numbers at each node
    mu_values=np.asin(1/M_values)
    #output creation
    output=np.column_stack((np.rad2deg(K_minus_values),np.rad2deg(K_plus_values),np.rad2deg(theta_values),np.rad2deg(nu_values),np.rad2deg(mu_values),M_values))
    return output
    

n_char=int(7) #number of characteristic lines to be drawn
n_nodes=node_volume_finder(n_char) #number of nodes in the characteristic mesh
x_0,y_0=0,1 #location of the first node
M_e=2.4 #exit Mach number
gamma=1.4 #specific heat ratio for exhaust gases
theta_s=np.deg2rad(0.375) # first characteristic line flow angle in radians
theta_e=PM_fun(M_e,gamma)/2 #exit flow angle in radians4
init_theta_vect=np.linspace(theta_s,theta_e,n_char) #initial theta values for the characteristic lines

test=kernel_node_value_calculator(n_char,theta_s,gamma,M_e)

print(test)