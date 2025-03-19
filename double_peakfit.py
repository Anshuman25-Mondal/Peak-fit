
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

os.chdir('C:/Users/mondalan/AnacondaProjects/19062021/CSV1/') #path
myfiles = [f for f in os.listdir() if f.endswith(".xy")]

xdata=[]
ydata=[]

for f in myfiles:
    file=open(f,"r+")
    contents=file.readlines()
    file.close
    contents=contents[23:]
    xdata.append(np.array([float(entry.split()[0]) for entry in contents]))
    ydata.append(np.array([float(entry.split()[1]) for entry in contents]))

starting_file=0
end_file=400 #can be used for handling multiple files together, use '1' for 1 file. 

wavelength=0.4845


twotheta = np.array(xdata[starting_file])
intensity = np.array(ydata[starting_file])
plt.figure(1)
plt.plot(intensity)
pts=np.array(plt.ginput(2))
plt.close(1)

x_exp=twotheta[int(pts[0,0]):int(pts[1,0])]
y_exp=intensity[int(pts[0,0]):int(pts[1,0])]

bk=(np.mean(intensity[int(pts[0,0]):int(pts[0,0]+3)])+np.mean(intensity[int(pts[1,0]-3):int(pts[1,0])]))/2

def pv(x,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10):
    y = a1 + a3*(1 - a2)*(np.sqrt(4*np.log(2))/(np.sqrt(np.pi)*a5))*np.exp(-(4*np.log(2)/a5**2)*(x-a4)**2)+ a3*a2*(2/np.pi)*(a5/(4*(x-a4)**2 + a5**2))+ a6*x + a8*(1 - a7)*(np.sqrt(4*np.log(2))/(np.sqrt(np.pi)*a10))*np.exp(-(4*np.log(2)/a10**2)*(x-a9)**2)+ a8*a7*(2/np.pi)*(a10/(4*(x-a9)**2 + a10**2))
    return y

plt.figure(2)
plt.plot(x_exp,y_exp)
pts2=np.array(plt.ginput(2))
plt.close(2)

a0=[bk,0.5,pts2[0,1],pts2[0,0],0.1,bk,0.5,pts2[1,1],pts2[1,0],0.25]


aopt,acov=curve_fit(pv,x_exp,y_exp,a0)

plt.figure(3)
plt.plot(x_exp,y_exp)
plt.plot(x_exp,pv(x_exp,*aopt))

results_array=np.zeros([len(myfiles),12])
results_array[starting_file]=np.array([*aopt,pts[0,0],pts[1,0]])

a_lattice=(wavelength/(2*np.sin(results_array[starting_file,3]*np.pi/2/180)))*np.sqrt(3)
V_lattice=((a_lattice)**3)/4
V0_au=16.9635;
K0_au=166.65;
Kprime_au=5.4823;
P= (3/2)*K0_au*((V0_au/V_lattice)**(7/3)-(V0_au/V_lattice)**(5/3))*(1+(3/4)*(Kprime_au-4)*((V0_au/V_lattice)**(2/3)-1))

print(P)

for p in range(starting_file+1,end_file):
    twotheta = np.array(xdata[p])
    intensity = np.array(ydata[p])
    
    x_exp = twotheta[int(results_array[p-1,10]):int(results_array[p-1,11])]
    y_exp = intensity[int(results_array[p-1,10]):int(results_array[p-1,11])]
    
    bk=(np.mean(intensity[int(results_array[p-1,10]):int(results_array[p-1,10]+3)])+np.mean(intensity[int(results_array[p-1,11]-3):int(results_array[p-1,11])]))/2
    a0 = [bk, results_array[p-1,1], results_array[p-1,2],results_array[p-1,3],results_array[p-1,4],bk, results_array[p-1,6],results_array[p-1,7],results_array[p-1,8],results_array[p-1,9]]
    
    aopt,acov=curve_fit(pv,x_exp,y_exp,a0)

    idx_low = (np.abs(twotheta - results_array[p-1,3])).argmin()
    idx_high = (np.abs(twotheta - aopt[3])).argmin()
    
    pts_0_moving = results_array[p-1,10]+idx_high-idx_low
    pts_1_moving = results_array[p-1,11]+idx_high-idx_low
    
    results_array[p]=np.array([*aopt,pts_0_moving,pts_1_moving])
    
plt.figure(4)
plt.plot(x_exp,y_exp)
plt.plot(x_exp,pv(x_exp,*aopt))
    
del a_lattice, V_lattice, P
a_lattice=wavelength/(2*np.sin(results_array[:end_file,3]*np.pi/2/180))*np.sqrt(3)
V_lattice=(a_lattice)**3/4
P= (3/2)*K0_au*((V0_au/V_lattice)**(7/3)-(V0_au/V_lattice)**(5/3))*(1+(3/4)*(Kprime_au-4)*((V0_au/V_lattice)**(2/3)-1))

plt.figure(5)
plt.plot(P)
