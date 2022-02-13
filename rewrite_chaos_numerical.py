#20210220 Chua's circuit
#論文：ROBUST OP AMP REALIZATION OF CHUA'S CIRCUIT
#使用1.9.12.13頁的資料來模擬

#論文：利用混沌系統偵測植物覺知人體之反應
#使用12頁的資料來模擬


#------------------------------------------------------------------------------
# import time
# time_start = time.time() #開始計時，了解程式跑多久
#------------------------------------------------------------------------------

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def pend(y, t, C1, C2, L, R, R1, R2, R3, R4, R5, R6, Bp):
# y is a vecter[V1, V2, I]
# V1對應論文中的v_C1
# V2對應論文中的v_C2
# I對應論文中的i_L
    V1, V2, I = y


# Differential equations (Michael Peter Kennedy(1994))
# Create dydt = (V1',V2',I'):
    dydt = [
         ((V2-V1)/R-V1/R1-(-1/R3-1/R1)*(abs(V1+R3/(R2+R3)*Bp)-abs(V1-R3/(R2+R3)*Bp))/2
         -V1/R4-(-1/R6-1/R4)*(abs(V1+R6/(R5+R6)*Bp)-abs(V1-R6/(R5+R6)*Bp))/2)/C1,

         ((V1-V2)/R+I)/C2,

         -(V2)/L]
    return dydt



# Parameter values (使用楊,林(2017) 第12頁的的元件)
C1=10*10**-9
C2=100*10**-9
L=18*10**-3
R1=220
R2=220
R3=2.2*10**3
R4=22*10**3
R5=22*10**3
R6=3.3*10**3
Bp=8.3

R=1800.0  #variable resistor


# Initial conditions
# V10 and V20 are the initial voltage; I0 are the initial current 
# y0 = [V10, V20, I0]
y0 = [0.0, 2.0, 0.0]


# ODE solver parameters
abserr = 1.0e-8 #error
relerr = 1.0e-6 #error
endtime=0.2 #stop time
numpoints = 100000 #number of points


# Create the time samples for the output of the ODE solver.
t = np.linspace(0,endtime,numpoints+1)


# Call the ODE solver.
sol = odeint(pend, y0, t, args=(C1, C2, L, R, R1, R2, R3, R4, R5, R6, Bp),atol=abserr, rtol=relerr)
#use args to pass the parameters to pend


#------------------------------------------------------------------------------
# save data (可以不用)
# with open('rewrite_chaos.dat', 'w') as f:
#     # Print & save the solution.
#     for t1, s in zip(t, sol):
#         print(t1, s[0], s[1], s[2], file = f)
#------------------------------------------------------------------------------


# Plot the solution that was generated
plt.plot(sol[10000:len(sol),0],sol[10000:len(sol),1],lw=0.3) #delete the first 10000 records
#plt.plot(t[10000:len(sol)],sol[10000:len(sol),1],lw=0.3)
plt.title('chaos')
plt.xlabel(r'$V_1$ (V)')
plt.ylabel(r'$V_2$ (V)')
plt.grid(True)
plt.savefig('rewrite_chaos.png', dpi=500)


#------------------------------------------------------------------------------
# 計時 (可以不用)
# time_end = time.time()    #結束計時
# time_c= time_end - time_start   #執行所花時間
# print('time cost', time_c, 's')
#------------------------------------------------------------------------------

plt.show()