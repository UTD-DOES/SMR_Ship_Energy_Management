import os,sys
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd



# Path stuff
sys.path.append(r"C:\Program Files\PTI\PSSE35\35.3\PSSPY37")
sys.path.append(r"C:\Program Files\PTI\PSSE35\35.3\PSSBIN")
sys.path.append(r"C:\Program Files\PTI\PSSE35\35.3\PSSLIB")
sys.path.append(r"C:\Program Files\PTI\PSSE35\35.3\EXAMPLE")
os.environ['PATH'] = (r"C:\Program Files\PTI\PSSE35\35.3\PSSPY37;" 
  + r"C:\Program Files\PTI\PSSE35\35.3\PSSBIN;" 
  + r"C:\Program Files\PTI\PSSE35\35.3\EXAMPLE;" + os.environ['PATH'])

 
#pssbindir  = r"C:\Program Files\PTI\PSSE35\PSSBIN"
#os.environ['PATH'] = pssbindir + ';' + os.environ['PATH']
#os.system("pssplt -inpdev test.idv")

import psse35
psse35.set_minor(3)


import psspy
import pssplot

import redirect
redirect.psse2py()

psspy.psseinit(1500)


import dyntools
import matplotlib.pyplot as plt
from psspy import _i, _f, _s, _o


psspy.case(r"""C:--Address File--\Grid_v33_B0.sav""")
#psspy.machine_chng_4(1,r"""1""",[_i,_i,_i,_i,_i,2,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],_s)


psspy.dyre_new([1,1,1,1],r"""C:\--Address File--\Modeling_v2.dyr""","","","")


psspy.change_plmod_con(2,r"""1""",r"""PSS2A""",11, 20.0)
psspy.change_plmod_con(3,r"""1""",r"""PSS2A""",11, 20.0)

psspy.plmod_remove(3,r"""1""",1)
psspy.plmod_remove(3,r"""1""",6)
psspy.plmod_remove(3,r"""1""",7)
psspy.plmod_remove(3,r"""1""",3)


psspy.add_wind_model(1,r"""1""",1,r"""REGCA1""",1,[0],[""],14,[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",1, 0.017)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",2, 10.0)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",3, 0.1)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",4, 0.05)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",5, 1.22)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",6, 1.2)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",7, 0.2)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",8, 0.05)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",9,-1.3)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",10, 0.02)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",12, 99.0)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",13,-99.0)
psspy.change_wnmod_con(1,r"""1""",r"""REGCA1""",14, 0.7)
psspy.add_wind_model(1,r"""1""",2,r"""REECC1""",5,[0,0,0,0,0],["","","","",""],45,[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",1,-99.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",2, 99.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",3, 0.01)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",4,-0.05)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",5, 0.05)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",6, 15.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",7, 0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",8,-0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",9, 1.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",10, 0.05)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",11, 0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",12,-0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",13, 1.1)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",14, 0.9)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",16, 1.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",18, 1.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",19, 0.017)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",20, 99.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",21,-99.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",22, 1.1)       #PMAX (pu), Max. power limit
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",23,-0.567)     #PMIN (pu), Min. power limit
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",24, 1.11)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",25, 0.017)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",27, 0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",28, 0.2)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",29, 0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",30, 0.5)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",31, 0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",32, 1.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",33, 0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",34, 0.2)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",35, 1.11)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",36, 0.5)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",37, 1.11)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",38, 0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",39, 1.11)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",40, 1.0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",41, 1.11)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",42, 999.0)     #T, battery discharge time (s)  (>0)
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",43, 0.5)       #SOCini (pu), Initial state of charge
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",44, 0.9)       #SOCmax (pu), Maximum allowable state of charge
psspy.change_wnmod_con(1,r"""1""",r"""REECC1""",45, 0.1)       #SOCmin (pu), Minimum allowable state of charge
psspy.change_wnmod_icon(1,r"""1""",r"""REECC1""",4,0)
psspy.add_wind_model(1,r"""1""",7,r"""REPCA1""",7,[0,0,0,0,0,0,0],["","","","","","",""],27,[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",1, 0.02)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",2, 0.0008)     #Kp, Reactive power PI control proportional gain (pu)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",3, 0.008)      #Ki, Reactive power PI control integral gain (pu)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",5, 0.05)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",10, 0.1)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",11,-0.1)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",14, 0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",15,-0.75)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",16, 1.0)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",18, 0.25)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",19,-0.00083)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",20, 0.00083)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",21, 99.0)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",22,-99.0)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",23, 1.0)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",24,-0.667)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",25, 0.1)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",26, 126.0)       #Ddn, reciprocal of droop for over-frequency conditions (pu)
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",27, 126.0)       #Dup, reciprocal of droop for under-frequency conditions (pu)
psspy.change_wnmod_icon(1,r"""1""",r"""REPCA1""",6,0)            ##RefFlag (flag for V or Q control): 0: Q control 1: voltage control
psspy.change_wnmod_icon(1,r"""1""",r"""REPCA1""",7,1)            ##Fflag (flag to disable frequency control): 1: Enable control 0: disable
psspy.change_wnmod_con(1,r"""1""",r"""REPCA1""",17, 0.0001)
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,2,0])
psspy.chsb(0,1,[-1,-1,-1,1,3,0])
psspy.chsb(0,1,[-1,-1,-1,1,4,0])
psspy.chsb(0,1,[-1,-1,-1,1,6,0])
psspy.chsb(0,1,[-1,-1,-1,1,25,0])
psspy.chsb(0,1,[-1,-1,-1,1,26,0])
psspy.chsb(0,1,[-1,-1,-1,1,12,0])
psspy.chsb(0,1,[-1,-1,-1,1,13,0])
psspy.chsb(0,1,[-1,-1,-1,1,14,0])
psspy.chsb(0,1,[-1,-1,-1,1,16,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[ 100.0,0.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[ 100.0,0.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[ 100.0,0.0,0.0, 100.0])


###### Reliable Operation and Import Load profile #######

path = r"C:\--Address File--\P_motor6.csv"

df = pd.read_csv(path)

loads = list(df.Load4)
#print (loads)



psspy.strt_2([0,1],r"""C:\--Address File--\Operation_B2025.out""")
psspy.run(0, 0.2,1,1,0)

st = 0.5
time_step = 0.05

##for LD in loads:
##    psspy.load_chng_6(6,r"""1""",[_i,_i,_i,_i,_i,_i,_i],[_f,_f, LD,_f,_f,_f,_f,_f],"")
##    st = st + time_step
##    psspy.run(0, st,1,1,0)



for idx, LD in enumerate(loads):
    # Apply load
    psspy.load_chng_6(6, r"""1""", [_i,_i,_i,_i,_i,_i,_i], [_f,_f, LD,_f,_f,_f,_f,_f], "")

    # Change Pref based on time index
    if 246 <= idx <= 352:
        psspy.change_plmod_con(2, '1', 'IEESGO', 1, 1.2)  # Bus 2, machine ID '1'
    elif 353 <= idx <= 573:
        psspy.change_plmod_con(2, '1', 'IEESGO', 1, 0.9)
    elif 574 <= idx <= 783:
        psspy.change_plmod_con(2, '1', 'IEESGO', 1, 1.0)
    elif 784 <= idx <= 945:
        psspy.change_plmod_con(2, '1', 'IEESGO', 1, 1.1)  # Bus 2, machine ID '1'
    elif 946 <= idx <= 1095:
        psspy.change_plmod_con(2, '1', 'IEESGO', 1, 0.9)
    elif idx >= 1096:
        psspy.change_plmod_con(2, '1', 'IEESGO', 1, 1.1)

    # Run simulation to next time step
    st += time_step
    psspy.run(0, st, 1, 1, 0)



#psspy.run(0, 3.0,1,1,0)

##########################################################



outfile = r"""C:\--Address File--\Operation_B2025.out"""
chnfobj = dyntools.CHNF(outfile)
short_title, chanid, chandata = chnfobj.get_data()


t = chandata['time']
P1 = chandata[5]              #based on the channel identifier set in the dynamic simulation
P2 = chandata[6]              #based on the channel identifier set in the dynamic simulation

freq = chandata[22]              #based on the channel identifier set in the dynamic simulation 27-35
Pload = chandata[17]
Volt = chandata[31]               # 31-40



plt.plot(t, P1, label = "Bat 1")
plt.plot(t, P2, label = "Gen 2")
plt.title('Power of Battery')
plt.xlabel('Time [s]')
plt.ylabel('Power [MW]')
plt.legend()
plt.show()






