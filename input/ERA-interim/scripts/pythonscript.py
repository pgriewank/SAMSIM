#Is used to cut ERA-interim data into the thrr hourly values fed into SAMSIM for atmoflux_flag = 2
#Is mostly there to turn the integrated values given my ERA_interim in to per time units
#Only works if the ERA-interim data is already sorted into ascci files called precip.grb.txt, lw.grb.txt, sw.grb.txt, and T2m.grb.txt 


import numpy as n

precip = n.loadtxt("./precip.grb.txt")

lenr=len(precip)


precip2 = precip*10000000
precip2[0:lenr:4]=precip[0:lenr:4]
precip2[1:lenr:4]=precip[1:lenr:4]-precip[0:lenr:4]
precip2[2:lenr:4]=precip[2:lenr:4]-precip[1:lenr:4]
precip2[3:lenr:4]=precip[3:lenr:4]-precip[2:lenr:4]
precip3 = precip2  / (3*3600) /3 


sw = n.loadtxt("./sw.grb.txt")
sw2 = sw*10000000
sw2[0:lenr:4]=sw[0:lenr:4]
sw2[1:lenr:4]=sw[1:lenr:4]-sw[0:lenr:4]
sw2[2:lenr:4]=sw[2:lenr:4]-sw[1:lenr:4]
sw2[3:lenr:4]=sw[3:lenr:4]-sw[2:lenr:4]
sw3 = sw2 / (3*3600)

lw  = n.loadtxt("./lw.grb.txt")
lw  = lw
lw2 = lw*10000000
lw2[0:lenr:4]=lw[0:lenr:4]
lw2[1:lenr:4]=lw[1:lenr:4]-lw[0:lenr:4]
lw2[2:lenr:4]=lw[2:lenr:4]-lw[1:lenr:4]
lw2[3:lenr:4]=lw[3:lenr:4]-lw[2:lenr:4]


lw3=lw2/(3*3600)



T2m = n.loadtxt("./T2m.grb.txt")
T2m2 = T2m - 273.15

#Add possible shortening here, otherwise leave cut=0
cut = 8*365/2

T2mf    = T2m2[cut:lenr]
precipf = precip3[cut:lenr]
lwf     = lw3[cut:lenr]
swf     = sw3[cut:lenr]


#OUTPUT

n.savetxt('T2m.txt.input',T2mf)
n.savetxt('precip.txt.input',precipf)
n.savetxt('flux_sw.txt.input',swf)
n.savetxt('flux_lw.txt.input',lwf)


