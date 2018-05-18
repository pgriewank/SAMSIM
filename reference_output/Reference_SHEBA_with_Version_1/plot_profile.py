#Standard plot script to plot profiles of three chosen variables in accordance with the finite volume nature of SAMSIM. 
#Steps:
#0. Copy script into folder with the output data
#2. Set dx to the right time step used in the model and in the time unit you want (e.g. in days or years
#3. Set timeunit for xlabel to agree with dx
#4. Set output name with format to save (jpg, png, pdf, or eps). Pdf and eps look better but can be very large and slow to load.
#5. Do you want the freeboard to by incorporated into your y-axis? set free_flag
#6. Set all 3 variable names, unit, and files
#7. Set the output times at which the profiles should be plotted in time (as many as you wish) 
#6. Run the script! (python plot_profile.py)
#7. Open the outputfile in the imageviewer of your choice.

#Loading modules and setting fonts


import numpy
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib

#rc('text', usetex=True)
#rc('font', size='10')
#rc('font', family='serif')

#Settings 
dx           = 3.5 
timeunit     = 'days'

#select output times for profiles 
time      = ([70,80,90,100,110,120,130])

var1name     = 'Bulk salinity'
var1unit     = '[g/kg]'
var1         = numpy.loadtxt("./dat_S_bu.dat")

var2name     = 'Temperature'
var2unit     = '[C]'
var2         = numpy.loadtxt("./dat_T.dat")

var3name     = 'Solid fraction'
var3unit     = '[fraction]'
var3         = numpy.loadtxt("./dat_psi_s.dat")

outputfile   = 'pic_profile'
outputformat = 'pdf' #e.g. png, jpg, pdf
free_flag    = 1     #1: freeboard is included, 0:freeboard is not included 


tlen      = len(time)





#Loading data

thick      = numpy.loadtxt("./dat_thick.dat")
freeboard  = numpy.loadtxt("./dat_freeboard.dat")

#Setting freeboard to zero if free_flag = 0
if   free_flag == 0:
  freeboard[:] = 0.

ylen = len(thick[0,:])
xlen = len(thick[:,0])


#Restructuring the data so it represents the finite volume nature of SAMSIM
depth_step=numpy.hstack((thick,thick))
var1_step =numpy.hstack((thick,thick))
var2_step =numpy.hstack((thick,thick))
var3_step =numpy.hstack((thick,thick))

ylen = len(thick[0,:])
xlen = len(thick[:,0])
i=0
j=0
while (i<xlen):
  while (j<ylen):
    depth_step[i,2*j]	  = -sum(thick[i,0:j])+freeboard[i]
    depth_step[i,2*j+1]	= -sum(thick[i,0:j])-thick[i,j]+freeboard[i]
    var1_step[i,2*j]	   = var1[i,j]
    var1_step[i,2*j+1]	 = var1[i,j]
    var2_step[i,2*j]	   = var2[i,j]
    var2_step[i,2*j+1]	 = var2[i,j]
    var3_step[i,2*j]	   = var3[i,j]
    var3_step[i,2*j+1]	 = var3[i,j]
    j=j+1
  i=i+1
  j=0





#######################################################
#Plotting
#######################################################

fig1=plt.figure(figsize=(9.,5.))
fsize=10.
lsize=9. 


#Set distances between subplots and margins
fig1.subplots_adjust(hspace=0.05)
fig1.subplots_adjust(left=0.07)
fig1.subplots_adjust(right=0.95)
fig1.subplots_adjust(top=0.90)
fig1.subplots_adjust(bottom=0.08)

#######################################################
#Plotting temperature
#######################################################
ax1 = fig1.add_subplot(131)
i=0
while (i<tlen):
  plt.plot(var1_step[time[i],:],depth_step[time[i],:])
  i = i+1
plt.ylabel(r'depth [m]')
plt.xlabel(r''+var1name+' '+var1unit)
#ax1.set_title('Temperature, liquid volume fraction, and bulk salinity')
#ax1.xaxis.set_ticklabels([])

#######################################################
#Plotting liquid fraction
#######################################################
ax2 = fig1.add_subplot(132)
i=0
while (i<tlen):
  plt.plot(var2_step[time[i],:],depth_step[time[i],:],label=str(time[i]*dx))
  i = i+1
plt.xlabel(r''+var2name+' '+var2unit)
ax2.yaxis.set_ticklabels([])

#Plotting legend
handles, labels = ax2.get_legend_handles_labels()
leg = plt.legend(handles,labels,ncol=tlen, mode='expand',handletextpad=1. , title = 'time in '+timeunit,bbox_to_anchor=(-1.0,1.02,3.,0.15),loc=3,borderaxespad=0.,handlelength=2.,prop={'size':lsize},numpoints = 1)
leg.get_frame().set_linewidth(0)

#######################################################
#Plotting salinity
#######################################################
ax3 = fig1.add_subplot(133)
i=0
while (i<tlen):
  plt.plot(var3_step[time[i],:],depth_step[time[i],:])
  i = i+1
plt.xlabel(r''+var3name+' '+var3unit)
ax3.yaxis.set_ticklabels([])


#Saving and exporting
plt.savefig('./'+outputfile+'.'+outputformat,dpi=1000)
#plt.close()

