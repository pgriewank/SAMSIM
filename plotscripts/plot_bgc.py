#This is the standard plotscript to print bulk and brine concentration of bgc tracers.
#It plots the brine and bulk concentration of a specific tracer and saves it as plot_bgc+tracernumber+format.

#Steps needed
#0. Copy script into folder with the output data
#1. Set which tracer you want to plot the settings, will be added to outputfilename
#2. Set dx to the right time step used in the model and in the time unit you want (e.g. in days or year
#3. Set timeunit for xlabel to agree with dx
#4. Choose format to save (jpg, png, pdf, or eps). Pdf and eps look better but can be very large and slow to load.
#5. Do you want the freeboard to by incorporated into your y-axis? set free_flag
#6. Run the script! (python plot_bgc.py)
#7. Open the outputfile in the imageviewer of your choice.
#8. Fine tune everything to fit your needs

#Warning: Plotscript will not run if bulk and brine concentration are always zero.


#Settings 
tracernumber = "01"  #Tracernumber of the tracer you want to plot, begins at 01
dx           = 1. 
timeunit     = '[days]'
outputfile   = 'pic_bgc'
outputformat = 'png' #e.g. png, jpg, pdf
free_flag    = 1     #1: freeboard is included, 0:freeboard is not included 

#Contour levels 
levelsbulk   = ([50, 100])
levelsbrine  = ([400,2000])

#Loading modules and setting fonts
import numpy
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', size='10')
rc('font', family='serif')

#Loading data
bgc_br     = numpy.loadtxt("./dat_bgc"+tracernumber+".br.txt")
bgc_bu     = numpy.loadtxt("./dat_bgc"+tracernumber+".bu.txt")
thick      = numpy.loadtxt("./dat_thick.dat")
snow       = numpy.loadtxt("./dat_snow.dat")
freeboard  = numpy.loadtxt("./dat_freeboard.dat")

#Setting freeboard to zero if free_flag = 0
if   free_flag == 0:
  freeboard[:] = 0.
 

#Restructuring the data so it can be ploted by pcolor
depth = thick*1.
depth_contour = thick*1.
Xgrid = thick*1.

ylen = len(thick[0,:])
xlen = len(thick[:,0])
Xaxis = numpy.arange(0,xlen*dx,dx)

i=0
j=0
ireal = 0.
while (i<xlen):
  while (j<ylen):
    depth[i,j]=-sum(thick[i,0:j])+freeboard[i]
    Xgrid[i,j]=ireal
    depth_contour[i,j]=-sum(thick[i,0:j])-thick[i,j]/2.+freeboard[i]
    j=j+1
  i=i+1
  j=0
  ireal=ireal+dx


#######################################################
#Plotting
#######################################################


fig1=plt.figure(figsize=(7.,4.))
fsize=10.
whatisup=freeboard+snow[:,0]
ymin=(depth.min()+freeboard.min())*1.03
ymax=whatisup.max()+depth.min()*-0.03

#Set distances between subplots and margins
fig1.subplots_adjust(hspace=0.08)
fig1.subplots_adjust(left=0.08)
fig1.subplots_adjust(right=1.03)
fig1.subplots_adjust(top=0.95)
fig1.subplots_adjust(bottom=0.1)

#######################################################
#Plotting bulk concentrations
#######################################################
ax1 = fig1.add_subplot(211, axisbg='grey')
zmin=0.
zmax=bgc_bu.max()
plt.pcolor(Xgrid,depth,bgc_bu,cmap='ocean_r',vmin=zmin,vmax=zmax)
c1 = plt.colorbar(pad=0.01)
c1.set_label(r'BULK content /kg')
plt.axis([Xgrid.min(), Xgrid.max(), ymin, ymax])
CS1 = plt.contour(Xgrid, depth_contour, bgc_bu, levelsbulk,colors='k')
plt.clabel(CS1, fontsize=9, inline=1,fmt='%1.0f')
ax1.fill_between(Xaxis[:],freeboard[:], snow[:,0]+freeboard[:,],facecolor='white',edgecolor='white')
plt.plot(Xaxis[:],freeboard[:],'k--')
plt.ylabel(r'depth [m]')
#ax1.set_title('Tracer II')
ax1.xaxis.set_ticklabels([])

#######################################################
#Plotting brine concentrations
#######################################################
ax2 = fig1.add_subplot(212,axisbg='grey')
fsize=10.
zmin=0.
zmax=bgc_br.max()
plt.pcolor(Xgrid,depth,bgc_br,cmap='RdPu',vmin=zmin,vmax=zmax)
plt.axis([Xgrid.min(), Xgrid.max(), ymin, ymax])
CS2 = plt.contour(Xgrid, depth_contour, bgc_br, levelsbrine,colors='k')
plt.clabel(CS2, fontsize=9, inline=1,fmt='%1.0f')
c2 = plt.colorbar(pad=0.01)
c2.set_label(r'BRINE content /kg')
ax2.fill_between(Xaxis[:],freeboard[:], snow[:,0]+freeboard[:,],facecolor='white',edgecolor='white')
plt.plot(Xaxis[:],freeboard[:],'k--')
plt.ylabel(r'depth [m]')
plt.xlabel(r'time '+timeunit)



#Saving and exporting
plt.savefig("./"+outputfile+'_'+tracernumber+"."+outputformat,dpi=1000)
plt.close()


