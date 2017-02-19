import numpy as np
from numpy import column_stack
from matplotlib import pyplot as plt 
#from matplotlib.pyplot import *
from matplotlib import animation
from math import pi, log
from scipy import integrate, interpolate
import os

# paths
workdir  = os.getcwd()
workdir  = workdir + '/..'
outdir   = workdir + '/../output/go/Ca1em01/Bo1ep00/Ma0ep00/dr5em02/dt1em05'
#outdir   = workdir + '/../output/go/Ca1em02/Bo1ep00/Ma0ep00/dr5em02/dt1em06'
#outdir   = workdir + '/../output/go/Ca1em03/Bo1ep00/Ma0ep00/dr5em02/dt1em07'
#outdir   = workdir + '/../output/go/Ca1em02/Bo1ep00/Ma0ep00/dr5em02/dt1em05'

outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma0ep00/tstop1-2/dr5em02/dt1em05'
outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma1ep00/tstop1-2/dr5em02/dt1em05'
outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma1ep01/tstop1-2/dr5em02/dt1em05'
outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma1ep02/tstop1-2/dr5em02/dt1em05'
outdir   = workdir + '/../output/stop/Ca2em02/Bo1ep00/Ma0ep00/tstop1-3/dr5em02/dt1em05'
outdir   = workdir + '/../src3/output'
#outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma1ep02/tstop1-2/dr5em02/dt1em05'
#outdir   = workdir + '/../output/stop/Ca1em01/Bo1ep00/Ma1ep02/tstop1-2/dr5em02/dt1em4'

# filenames
tfile    = '/t.txt'
rfile    = '/r.txt'
hfile    = '/h.txt'
gfile    = '/g.txt'
ffile    = '/f.txt'
pfile    = '/p.txt'
qfile    = '/q.txt'
vfile    = '/vs.txt'

# load data
tdata    = np.loadtxt(outdir + tfile, unpack=True,skiprows=0)
rdata    = np.loadtxt(outdir + rfile, unpack=True,skiprows=0)
hdata    = np.loadtxt(outdir + hfile, unpack=True,skiprows=0)
gdata    = np.loadtxt(outdir + gfile, unpack=True,skiprows=0)
fdata    = np.loadtxt(outdir + ffile, unpack=True,skiprows=0)
pdata    = np.loadtxt(outdir + pfile, unpack=True,skiprows=0)
qdata    = np.loadtxt(outdir + qfile, unpack=True,skiprows=0)
vdata    = np.loadtxt(outdir + vfile, unpack=True,skiprows=0)

# number of space and time points
J1 = len(tdata[:,0])
N1 = len(tdata[0,:])
J  = J1 - 1
N  = N1 - 1

# set up the figure, the axis, and the plot element to be animated
fig = plt.figure(figsize=(15,8))
ax1 = fig.add_subplot(111)
line1a, = ax1.plot([], [], 'b-')
line1b, = ax1.plot([], [], 'k-')
time1   = ax1.text(0.85, 0.05, '', transform=ax1.transAxes)

ax1.set_xlim( 0, 4)
ax1.set_ylim(-1.5, 1.5)
ax1.set_xlabel('$\sigma / \sqrt{a b}$',fontsize=18)
ax1.set_ylabel('$h/b$',fontsize=18)

# initialization function: plot the background of each frame
def init() :
	line1a.set_data([], [])
	line1b.set_data([], [])
	time1 .set_text('')
	return line1a, line1b, time1

# animation function, to be called sequentially
def animate(i) :
	t   =  tdata[0,i]
	x   =  rdata[:,i]
	y1a =  hdata[:,i]
	y1b = -fdata[:,i]

	line1a.set_data(x, y1a)
	line1b.set_data(x, y1b)

	time1 .set_text("$t U/b$ = " + str(t))

	# shift axes on ax1
	shift = True
	if (shift) :
		# time range
		i0 = 40
		t0 = 4.0
		t1 = 20.0
		dt = t1 - t0

		# h plot
		xmin0, xmax0 = ax1.get_xlim()
		ymin0, ymax0 = ax1.get_ylim()
		xmax1 =  0.8
		dxmax = (xmax1 - xmax0)/dt
		ymin1 =  0.0
		ymax1 =  0.4
		dymin = (ymin1 - ymin0)/dt
		dymax = (ymax1 - ymax0)/dt
		if i > i0 and xmax0 > xmax1 :
			xmax = xmax0 + (i - i0)*dxmax
			ymin = ymin0 + (i - i0)*dymin
			ymax = ymax0 + (i - i0)*dymax
			ax1.set_xlim( xmin0, xmax)
			ax1.set_ylim( ymin , ymax)

	print "t = " + str(t)
	return line1a, line1b, time1

# call the animator. blit=True means only re-draw the parts that have changed
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=N1, interval=45, repeat=True, blit=False)

# format and save movie
#fileName = 'hg_go_Ca1em01_Bo1ep00_Ma0ep00_dr5em02_dt1em05.mp4'
#fileName = 'hg_go_Ca1em02_Bo1ep00_Ma0ep00_dr5em02_dt1em06.mp4'
#fileName = 'hg_go_Ca1em03_Bo1ep00_Ma0ep00_dr5em02_dt1em07.mp4'
#fileName = 'hg_stop_Ca1em03_Bo1ep00_Ma0ep00_tstop_1-2_dr5em02_dt1em05.mp4'
#fileName = 'hg_stop_Ca1em03_Bo1ep00_Ma1ep00_tstop_1-2_dr5em02_dt1em05.mp4'
#fileName = 'hg_stop_Ca1em03_Bo1ep00_Ma1ep01_tstop_1-2_dr5em02_dt1em05.mp4'
fileName = 'h_vs_t_Ca0-1_Bo1_d1-3_Ma0.mp4'
Writer = animation.writers['ffmpeg']
writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)
anim.save(fileName, writer=writer)
#writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)
#anim.save('im.mp4', writer=writer)


plt.show()
