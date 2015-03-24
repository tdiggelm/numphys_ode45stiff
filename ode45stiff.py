import matplotlib.pyplot as plt
from numpy import diff, finfo, double
from ode45 import ode45
from numpy import *
from time import time
from pylab import *
from scipy.optimize import fsolve

rhs = lambda y: 500 * y**2 * (1-y)
fun = lambda t, y: rhs(y)

def integrate_EE(y0, xStart, xEnd, steps, flag=False):
    r"""Integrate ODE with explicit Euler method
    Input: y0     ... initial condition
           xStart ... start x
           xEnd   ... end   x
           steps  ... number of steps (h = (xEnd - xStart)/N)
           flag   ... flag == False return complete solution: (phi, phi', t)
                      flag == True  return solution at endtime only: phi(tEnd)
    Output: x ... variable
            y ... solution
    """
    x = zeros(steps)
    y = zeros((steps, size(y0)))
    ###########################################################
    #                                                         #
    # TODO: Implementieren Sie hier die explizite Euler Regel #
    #       zur integration der funktion y(x).                #
    #                                                         #
    ###########################################################

    h = double(xEnd)/steps
    y[0,:] = y0

    for k in xrange(steps-1):
        y[k+1] = y[k] + h * rhs(y[k])
        x[k+1] = (k+1)*h

    if flag:
        return x[-1], y[-1][:]
    else:
        return x, y


def integrate_IE(y0, xStart, xEnd, steps, flag=False):
    r"""Integrate ODE with implicit Euler method
    Input: y0     ... initial condition
           xStart ... start x
           xEnd   ... end   x
           steps  ... number of steps (h = (xEnd - xStart)/N)
           flag   ... flag == False return complete solution: (phi, phi', t)
                      flag == True  return solution at endtime only: phi(tEnd)
    Output: x ... variable
            y ... solution
    """
    x = zeros(steps)
    y = zeros((steps, size(y0)))
    ###########################################################
    #                                                         #
    # TODO: Implementieren Sie hier die implizite Euler Regel #
    #       zur integration der funktion y(x).                #
    #                                                         #
    ###########################################################
    h = double(xEnd)/steps
    y[0,:] = y0

    for k in xrange(steps-1):
        F = lambda x: x - y[k] - h * rhs(x)
        y[k+1] = fsolve(F, y[k] + h * rhs(y[k]))
        x[k+1] = (k+1)*h

    if flag:
        return x[-1], y[-1][:]
    else:
        return x, y


def integrate_IM(y0, xStart, xEnd, steps, flag=False):
    r"""Integrate ODE with implicit midpoint rule
    Input: y0     ... initial condition
           xStart ... start x
           xEnd   ... end   x
           steps  ... number of steps (h = (xEnd - xStart)/N)
           flag   ... flag == False return complete solution: (phi, phi', t)
                      flag == True  return solution at endtime only: phi(tEnd)
    Output: x ... variable
            y ... solution
    """
    x = zeros(steps)
    y = zeros((steps, size(y0)))
    #################################################################
    #                                                               #
    # TODO: Implementieren Sie hier die implizite Mittelpunktsregel #
    #       zur integration der funktion y(x).                      #
    #                                                               #
    #################################################################
    h = double(xEnd)/steps
    y[0,:] = y0

    for k in xrange(steps-1):
        F = lambda x: x - y[k] - h * rhs(0.5*(x + y[k]))
        y[k+1] = fsolve(F, y[k] + h * rhs(y[k]))
        x[k+1] = (k+1)*h

    if flag:
        return x[-1], y[-1][:]
    else:
        return x, y


y0 = array([0.01])

# Compute
ts = 0
te = 1
nrsteps = 200

starttime = time()
t_ee, y_ee = integrate_EE(y0, ts, te, nrsteps, False)
endtime = time()
print('EE needed %f seconds' % (endtime-starttime))

starttime = time()
t_ee, y_ie = integrate_IE(y0, ts, te, nrsteps, False)
endtime = time()
print('IE needed %f seconds' % (endtime-starttime))

starttime = time()
t_ee, y_im = integrate_IM(y0, ts, te, nrsteps, False)
endtime = time()
print('IM needed %f seconds' % (endtime-starttime))

# Plot
fig = figure(figsize=(12,8))
ax = fig.gca()
ax.set_aspect("equal")

ax.scatter(t_ee, y_ee, marker="+", color="red", label="EE", s=30)
ax.scatter(t_ee, y_ie, marker="x", color="green", label="IE", s=30)
ax.scatter(t_ee, y_im, marker="v", color="blue", label="IM", s=30)

ax.grid(True)
ax.legend(loc="lower right")
ax.set_xlabel("x")
ax.set_ylabel("y")
savefig("ode45stiff.pdf")
show()

"""
def ode45stiff():


	# define the differential equation
	fun = lambda t,x: 500*x**2*(1-x)
	# define the time span
	tspan = (0.0,1.0)
	# get total time span
	L = diff(tspan)
	# starting value
	y0 = 0.01

	# make the figure
	fig = plt.figure()

	# exact solution
	# set options
	eps = finfo(double).eps
	options = {'reltol':10*eps,'abstol':eps,'stats':'on'}
	tex,yex = ode45(fun,(0,1),y0,**options)
	# plot 'exact' results
	plt.plot(tex,yex,'--',color=[0,0.75,0],linewidth=2,label='$y(t)$')
	plt.hold(True)

	# ODE45
	# set options
	options = {'reltol':0.01,'abstol':0.001,'stats':'on'}
	# get the solution using ode45
	[t,y] = ode45(fun,(0,1),y0,**options)
	# plot the solution
	plt.plot(t,y,'r+',label='ode45')

	# set label, legend, ....
	plt.xlabel('t',fontsize=14)
	plt.ylabel('y')
	plt.title('ode45 for d_ty = 500 y^2(1-y)')
	plt.show()
	# write to file
	#	plt.savefig('../PICTURES/logode451.eps')

	# now plot stepsizes
	fig = plt.figure()
	plt.title('ode45 for d_ty = 500 y^2(1-y)')
	ax1 = fig.add_subplot(111)
	ax1.plot(tex,yex[:,2], 'r--', linewidth=2)
	ax1.set_xlabel('t')
	ax1.set_ylabel('y(t)')

	ax2 = ax1.twinx()
	ax2.plot(0.5*(t[:-1]+t[1:]),diff(t), 'r-')
	ax2.set_ylabel('stepsize')

	plt.show()

	# write to file
	#plt.savefig('../PICTURES/logode452.eps')


if __name__ == '__main__':
	print 'warning: this function takes a long time to complete! (not really worth it)\n'
	ode45stiff()
"""
