#  Routine to calibrate the image frames of Halpha for estimation of
#  absolute flux of net emission line by using line and continuum image
#  frames of standard objects and program objects.

import os.path
import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline

def readfile(filename):
    ''' Reads a file.
    Input: Filename.
    Returns: x, y, xmin, xmax, number of points, filename.
    '''
    if os.path.isfile(filename) != True:
        print("File not in this directory!")
    else:
        file = np.loadtxt(filename, unpack=True)
        x = file[0]
        y = file[1]
        xmax = np.amax(x)
        xmin = np.amin(x)
        npts = len(file[0])
        return(x, y, xmin, xmax, npts, filename)

def rWavelength(x, vel):
    ''' Corrects for systemic velocity.
    Input: Observed wavelength (number) and velocity in km/s.
    Output: Corrected wavelengh (number).
    '''
    c = 299792458
    z = vel / c
    newW = x / (1 + z)
    return newW

def interp(x, y):
    ''' Interpolates a function fitting a spline y = spl(x).
    Input: x, y (arrays).
    Output: Interpolated function.
    '''
    f = InterpolatedUnivariateSpline(x, y)
    return(f)

def integral(func, a, b):
    ''' Method of InterpolatedUnivariateSpline to calculate the integral of the interpolated function.
    Input: Function to be integrated and interval.
    Output: Result of the integral.
    '''
    I = func.integral(a, b)
    return I

# Data input
Xl = 1.38
print('Line filter air masses (program objects):', Xl)

Xsl = 1.50
print('Line filter air masses (standard objects):', Xsl)

Xc = 1.34
print('Air masses for the continuum filter (program):', Xc)

Xsc = 1.54
print('Air masses for the continuum filter (standard):', Xsc)

kpl = 0.14
print('mag/air mass for line:', kpl)

kpc = 0.23
print('mag/air mass for continuum:', kpc)

fsl = 294.3
print('Raw fluxes of the standard with sky subtracted (line):', fsl)

fsc = 2565.5
print('Raw fluxes of the standard with sky subtracted (continuum):', fsc)

nline = 1
print('Number of lines found in the line filter range, its rest wavelengths and fractional contribution (sum=1.0):', nline)

rfWave = [6563,1]
print('Rest wavelength (Angstroms) of each line as well as its fractional contibution (sum = 1.0):', rfWave)

vsys = 45
print('Systemic velocity of the galaxy:', vsys)

texpl = 2400
print('Exposure times of program frames (line):', texpl)

texpc = 600
print('Exposure times of program frames (continuum):', texpc)

skyl = 67.9
print('Sky background of the program frames in counts/pixel (line):', skyl)

skyc = 123.8
print('Sky background of the program frames in counts/pixel (continuum):', skyc)

# Line file
line = readfile('/Users/tuilaziliotto/Documents/GitHub/calibrate_line/6568.dat')
print('Line filter file:', line[5])
linex = line[0]
liney = line[1]
xlmax = line[3]
xlmin = line[2]

# Continuum file
continuum = readfile('/Users/tuilaziliotto/Documents/GitHub/calibrate_line/6092.dat')
print('Continuum filter file:', continuum[5])
xcmax = continuum[3]
xcmin = continuum[2]
contx = continuum[0]
conty = continuum[1]

# Standard file
standard = readfile('/Users/tuilaziliotto/Documents/GitHub/calibrate_line/f34f.dat')
print('Standard flux file:', standard[5],'\n')
standx = standard[0]
standy = standard[1]
xmin = standard[2]
xmax = standard[3]

# Interpolations
lineInterp = interp(linex,liney)
standInterp = interp(standx,standy)
contInterp = interp(contx,conty)

# Plotting the data and interpolations
#plt.figure()
#plt.subplot(121)
#linexNew = np.arange(xlmin,xlmax,0.1)
#lineyNew = lineInterp(linexNew)
#plt.plot(linex, liney, 'o')
#plt.plot(linexNew, lineyNew)
#plt.title('Interpolation for line filter')
#plt.ylabel('Transmission (%)')
#plt.xlabel(r'Wavelength ($\AA$)')
#plt.subplot(122)
#contxNew = np.arange(xcmin,xcmax,0.1)
#contyNew = contInterp(contxNew)
#plt.plot(contx, conty, 'o')
#plt.plot(contxNew, contyNew)
#plt.title('Interpolation for continuum filter')
#plt.ylabel('Transmission (%)')
#plt.xlabel(r'Wavelength ($\AA$)')
#plt.subplots_adjust(top = 0.7,bottom = 0.3,left = 0.10,hspace = 0.9,wspace = 0.5)
#plt.show()
#standxNew = np.arange(xmin,xmax,0.1)
#standyNew = standInterp(standxNew)
#plt.plot(standx, standy, 'o')
#plt.plot(standxNew, standyNew)
#plt.title('Interpolation for standard')
#plt.ylabel(r'Flux ($erg/s/cm^2/\AA$)')
#plt.xlabel(r'Wavelength ($\AA$)')
#plt.show()

# Continuum filter integration
contInt = integral(contInterp,xcmin,xcmax)

# Line filter integration
lineInt = integral(lineInterp,xlmin,xlmax)

# Standard flux file integration
standInt = integral(standInterp,xmin,xmax)

# Integration of continuum filter * standard flux
auxContFunc = lambda x: contInterp(x) * standInterp(x)
contFluxInt = integrate.quad(auxContFunc,xcmin,xcmax,epsabs=1.49e-11)

# Integration of line filter * standard flux
auxLineFunc = lambda x: lineInterp(x) * standInterp(x)
lineFluxInt = integrate.quad(auxLineFunc,xlmin,xlmax,epsabs=1.49e-11)

# Q
factorQ = 10 ** ( 0.4 * ( kpc * ( Xc - Xsc ) ) )
Q = (factorQ / fsc) * (contFluxInt[0] / 1) * (lineInt / contInt) # changed
# lineFluxInt[0] to 1
print('Q =', Q)

# P
factorP = 10 ** ( 0.4 * ( kpl * ( Xl - Xsl ) ) )
P = (factorP / fsl) * lineFluxInt[0]
print('P =', P)

# R
fwhm = 1
g = fwhm / (2 * np.sqrt(np.log(2)))
factorR = 1 / ( np.sqrt( np.pi * g ) )
funcR = lambda x: lineInterp(x) * np.exp( - ( ( x - rWavelength(x,vsys) )/g )**2 )
RInt = integrate.romberg(funcR,xlmin,xlmax)
R = factorR * RInt
print('R =', R,'\n')

alpha = P / ( R * texpl )
beta = Q * texpl / ( P * texpc )
gamma = -skyl + ( skyc * texpl * Q  / ( texpc * P ) )

print('F = ', alpha, '(I(line) -', beta, 'I(continuum) +', gamma)
