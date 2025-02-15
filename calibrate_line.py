#  Routine to calibrate the image frames of Halpha for estimation of
#  absolute flux of net emission line by using line and continuum image
#  frames of standard objects and program objects.

import os.path
import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from astropy.io import fits
import time
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.visualization import make_lupton_rgb

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

def obsWavelength(x, vel):
    ''' Corrects for systemic velocity.
    Input: Rest wavelength and velocity in km/s.
    Output: Observed wavelengh.
    '''
    c = 299792.458
    z = vel / c
    lambda_obs = x * (1 + z)
    return lambda_obs

def blackbody(wavelength, T):
    ''' Planck's law, which describes the electromagnetic radiation
    emitted by a black body in thermal equilibrium at a given temperature T.
    (Atomic units, and wavelength in Angstroms)
    '''
    h = 4.136e-15 # eV/s
    kb = 8.617e-5 # eV/K
    c = 2.99792458e18 # A/s
    B = 2 * h * c**2 / (wavelength**5 * (np.exp(h*c / (wavelength*kb*T)) - 1))
    return B

def heatmap(u,g,i,w,x,y):
    ''' Fits the Planck's law to each pixel of a given image.
    Input: flux in u, g and i bands, with their effective wavelengths. Dimension of the image (x,y).
    Output: array with the temperatures given by the fit.
    '''
    start = time.time()
    f = np.empty(shape=(x,y), dtype=float)
    for i in range(0,x):
        for j in range(0,y):
            fluxes = np.array([u_s[i][j], g_s[i][j], i_s[i][j]])
            if (fluxes != 0).all():
                popt, pcov = curve_fit(blackbody, wavelength, fluxes, bounds=(0,1e5))
                T = popt
            else:
                T = 0
            f[i][j] = T
    end = time.time()
    print(end - start,'s')
    return f

if __name__ == "__main__":
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

    # Line filter
    path = '/Users/ziliotto/Documents/GitHub/calibrate_line'
    line = readfile(os.path.join(path, '6568.dat'))
    print('Line filter file:', line[5])
    linex = line[0]
    liney = line[1]
    xlmax = line[3]
    xlmin = line[2]

    # Continuum filter
    continuum = readfile(os.path.join(path, '6092.dat'))
    print('Continuum filter file:', continuum[5])
    xcmax = continuum[3]
    xcmin = continuum[2]
    contx = continuum[0]
    conty = continuum[1]

    # Standard
    standard = readfile(os.path.join(path, 'f34f.dat'))
    print('Standard flux file:', standard[5],'\n')
    standx = standard[0]
    standy = standard[1]
    xmin = standard[2]
    xmax = standard[3]

    # Interpolations
    lineInterp = interp(linex,liney)
    standInterp = interp(standx,standy)
    contInterp = interp(contx,conty)

    # Continuum function (for the case of constant continuum. Use blackbody function in the variable case)
    def continuum():
        ''' Returns the flux of the continuum.
        '''
        funcQ1 = lambda x: standInterp(x) * contInterp(x)
        funcQ2 = lambda x: contInterp(x)
        fc = (1 / fsc) * 10 ** ( 0.4 * kpc * ( Xc - Xsc ) ) * ( integrate.quad(funcQ1, xcmin, xcmax)[0] / integrate.quad(funcQ2, xcmin, xcmax)[0] )
        return fc

    # A
    funcA = lambda x: continuum() * lineInterp(x)
    A = integrate.quad( funcA, xlmin, xlmax, epsabs=1.49e-11 )[0]
    print('A =', A)

    # P
    funcP = lambda x: standInterp(x) * lineInterp(x) # phi_l
    P = (1 / fsl) * 10 ** ( 0.4 * kpl * ( Xl - Xsl ) ) * integrate.quad( funcP, xlmin, xlmax, epsabs=1.49e-11 )[0]
    print('P =', P)

    # Q
    funcQ1 = lambda x: standInterp(x) * contInterp(x) # phi_c
    funcQ2 = lambda x: continuum() * contInterp(x) # phi_c
    Q = (1 / fsc) * 10 ** ( 0.4 * kpc * ( Xc - Xsc ) ) * ( integrate.quad(funcQ1, xcmin, xcmax, epsabs=1.49e-11)[0] / integrate.quad(funcQ2, xcmin, xcmax, epsabs=1.49e-11)[0] )
    print('Q =', Q)

    # R
    fwhm = 1
    g = fwhm / (2 * np.sqrt(np.log(2)))
    cte = 1. / ( np.sqrt( np.pi ) * g )
    waveobs = obsWavelength(6563., vsys)
    funcR = lambda x: lineInterp(x) * np.exp(-(( x - waveobs )/g )**2)
    intR = integrate.quad(funcR, xlmin, xlmax)
    R = intR[0] * cte
    print('R =', R)

    Fl = print('F(line) =', ( 1 / R ), '* (', A, '- ( (I(line) - ', skyl,') / (I(continuum) - ', skyc, ') ) * ( ',texpc*P/(texpl*Q), ') ) ')
    print(Fl)

    u_image_file = fits.open('/Users/ziliotto/Downloads/stacks/ss_fornax_tile1_u_long_ALIGNi.003.fits')
    g_image_file = fits.open('/Users/ziliotto/Downloads/stacks/ss_fornax_tile1_g_long_ALIGNi.003.fits')
    #r_image_file = fits.open('/Users/ziliotto/Downloads/stacks/')
    i_image_file = fits.open('/Users/ziliotto/Downloads/stacks/ss_fornax_tile1_i_long.003.fits')

    u_gain = u_image_file['PRIMARY'].header['GAIN']
    g_gain = g_image_file['PRIMARY'].header['GAIN']
    i_gain = i_image_file['PRIMARY'].header['GAIN']

    u_image_data = fits.getdata('/Users/ziliotto/Downloads/stacks/ss_fornax_tile1_u_long_ALIGNi.003.fits', ext=0)
    g_image_data = fits.getdata('/Users/ziliotto/Downloads/stacks/ss_fornax_tile1_g_long_ALIGNi.003.fits', ext=0)
    i_image_data = fits.getdata('/Users/ziliotto/Downloads/stacks/ss_fornax_tile1_i_long.003.fits', ext=0)

    flux_u = u_image_data / u_gain
    flux_g = g_image_data / g_gain
    flux_i = i_image_data / i_gain

    # Effective wavelengths for u, g, i bands (Angstrom)
    # Available here: http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=CTIO&gname2=DECam&asttype=
    u_w = 3580.78
    g_w = 4773.99
    r_w = 6444.80
    i_w = 7858.77

    wavelength = np.array([u_w,g_w,i_w])

    # Example for an individual galaxy (time of execution ~170s)
    u_s = flux_u[12220:12380,8840:9000]
    g_s = flux_g[12220:12380,8840:9000]
    i_s = flux_i[12220:12380,8840:9000]
    #print(u_s.shape)

    #image = make_lupton_rgb(i_s, g_s, u_s, Q=10, stretch=0.5)
    #plt.grid(None)
    #plt.imshow(image)
    f = heatmap(u_s,g_s,i_s,wavelength,u_s.shape[0],u_s.shape[1])
    plt.imshow(f,cmap='gray')
    plt.grid(None)
    plt.colorbar(label='Temperature (K)')
    plt.show()
