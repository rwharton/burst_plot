import numpy as np
import matplotlib.pyplot as plt


def read_pulsars(psrfile, freq):
    """
    Read width and rlum for pulsars 
    freq in GHz
    """
    ws = []
    lums = []
    with open(psrfile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", ' ', '\n']:
                continue
            else: pass

            cols = line.split()
            w50 = float(cols[1])
            lum = float(cols[2]) 

            ws.append(w50)
            lums.append(lum)

    ws = np.array(ws)
    lums = np.array(lums) 

    # convert from ms to s
    ws = ws * 1e-3

    # convert from mJy kpc^2 to Jy kpc^2
    lums = lums * 1e-3

    wv = ws * freq

    return wv, lums


def check_rval(val):
    if val.strip() == '--':
        return -999
    else:
        return float(val)


def read_rrats(infile):
    fluxes = []
    freqs  = []
    ds     = []
    widths = []

    with open(infile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", ' ', '\n']:
                continue
            else: pass

            cols = line.split()

            d = check_rval(cols[11])

            if d == -999:
                continue
            else: pass

            s140  = check_rval(cols[12])
            s350  = check_rval(cols[13])
            s1400 = check_rval(cols[14])

            w140  = check_rval(cols[15])
            w350  = check_rval(cols[16])
            w1400 = check_rval(cols[17])

            if (s140 > -999) and (w140 > -999):
                flux = s140
                w    = w140
                freq = 0.140
            
            elif (s350 > -999) and (w350 > -999):
                flux = s350
                w    = w350
                freq = 0.350
            
            elif (s1400 > -999) and (w1400 > -999):
                flux = s1400
                w    = w1400
                freq = 1.400

            else: continue

            ds.append(d)
            fluxes.append(flux)
            widths.append(w)
            freqs.append(freq)

        ds = np.array(ds)
        fluxes = np.array(fluxes)
        widths = np.array(widths)
        freqs = np.array(freqs)

        # Get lums in Jy kpc^2
        lums = (1e-3 * fluxes) * ds**2.0

        # Get freq * width in GHz s
        wv = freqs * 1e-3 * widths

        return wv, lums

def check_fval(val):
    val = val.strip()
    val = val.strip('\"')
    if val.strip() == "null":
        return -999
    else:
        return float(val)

def read_frbs(infile):
    freqs  = []
    widths = []
    fluxes = []
    redshifts = []
    
    with open(infile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", ' ', '\n']:
                continue
            else: pass

            cols = line.split(',')

            freq = check_fval(cols[4])
            width = check_fval(cols[6])
            flux  = check_fval(cols[8])
            z     = check_fval(cols[9])

            if freq <= 0:
                continue
            if width <= 0:
                continue
            if flux <= 0:
                continue

            freqs.append(freq)
            widths.append(width)
            fluxes.append(flux)
            redshifts.append(z)

    freqs = np.array(freqs)
    widths = np.array(widths)
    fluxes = np.array(fluxes)
    redshifts = np.array(redshifts)

    # freqs to GHz
    freqs *= 1e-3

    # width to sec
    widths *= 1e-3

    return freqs, widths, fluxes, redshifts


def read_frb180916(infile, d=1.49e5):
    spks = []
    ws   = []
    freqs = []

    with open(infile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", " ", "\n"]:
                continue
            else: pass

            cols = line.split()
            
            spk = float(cols[0])
            w   = float(cols[1])
            freq = float(cols[2])

            spks.append(spk)
            ws.append(w)
            freqs.append(freq)

    spks = np.array(spks)
    ws = np.array(ws) * 1e-3
    freqs = np.array(freqs)

    wvs = ws * freqs
    lums = spks * d**2.0

    return wvs, lums

def read_frb121102(infile, d=950e3):
    spks = []
    ws   = []
    freqs = []

    with open(infile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", " ", "\n"]:
                continue
            else: pass

            cols = line.split()
            
            spk = float(cols[0])
            w   = float(cols[1])
            freq = float(cols[2])

            spks.append(spk)
            ws.append(w)
            freqs.append(freq)

    spks = np.array(spks)
    ws = np.array(ws) * 1e-3
    freqs = np.array(freqs)

    wvs = ws * freqs
    lums = spks * d**2.0

    return wvs, lums


def read_B0540(infile, d=50):
    """
    Get burst info
    """
    spks = []
    ws   = []
    freqs = []

    with open(infile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", " ", "\n"]:
                continue
            else: pass

            cols = line.split()
            
            freq = float(cols[0])
            w   = float(cols[1])
            spk_obs = float(cols[5])
            fluence = float(cols[6])
            spk_int = fluence / w

            spks.append(spk_obs)
            ws.append(w)
            freqs.append(freq)

    spks = np.array(spks)
    ws = np.array(ws) * 1e-3
    freqs = np.array(freqs)

    wvs = ws * freqs
    lums = spks * d**2.0

    return wvs, lums
    



def make_full_plot():
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    # Add pulsars 
    wv400, lums400 = read_pulsars('pulsars_400.txt', 0.4)
    wv1400, lums1400 = read_pulsars('pulsars_1400.txt', 1.4)
    ax.plot(wv400, lums400, 'k.')
    ax.plot(wv1400, lums1400, 'k.')

    # Add RRATs
    wvr, lumsr = read_rrats('rrats.txt')
    ax.plot(wvr, lumsr, 'r.')

    # B1937+21 GPs -- McKee et al. 2018
    w_m18 = 0.4e-6 # sec
    nu_m18 = 1.4 # GHz 
    d_m18 = 3.27 # kpc
    S_lo_m18 = 5 # Jy
    S_hi_m18 = 125 # Jy
    wv_m18 = w_m18 * nu_m18
    ax.plot(np.ones(2) * wv_m18, 
            np.array([ S_lo_m18, S_hi_m18 ]) * d_m18**2.0,
            c='b', lw=5)
    
    
    # Crab GPs 800 MHz -- Lundgren et al. 1995
    w_c1 = 1.5e-3 # sec
    nu_c1 = 0.8 # GHz 
    d_crab = 2.0 # kpc
    S_lo_c1 = 200 # Jy
    S_hi_c1 = 2000 # Jy
    wv_c1 = w_c1 * nu_c1
    ax.plot(np.ones(2) * wv_c1, 
            np.array([ S_lo_c1, S_hi_c1 ]) * d_crab**2.0,
            c='g', lw=5)
    
    # Crab GPs 1400 MHz -- Karuppusamy et al. 2010
    w_c2 = 10e-6 # sec
    nu_c2 = 1.4 # GHz 
    d_crab = 2.0 # kpc
    S_lo_c2 = 50 # Jy
    S_hi_c2 = 1000 # Jy
    wv_c2 = w_c2 * nu_c2
    ax.plot(np.ones(2) * wv_c2, 
            np.array([ S_lo_c2, S_hi_c2 ]) * d_crab**2.0,
            c='g', lw=5)
    
    # Crab GPs 1700 MHz -- Majid et al. 2011
    w_c3 = 2e-6 # sec
    nu_c3 = 1.7 # GHz 
    d_crab = 2.0 # kpc
    S_lo_c3 = 1e3 # Jy
    S_hi_c3 = 1e5 # Jy
    wv_c3 = w_c3 * nu_c3
    ax.plot(np.ones(2) * wv_c3, 
            np.array([ S_lo_c3, S_hi_c3 ]) * d_crab**2.0,
            c='g', lw=5)


    # Crab nanoshot -- Hankins & Eilek 2007
    w_ns = 0.4e-9 # sec
    S_ns = 2.3e6 # Jy
    d_ns = 2.0 # kpc
    nu_ns = 9.25 # GHz
    ax.plot(w_ns * nu_ns, S_ns * d_ns**2.0, marker='o', 
            color='g', mec='g', ms=10)
    
    # SGR 1935+2154  -- STARE2
    w_sgr1 = 6e-4  # sec
    S_sgr1 = 2.5e6 # Jy
    d_sgr1 = 10 # kpc
    nu_sgr1 = 1.53 # GHz
    ax.plot(w_sgr1 * nu_sgr1, S_sgr1 * d_sgr1**2.0, marker='o', 
            color='k', ms=10, mec='k')
    
    # SGR 1935+2154  -- CHIME
    w_sgr2 = 2e-3  # sec
    S_sgr2 = 1.1e5 # Jy
    d_sgr2 = 10 # kpc
    nu_sgr2 = 0.6 # GHz
    ax.plot(w_sgr2 * nu_sgr2, S_sgr2 * d_sgr2**2.0, marker='o', 
            color='k', ms=10, mec='k')
    
    # SGR 1935+2154  -- Westerbork
    w_sgr3 = np.array([0.86e-3, 0.96e-3])  # sec
    S_sgr3 = np.array([129.0, 25.0])  # Jy
    d_sgr3 = 10 # kpc
    nu_sgr3 = 1.3 # GHz
    ax.plot(w_sgr3 * nu_sgr3, S_sgr3 * d_sgr3**2.0, marker='o', 
            ls='', color='k', ms=10, mec='k')
    
    # SGR 1935+2154  -- FAST
    w_sgr4 = 2e-3  # sec
    S_sgr4 = 30e-3 # Jy
    d_sgr4 = 10 # kpc
    nu_sgr4 = 1.3 # GHz
    ax.plot(w_sgr4 * nu_sgr4, S_sgr4 * d_sgr4**2.0, marker='o', 
            color='k', ms=10, mec='k')

    sgr_y = np.linspace( 0.5 * S_sgr4 * d_sgr4**2., 
                         1.5 * S_sgr1 * d_sgr1**2.0, 100)
    ax.fill_betweenx(sgr_y, x1=0.5*8e-4, x2=1.5*3e-3, color='k', alpha=0.1)
    
    
    
    # FRB 121102 -- Gajjar + Spitler + Scholz
    wv_f0, lum_f0 = read_frb121102('frb121102.txt')
    ax.plot(wv_f0, lum_f0, marker='o', ls='', color='orange')


    # FRB 180916 -- Marthi+2021 and Marcote+2020
    wv_f1, lum_f1 = read_frb180916('frb180916.txt')
    ax.plot(wv_f1, lum_f1, 'co') 


    askap_c = 'purple'
    askap_mec = 'k'


    # FRB 190523 -- Ravi et al. 2019
    w_f2 = 0.4e-3 # sec
    S_f2 = 700 # Jy
    d_f2 = 4.0e6 # kpc (d_l = 4.0 Gpc for z = 0.66)
    nu_f2 = 1.4 # GHz
    ax.plot( w_f2 * nu_f2, S_f2 * d_f2**2.0, marker='s', 
             color=askap_c, mec=askap_mec, ms=10)


    # FRB 190711 -- Macquart et al. 2020
    w_f3 = 6.5e-3 # sec
    S_f3 = 5.2 # Jy
    d_f3 = 3.0e6 # kpc (d_l = 3.0 Gpc for z = 0.52)
    nu_f3 = 1.3 # GHz
    ax.plot( w_f3 * nu_f3, S_f3 * d_f3**2.0, marker='o', 
             color=askap_c, mec=askap_mec, ms=10)


    # FRB 181112 -- Prochaska et al. 2019
    w_f4 = 2.1e-3 # sec
    S_f4 = 12.4 # Jy
    d_f4 = 2.7e6 # kpc (d_l = 2.7 Gpc for z = 0.48)
    nu_f4 = 1.3 # GHz
    ax.plot( w_f4 * nu_f4, S_f4 * d_f4**2.0, marker='o', 
             color=askap_c, mec=askap_mec, ms=10)


    # FRB 190611 -- Macquart et al. 2020
    w_f5 = 2e-3 # sec
    S_f5 = 5 # Jy
    d_f5 = 2.0e6 # kpc (d_l = 2.0 Gpc for z = 0.38)
    nu_f5 = 1.3 # GHz
    ax.plot( w_f5 * nu_f5, S_f5 * d_f5**2.0, marker='o', 
             color=askap_c, mec=askap_mec, ms=10)

    
    # FRB 180924 -- Bannister et al. 2019
    w_f6 = 1.3e-3 # sec
    S_f6 = 12.3 # Jy
    d_f6 = 1.7e6 # kpc (d_l = 1.7 Gpc for z = 0.32)
    nu_f6 = 1.3 # GHz
    ax.plot( w_f6 * nu_f6, S_f6 * d_f6**2.0, marker='o', 
             color=askap_c, mec=askap_mec, ms=10)

    
    # FRB 190102 -- Macquart et al. 2020
    w_f7 = 1.7e-3 # sec
    S_f7 = 8.2 # Jy
    d_f7 = 1.5e6 # kpc (d_l = 1.5 Gpc for z = 0.29)
    nu_f7 = 1.3 # GHz
    ax.plot( w_f7 * nu_f7, S_f7 * d_f7**2.0, marker='o', 
             color=askap_c, mec=askap_mec, ms=10)

    
    # FRB 190608 -- Macquart et al. 2020
    w_f8 = 6e-3 # sec
    S_f8 = 4.3 # Jy
    d_f8 = 560e3 # kpc (d_l = 560 Mpc for z = 0.12)
    nu_f8 = 1.3 # GHz
    ax.plot( w_f8 * nu_f8, S_f8 * d_f8**2.0, marker='o', 
             color=askap_c, mec=askap_mec, ms=10)

    
    # M81 CHIME -- Bhardwaj et al 2021
    w_m81c = 2e-4  # sec
    S_m81c = 10 # Jy
    d_m81c = 3.6e3 # kpc
    nu_m81c = 0.6 # GHz
    ax.plot(w_m81c * nu_m81c, S_m81c * d_m81c**2.0, marker='o', 
            color='r', ms=10, mec='r')

    
    # M81 R -- BURST 
    w_m81b = 33e-6  # sec
    S_m81b = 59 # Jy
    d_m81b = 3.6e3 # kpc
    nu_m81b = 2.3 # GHz
    ax.plot(w_m81b * nu_m81b, S_m81b * d_m81b**2.0, marker='x', 
            color='r', ms=15, mec='r', mew=3)
    
    # M81 R NANOSHOT
    w_m81 = 100e-9  # sec
    S_m81 = 270 # Jy
    d_m81 = 3.6e3 # kpc
    nu_m81 = 2.3 # GHz
    ax.plot(w_m81 * nu_m81, S_m81 * d_m81**2.0, marker='x', 
            color='r', ms=15, mec='r', mew=3)

    # Get B0540 bursts
    wv_B, lum_B = read_B0540('B0540.txt')
    ax.plot(wv_B, lum_B, ls='', mew=3, marker='+', color='LimeGreen', 
            ms=15)


    # Lines of constant Tb
    tb_pows = np.arange(20, 70, 5)
    tb_const = 2.9
    vwx = 10**np.linspace(-10, 10, 100)

    # Thick line for 10^12
    ax.plot(vwx, tb_const * vwx**2.0 * 10**(12-18.), lw=4, c='k')

    # Thinner lines for rest
    for tb_pow in tb_pows:
        ax.plot(vwx, tb_const * vwx**2.0 * 10**(tb_pow-18.), 
                lw=1, ls='-', c='0.7', zorder=-1)

    # Uncertainty Limit 
    uys = 10**np.linspace(-15, 20, 50)
    ax.fill_betweenx(uys, 1e-9, x2=1e-15, color='0.8', zorder=1)

    ax.text(0.05, 0.5, "Uncertainty Principle", fontsize=20, 
            ha='center', va='center', transform=ax.transAxes, 
            rotation=90)


    # Plot labels
    ax.text(3e-3, 1e-6, "Pulsars", fontsize=20, color='k', 
            ha='center', va='center', rotation=0)

    ax.text(5e-2, 1e2, "RRATs", fontsize=20, color='r', 
            ha='center', va='center', rotation=0)
    
    ax.text(1e-5, 4e5, "Crab\nGPs", fontsize=18, color='g', 
            ha='center', va='center', rotation=0)
    
    ax.text(1e-8, 3e5, "Crab\nNanoshot", fontsize=18, color='g', 
            ha='center', va='center', rotation=0)
    
    ax.text(5e-7, 1e0, "B1937+21 \n GPs", fontsize=20, color='b', 
            ha='center', va='center', rotation=0)

    #ax.text(2e-7, 3e11, "This Work", fontsize=20, color='r', 
    #        ha='center', va='center', rotation=0)
    
    ax.text(2e-6, 1e8, "200120E", fontsize=20, color='r', 
            ha='center', va='center', rotation=0)
    
    ax.text(8e-5, 1e7, "CHIME", fontsize=14, color='r', 
            ha='center', va='center', rotation=0)
    
    ax.text(2e-7, 3e10, "Nanoshot", fontsize=16, color='r', 
            ha='center', va='center', rotation=0)
    
    ax.text(1e-5, 3e9, "Burst", fontsize=16, color='r', 
            ha='center', va='center', rotation=0)
    
    ax.text(2e-1, 1e5, "1935+2154", fontsize=20, color='k', 
            ha='center', va='center', rotation=0)
    
    ax.text(3e-5, 8e15, "190523", fontsize=20, color='purple', 
            ha='center', va='center', rotation=0)
    
    ax.text(5e-2, 1e9, "180916", fontsize=20, color='darkcyan', 
            ha='center', va='center', rotation=0)
    
    ax.text(0.2, 1e11, "121102", fontsize=20, color='darkorange', 
            ha='center', va='center', rotation=0)
    
    ax.text(0.2, 1e15, "ASKAP\nFRBs", fontsize=20, color='purple', 
            ha='center', va='center', rotation=0)

    # Label Temps
    ax.text(1, 2e-5, r"$T_{\rm b} = 10^{12}$ K", fontsize=14, color='0.1', 
            ha='center', va='center', rotation=40)
    
    ax.text(5e-5, 3e-6, r"$10^{20}$ K", fontsize=14, color='0.1', 
            ha='center', va='center', rotation=40)
    
    ax.text(6e-7, 6e-5, r"$10^{25}$ K", fontsize=14, color='0.1', 
            ha='center', va='center', rotation=40)
    
    ax.text(2e-8, 4e-3, r"$10^{30}$ K", fontsize=14, color='0.1', 
            ha='center', va='center', rotation=40)

    ax.text(4e-9, 2e1, r"$10^{35}$ K", fontsize=14, color='0.1', 
            ha='center', va='center', rotation=40)
    
    ax.text(1.5e-7, 3e14, r"$10^{45}$ K", fontsize=14, color='0.1', 
            ha='center', va='center', rotation=40)
    
    ax.text(4e-9, 2e16, r"$10^{50}$ K", fontsize=14, color='0.1', 
            ha='center', va='center', rotation=40)

    # Make secondary axis
    ax2 = ax.twinx()

    # Set fig lims
    #ax.set_xlim(1e-10, 1e10)
    #ax.set_ylim(1e-10, 1e16)

    lum_const = 4 * np.pi * 9.5e19
   
    ax1_ylo = 1e-7
    ax1_yhi = 1e18
    ax2_ylo = lum_const * ax1_ylo
    ax2_yhi = lum_const * ax1_yhi
    
    ax.set_xlim(1e-10, 10)


    # set log scale
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    
    ax2.set_yticks( 10.0**np.arange(13, 41, 2))
    ax2.set_ylim(ax2_ylo, ax2_yhi)
    ax2.minorticks_off()
    
    ax.set_yticks( 10.0**np.arange(-7, 19, 2))
    ax.set_ylim(ax1_ylo, ax1_yhi)

    # Axis labels
    ax.set_xlabel(r"$\nu$W (GHz s)", fontsize=18)
    ax.set_ylabel(r"$S_{\rm pk} \, d^2$" + " " +
                  r"$({\rm Jy \, \, kpc}^2)$", fontsize=18)
    
    ax2.set_ylabel(r"$L_{\nu, {\rm iso}}$" + " " +
                  r"$({\rm erg \,\, s}^{-1} \, {\rm Hz}^{-1})$", 
                  fontsize=18, rotation=270, labelpad=30)


    plt.show()
    
    return



