import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import jn

def strainGCtest():
    masterSum = np.zeros(99)
    masterSum = np.array([masterSum])
    masterfreq = np.zeros(100)
    masterfreq = np.array([masterfreq])
    yay = 100
    while yay > 0: 
    #everything is in cgs units
        G = 6.67e-8
        c = 3.e8
        m0 = 1.99e33
    #SAm is the mass of sagitarius A*
        SAm = (4.*(10.**6))*m0
    #dist is the distance from the earth to the center of the galaxy
        dist = 2.5e22
    #need to get the master array from my galactic center realization
        numpoints = 1000
    #r0 is the innermost stable circular orbit in pc
        r0 = (3.2*10**-19)*(6.*G*SAm)/(c**2.)
    #for the purposes of this code, R is the periapsis distance
    #    R = np.linspace(r0,1,numpoints)
        start = np.log10(r0)
        stop = np.log10(1.1)
        R = np.logspace(start,stop,numpoints)
    #tau is to transform strain into characteristic strain
        tau = np.sqrt(5*np.pi*10.**7)
    #encMass is the mass enclosed within a certain radius shell in pc
        encMass = []
        for i in range (0,numpoints):
            if R[i] <= (0.25):
                dummy = 1.54e5*(R[i]/.25)**1.7 + 9.e4*(R[i]/.25)**3
            else:
                dummy = 2.2e5*(((R[i]/.25)**1.2)-1) + 1.54e5*(R[i]/.25)**1.7 + 9.e4*(R[i]/.25)**3
            encMass.append(dummy)
        massRadBin = []
        massRadBin.append(encMass[0])

    #now encMass is the enclosed mass at a certain radial bin
    #need to have mass for only that shell
        for i in range (0,numpoints-1):
            massRadBinDum = encMass[i+1] - encMass[i]
            massRadBin.append(massRadBinDum)
    #want to make sure all of my values are of the same type...
        massRadBin = np.array(massRadBin)

        massRadBinTest = []
        massRadBinTest.append(encMass[0])
        for i in range (0,numpoints-1):
            massRadBinTestDum = encMass[i+1] - encMass[i]
            massRadBinTest.append(massRadBinTestDum)

        massRadBinTest = np.array(massRadBinTest)

        def encMassInt(m):
            return m**(-.6)
        massIntTot = quad(encMassInt, 0.1, 100)[0]
    #k is my normalization constant for the IMF

        k = massRadBin/massIntTot
        def Nint(m):
            return m**(-1.6)
        Ninty = quad(Nint,0.1,100)[0]

        nsint = quad(Nint,8,20)[0]
        wdint = quad(Nint,1.1,7.9)[0]
        bhint = quad(Nint,19.9,100)[0]
        starint = quad(Nint,0.1,1)[0]

        nsintfrac = nsint/Ninty
        wdintfrac = wdint/Ninty
        bhintfrac = bhint/Ninty
        starintfrac = starint/Ninty
    #N is the number of stars in each radial bin
        N = k*Ninty
        #and now for my realization loop
        count = 0
        starMass = []
        starMassTot = []
        normConst = 1./(((100**-.6)-(0.1**-.6))/-.6)
        norm1 = 1./Ninty
        bran = []
        aran = []
        semil = []
        ecc = []
        totalCount = 0
        radi = []
        ns = 0
        bh = 0
        wd = 0
        star = 0
        binMass = np.zeros(numpoints)
        binMassCount = []
        for i in range(0,numpoints):
            epsilon = abs(massRadBin[i]-binMass[i])
            testCount = 0
            rad = R[i]
            while epsilon > 0.1:
                ran = np.random.uniform(0.0,1.0)
                e = np.random.uniform(0.0,1.2)
                if e < .9999:
                    a = rad/(1-e)
                elif e >= .9999 and e < 1.01:
                    a = rad
                    e = .9999
                else:
    #a for this is periapsis
                    a = rad/(e-1)
                starMassDum = (((-.6*ran)/norm1)+(0.1**-.6))**(-10./6.)
                starMassDummy = (((-.6*ran)/norm1)+(0.1**-.6))**(-10./6.)
    #            Starmassdum = (Ran/normConst)**(-10./16.)
    #rstar is in solar radii, need to convert to parsecs to compare
                if starMassDum >= 8 and starMassDum < 20:
                    starMassDum = 1.4
                    ns = ns + 1
                    rstar = (3.2e-14)
                    rtidals = (rstar*np.power(((SAm/m0)/starMassDum),1./3.))
                    if rtidals > rad:
                        ns = ns - 1 
                elif starMassDum > 20:
                    starMassDum = 10
                    bh = bh + 1
    #                print 'bh'
                    rtidals = (3.2e-19*(6.*G*SAm)/(c**2.))
                elif starMassDum > 1 and starMassDum < 8:
                    starMassDum = .6
                    wd = wd + 1
                    rstar = (2.26e-10)
                    rtidals = (rstar*np.power(((SAm/m0)/starMassDum),1./3.))
                    if rtidals > rad:
                        wd = wd -1
                else:
                    starMassDum = starMassDum
    #need to convert from solar radii to pc
                    rstar = (((1./4.4e7)*np.power(starMassDum,9./19.)))
                    star = star + 1
                    rtidals = (rstar*np.power(((SAm/m0)/starMassDum),1./3.))
                    if rtidals > rad:
                        star = star - 1
    #            print rtidals
                massRadBinTest[i-1] = massRadBinTest[i-1] - starMassDummy
                if massRadBinTest[i-1] < 0:
                    massRadBinTest[i-1] = massRadBinTest[i-1] + starMassDummy
    #                print 'doody'
                    if epsilon < 0.2:
                        epsilon = 0
                    elif epsilon - 0.2 < .1:
                        epsilon = 0
                elif rtidals > rad:
                    print rtidals
                    massRadBinTest[i-1] = massRadBinTest[i-1] + starMassDum
                else:
                    starMass.append(starMassDum)
                    count = count + 1
                    binMass[i-1] = binMass[i-1] + starMassDummy
                    epsilon = abs(massRadBin[i-1] - binMass[i-1])
                    testCount = testCount + 1
                    aran.append(a)
                    ecc.append(e)
                    radi.append(rad)
    #                print starMassDum, epsilon
            binMassCount.append(testCount)


        aran = np.array(aran)
        ecc = np.array(ecc)
        starMass = np.array(starMass)
        radi = np.array(radi)
    #parsecs to AU
        aranAU = aran*206265.
    
    #    plt.plot(ecc,np.log10(aranAU),'b.')
    #    plt.xlabel('Eccentricity')
    #    plt.ylabel('log(Semi-Major Axis) AU')
    #    plt.show()

    #    fig,ax1 = plt.subplots()
    #    ax1.plot(ecc,np.log10(aranAU),'b.')
    #    ax1.set_xlabel('Eccentricity')
    #    ax1.set_ylabel('log(Semi-Major Axis AU)')
    #    for tl in ax1.get_yticklabels():
    #        tl.set_color('b')

    #    ax2 = ax1.twinx()
    #    ax2.plot(ecc,np.log10(radi*206265.),'r.')
    #    ax2.set_ylabel('log(Periapsis Distance AU)')
    #    for tl in ax2.get_yticklabels():
    #        tl.set_color('r')
    #    plt.show()


        wdfrac = float(wd)/float(count)
        nsfrac = float(ns)/float(count)
        bhfrac = float(bh)/float(count)
        starfrac = float(star)/float(count)

        nserror = (nsfrac-nsintfrac)/(nsintfrac)*100.
        wderror = (wdfrac-wdintfrac)/(wdintfrac)*100.
        bherror = (bhfrac-bhintfrac)/(bhintfrac)*100.
        starerror = (starfrac-starintfrac)/(starintfrac)*100.

#        print 'Neutron Star Error', 'White Dwarf Error', 'Black Hole Error', 'Star Error'
#        print nserror,wderror,bherror,starerror

    #    numbaTest = binMassCount/N[:25]
    #    radTest = R[:25]
    #    plt.plot(radTest,numbaTest,'c.')
    #    plt.ylabel('N(loop)/N(Theory)')
    #    plt.xlabel('Radius')
    #    plt.show()



    #    m = np.linspace(0.1,100,1000)
    #    dndm = norm1*m**(-1.6)
    #    plt.hist(starMass, bins=1000)
    #    plt.plot(m,dndm)
    #    plt.xscale('log')
    #    plt.yscale('log')
    #    plt.xlabel('Mass (Solar Masses)')
    #    plt.ylabel('dn/dm')
    #    plt.show()

        arancm = aran*3.1e18
        radicm = radi*3.1e18

    #if a is in au and t is in years, periods easier to calculate (4pi/G = 1)
        period = []
        periodtest = []
        periodtesthyp = []
        etest = []
        theta = np.linspace(0,2*np.pi,360)
        def semilarea(theta):
            return 1./np.power((b**2/aranAU[i])/(1.-e*np.cos(theta)),2.)
        def posint(the):
            return ahyp*(1.-ecc[i]**2)/(1. + ecc[i]*np.cos(the)) 
    #semil is the semilatus rectum
        for i in range(0,count):
            if ecc[i] >= 1:
    #.99999999
    #            e = .99999
    #            b = aranAU[i]*np.sqrt(1.-e**2.)
    #            res = (aranAU[i]*(1.-e**2.))/(1.+e*np.cos(theta))
    #            semil = aranAU[i]*(ecc[i]**2.-1.)
    #            blardy = quad(semilarea,1.6,4.7)[0]
    #            wholeblardy = quad(semilarea,0,(2.*np.pi))[0]
    #            areac = wholeblardy-blardy
    #            areae = math.pi*aranAU[i]*b
    #            taue = np.power((np.power(aranAU[i],3.)/(SAm/m0 + starMass[i])),1./2.)
    #            tauc = taue*(areac/areae)
    #            periodumdum = tauc*np.pi*10**7. 
                rp = radicm[i]
                ahyp = rp/(ecc[i]-1.)
                r = 2.*posint(np.pi/2.)
    #            nu0 = np.arccos((ahyp*(1.-ecc[i]**2) + rp)/(ecc[i]*(-rp)))
    #            test = (ahyp*(1.-ecc[i]**2) + rp)/(ecc[i]*(-rp))
    #theta is the angle offset from peri (=0)
                F0 = np.arccosh((ecc[i] + np.cos(0))/(1.+ecc[i]*np.cos(0)))
                nu = np.arccos((ahyp*(1.-ecc[i]**2) - r)/(ecc[i]*r))
                testy = (ahyp*(1.-ecc[i]**2) - r)/(ecc[i]*r)
                F = np.arccosh((ecc[i] + np.cos(nu))/(1.+ecc[i]*np.cos(nu)))
                periodum = np.sqrt((ahyp)**3/(G*SAm))*((ecc[i]*np.sinh(F)-F)-(ecc[i]*np.sinh(F0)-F0))
    #            the = np.pi/2.
    #            E = np.arccos((e+np.cos(the))/(1.+(e*np.cos(the))))
    #            x = ahyp*(np.cos(E)-e)
    #            periodumtest = 2.*np.arcsin(np.sqrt(x)-np.sqrt(x*(1.-x)))
    #            periodtest.append(periodumdum)
                periodum = np.power((ecc[i]-1),3./2.)*periodum
    #            periodtesthyp.append(periodum)
    #            etest.append(ecc[i])
            elif ecc[i] == .9999:
                periodum = 2.*np.power(((2.*rp**3)/(G*(starMass[i]*m0+SAm))),1./2.)
            else:
                periodum = np.power((np.power(aranAU[i],3.)/(SAm/m0 + starMass[i])),1./2.)
    #            periodum = np.power((1.-ecc[i]),(3./2.))*periodum
    #            if ecc[i]< 1 and 1 - ecc[i] < .01:
    #                periodum = np.power((1-ecc[i]),3./2.)*periodum
            period.append(periodum)

    
        period = np.array(period)
#        periodtest = np.array(periodtest)
#        periodtesthyp = np.array(periodtesthyp)

     #   print 'Approx', 'Real', 'Eccentricity'
     #   print periodtest[:10], periodtesthyp[:10], etest[:10]

    #    periodtest = np.array(periodtest)
    #    error = (-periodtest+period)/period
    #converting years to seconds
        period = period*(3.14*10.**7)
#        print 'Maximum Period (s), Minimum Period (s)'
#        print max(period), min(period)
    #gravitational wave frequency is 2*f_orb, but this formalism calls for the orbital frequency
        freq = (2./period)
    #    t = 5*np.pi*10.**7

    #    plt.plot(freq,ecc,'m.')
    #    plt.xlabel('frequency (Hz)')
    #    plt.ylabel('Eccentricity')
    #    plt.xscale('log')
    #    plt.show()

    #multiplying by m0 because starMass is in solar masses
        mu = ((starMass*m0)*SAm)/((starMass*m0)+SAm)
    #this quadrupole tensor is for characteristic strain
        Qxx = 6*mu*np.pi**2*((G*starMass*m0*freq)/(4*np.pi))**(2./3.)*((np.sin(2*np.pi*freq)**2)-(np.cos(2*np.pi*freq)**2))
        Qxy = -12*mu*np.pi**2*((G*starMass*m0*freq)/(4*np.pi))**(2./3.)*(np.cos(2*np.pi*freq)*np.sin(2*np.pi*freq))

        Q = ((Qxx**2)+(Qxy**2))**.5   

        h = (2.*G*Q)/((np.power(c,4))*dist)

        chirpm = (G*np.power(SAm,(3./5.))*np.power(starMass*m0,3./5.))/np.power((SAm+starMass*m0),(1./5.))
#        n=1.
#        indh = []
#        gne = []
#        for i in range(0,len(ecc)):
#            jne = jn(n,n*ecc[i])
#            jn_2 = jn(n-2,n*ecc[i])
#            jn_1 = jn(n-1,n*ecc[i])
#            jn1 = jn(n+1,n*ecc[i])
#            jn2 = jn(n+2,n*ecc[i])
#            Bn = jn_2 - (2*ecc[i]*jn_1) + ((2./n)*jne + 2.*ecc[i]*jn1) - jn2
#            An = jn_2 - (2.*jne) + jn2
#            if ecc[i]<1:
#                gnedumb = (n**4./32.)*(Bn**2 + (1.-ecc[i]**2)*An**2 + (4./(3*n**2))*jn1)
#            else:
#                gnedumb = (n**4./32.)*(Bn**2 + (ecc[i]**2-1.)*An**2 + (4./(3*n**2))*jn1)
#            gne.append(gnedumb)
#            indh.append(i)
#        gne = np.array(gne)
#        htest = (G/(np.power(c,4.)))*2.*np.sqrt(32./5.)*(chirpm**(5./3.)/(n*dist))*np.power((2*np.pi*freq),2./3.)*np.sqrt(gne)
        tau = np.pi*5.*10.**7
        maxtest = []
        testfreq = []

        for i in range(0,len(period)):
            maxtestdum = max(period[i],tau)
            maxtest.append(maxtestdum)
        maxtest = np.array(maxtest)
#        htest = htest*np.sqrt(freq*maxtest)
#        error = ((h-htest)/h)*100.
        n = 2.
        htesttest = []
        for i in range(0,len(ecc)):
            if ecc[i] == .9999:
                ecc[i] = 0.9
            if ecc[i] < 1:
                gnetest = (1.+(73/24)*ecc[i]**2. + (37/96)*ecc[i]**4.)/((1-ecc[i]**2.)**(7./2.))
            else:
                gnetest = (1.+(73./24.)*(ecc[i]**2)+(37/96)*(ecc[i]**4.))/((ecc[i]**2.-1)**(7./2.))
            htesttestdum = (G/np.power(c,4.))*2.*np.sqrt(32./5.)*(chirpm[i]**(5./3.)/(n*dist))*np.power((2.*np.pi*freq[i]),2./3.)*np.sqrt(gnetest)
            htesttest.append(htesttestdum)

        htesttest = htesttest*np.sqrt(freq*maxtest)

#        for i in range(0,len(htesttest)):
#            if htesttest[i] >= 10.**-15:
#                print ecc[i], freq[i]

        htesttest = np.array(htesttest)
#        testErrors = ((htest-htesttest)/htest)*100.
#        circError = ((h-htesttest)/h)*100.
#        print testErrors
#        print circError

    #So this error is about 99%. The circular binary approximation grossly overestimates
    #the strength of the gravitational wave strain coming from these eccentric binaries.

    #    plt.plot(np.log10(h),ecc,'k.')
    #    plt.xlabel('log(h)')
    #    plt.ylabel('Eccentricity')
    #    plt.show()
        
    #next step is to bin the frequencies based on their range in logspace

        startfreq = float(np.log10(min(freq)))
        stopfreq = float(np.log10(max(freq)))
        freqBins = np.logspace(startfreq,stopfreq,100)

    #the next section sorts the frequencies from least to greatest values and then sorts the corresponding strains

        inds = freq.argsort()
        sortfreq = freq[inds]
    #sorting characteristic strain with respect to corresponding frequencies
        sorth = htesttest[inds]

    #this double loop is meant to retrieve the strain corresponding to a freqeuncy range defined by the array
    #freqBins, and add up all the strain in that frequency bin

        sumsorth = 0
        sumh = []
        k = 0
        numbin = []
        for i in range(1,len(freqBins)):
            for j in range(0,len(sortfreq)):
                if sortfreq[j] <= freqBins[i] and sortfreq[j] >= freqBins[i-1]:
                    sumsorth = sumsorth + sorth[j]
                    k = k+1
            sumh.append(sumsorth)
            numbin.append(k)
            k = 0
            sumsorth = 0

    #    plt.plot(np.log10(htest),ecc,'m.')
    #    plt.xlabel('log(h)')
    #    plt.ylabel('Eccentricity')
    #    plt.show()

    #    plt.plot(np.log10(freq), ecc, 'r.')
    #    plt.xlabel('log(f) [Hz]')
    #    plt.ylabel('Eccentricity')
    #    plt.show()

    #    plt.plot(np.log10(freq),radi,'k.')
    #    plt.xlabel('log(f) [Hz]')
    #    plt.ylabel('Radius [pc]')
    #    plt.show()
        print len(sumh)
        print len(freqBins)
        yay = yay - 1.
        print yay
#        sumh = np.array([sumh])
#        freqBins = np.array([freqBins])
        masterSum = np.concatenate((masterSum,np.array([sumh])),axis = 0)
        masterfreq = np.concatenate((masterfreq,np.array([freqBins])),axis = 0)
    masterSum = np.array(masterSum)
    masterfreq = np.array(masterfreq)
    stdDev = np.std(masterSum[1:],axis=0)
    stdDevfreq = np.std(masterfreq[1:],axis=0)
    print stdDevfreq
    plt.plot(np.log10(freqBins[1:]), np.log10(sumh), 'c.')
    plt.xlabel('log(f) [Hz]')
    plt.errorbar(np.log10(freqBins[1:]),np.log10(sumh),xerr=np.log10(stdDevfreq[1:]),yerr = np.log10(stdDev))
    plt.ylabel('log(h_c)')
    x,y = np.loadtxt('forgaby_ipta.txt',unpack='True',delimiter=',')
    plt.plot(x,y,'k-')
    plt.savefig('meanStrainTest.png')
#    plt.show()



strainGCtest()
