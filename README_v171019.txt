2017 August 18 - WCG Ho

Same as v170818 but extended logTeff grid to 6.8.

35 = NlogTeff, number of effective temperatures
11 = Nlogg, number of surface gravities
166 = NlogEkT, number of photon energies/kTeff
67 = Nmu, number of angles

Column 1: log(E/kTeff) = -1.30,-1.28,...,2.00
Column 2: mu = [1,1.0e-6], where mu = cos(theta),
  every 1.5 degrees and extra three points near 0 and three near 90 deg.
Column 3: specific intensity log(Inu/Teff^3),
  where Inu = [erg s^-1 cm^-2 Hz^-1 ster^-1] and Teff = [K];
  values set to -60 for log(E/kTeff) outside model ranges.
Column 4: logTeff(K) = 5.10,5.15,5.20,...,6.50,6.55,...,6.80
Column 5: logg(cm s^-2) = 13.7,13.8,...,14.7

4281970 rows = 35*11*166*67 = NlogTeff*Nlogg*NlogEkT*Nmu
Order of table:
Row       1: (logTeff,logg,logE/kTeff)=(5.10,13.7,-1.30) and mu=1
...
Row      67: (logTeff,logg,logE/kTeff)=(5.10,13.7,-1.30) and mu=1e-6
Row      68: (logTeff,logg,logE/kTeff)=(5.10,13.7,-1.28) and mu=1
...
Row   11122: (logTeff,logg,logE/kTeff)=(5.10,13.7,+2.00) and mu=1e-6
Row   11123: (logTeff,logg,logE/kTeff)=(5.10,13.8,-1.30) and mu=1
...
Row  122342: (logTeff,logg,logE/kTeff)=(5.10,14.7,+2.00) and mu=1e-6
Row  122343: (logTeff,logg,logE/kTeff)=(5.15,13.7,-1.30) and mu=1
...
Row 3547918: (logTeff,logg,logE/kTeff)=(6.50,14.7,+2.00) and mu=1e-6
...
Row 4281970: (logTeff,logg,logE/kTeff)=(6.80,14.7,+2.00) and mu=1e-6
