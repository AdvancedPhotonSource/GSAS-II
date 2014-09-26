# define some default instrument parameter files
# just like GSAS, sigh
defaultIparm_lbl = []
defaultIparms = []
defaultIparm_lbl.append('CuKa lab data')
defaultIparms.append({
    'INS   HTYPE ':'PXC ',
    'INS  1 ICONS':'  1.540500  1.544300       0.0         0       0.7    0       0.5   ',
    'INS  1PRCF1 ':'    3    8      0.01                                                ',
    'INS  1PRCF11':'   2.000000E+00  -2.000000E+00   5.000000E+00   0.000000E+00        ',
    'INS  1PRCF12':'   0.000000E+00   0.000000E+00   0.150000E-01   0.150000E-01        ',
    })
defaultIparm_lbl.append('0.6A synch')
defaultIparms.append({
    'INS   HTYPE ':'PXC ',
    'INS  1 ICONS':'  0.600000  0.000000       0.0         0      0.99    0       0.5   ',
    'INS  1PRCF1 ':'    3    8      0.01                                                ',
    'INS  1PRCF11':'   1.000000E+00  -1.000000E+00   0.300000E+00   0.000000E+00        ',
    'INS  1PRCF12':'   0.000000E+00   0.000000E+00   0.100000E-01   0.100000E-01        ',
    })
defaultIparm_lbl.append('1.5A CW neutron data')
defaultIparms.append({
    'INS   HTYPE ':'PNC',
    'INS  1 ICONS':'   1.54020   0.00000   0.04000         0',
    'INS  1PRCF1 ':'    3    8     0.005',
    'INS  1PRCF11':'   0.239700E+03  -0.298200E+03   0.180800E+03   0.000000E+00',
    'INS  1PRCF12':'   0.000000E+00   0.000000E+00   0.400000E-01   0.300000E-01',
    })
defaultIparm_lbl.append('10m TOF backscattering bank')
defaultIparms.append({
    'INS   FPATH1':'      9.00',
    'INS   HTYPE ':'PNT',
    'INS  1 ICONS':'   5000.00      0.00      0.00',
    'INS  1BNKPAR':'    1.0000   150.000',       
    'INS  1PRCF1 ':'    1    8   0.01000',
    'INS  1PRCF11':'   0.000000E+00   5.000000E+00   3.000000E-02   1.000000E-03',
    'INS  1PRCF12':'   0.000000E+00   4.000000E+01   0.000000E+00   0.000000E+00',        
    })
defaultIparm_lbl.append('10m TOF 90deg bank')
defaultIparms.append({
    'INS   FPATH1':'      9.00',
    'INS   HTYPE ':'PNT',
    'INS  1 ICONS':'   3500.00      0.00      0.00',
    'INS  1BNKPAR':'    1.0000    90.000',       
    'INS  1PRCF1 ':'    1    8   0.01000',
    'INS  1PRCF11':'   0.000000E+00   5.000000E+00   3.000000E-02   4.000000E-03',
    'INS  1PRCF12':'   0.000000E+00   8.000000E+01   0.000000E+00   0.000000E+00',        
    })
defaultIparm_lbl.append('63m POWGEN 90deg bank')
defaultIparms.append({
    'INS   FPATH1':'     60.00',
    'INS   HTYPE ':'PNT',
    'INS  1 ICONS':'  22585.80      0.00      0.00',
    'INS  1BNKPAR':'     3.169    90.000',       
    'INS  1PRCF1 ':'    1    8   0.01000',
    'INS  1PRCF11':'   0.000000E+00   1.000000E+00   3.000000E-02   4.000000E-03',
    'INS  1PRCF12':'   0.000000E+00   8.000000E+01   0.000000E+00   0.000000E+00',        
    })
