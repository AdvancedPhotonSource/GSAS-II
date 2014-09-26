#!/usr/bin/env python
# -*- coding: utf-8 -*-
#defaultIparms
########### SVN repository information ###################
# $Date: 2014-09-26 10:52:52 -0500 (Fri, 26 Sep 2014) $
# $Author: vondreele $
# $Revision: 1507 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/defaultIparms.py $
# $Id: GSASII.py 1507 2014-09-26 15:52:52Z vondreele $
########### SVN repository information ###################
'''
Define some default instrument parameters:
Format for each is a list of strings finished with a '\n'.
Begin with '#GSAS-II...' as the reader routine checks this.
Each line can be comprised of a block of ';' delimited name:value pairs.
All instrument parameters must be included; even those = 0.
Use a GSAS-II instprm file as a source for the entries.
For a new entry:
    Append a useful name to defaultIparms_lbl.
    Append the list of lines to defaultIparms.
See examples below.
'''

defaultIparm_lbl = []
defaultIparms = []

defaultIparm_lbl.append('CuKa lab data')
defaultIparms.append([
    '#GSAS-II instrument parameter file for lab CuKa data\n',
    'Type:PXC\n',
    'Lam1:1.5405;Lam2:1.5443;Zero:0.0;Polariz.:0.7;Azimuth:0.0;I(L2)/I(L1):0.5\n',
    'U:2.0;V:-2.0;W:5.0;X:0.0;Y:0.0;SH/L:0.002\n',
])

defaultIparm_lbl.append('APS 30keV 11BM')
defaultIparms.append([
    '#GSAS-II instrument parameter file APS 11BM @ 30keV\n', 
    'Type:PXC\n',
    'Lam:0.413263;Polariz.:0.99;Azimuth:0.0;Zero:0.0\n',
    'U:1.163;V:-0.126;W:0.063;X:0.0;Y:0.0;SH/L:0.002\n',
])

defaultIparm_lbl.append('0.7A synchrotron data')
defaultIparms.append([
    '#GSAS-II instrument parameter file 0.7A synchrotron data\n',
    'Type:PXC\n',
    'Lam:0.69968;Zero:0.0;Polariz.:0.99;Azimuth:0.0\n',
    'U:5.9840407759;V:-1.28771353531;W:0.118521878603\n',
    'X:-0.0977791308891;Y:4.40147397286;SH/L:0.0264356231583\n',
])

defaultIparm_lbl.append('1.9A ILL D1A CW data')
defaultIparms.append([
    '#GSAS-II instrument parameter file\n',
    'Type:PNC\n',
    'Lam:1.909;Zero:0.0;Polariz.:0.0;Azimuth:0.0\n',
    'U:257.182710995;V:-640.525145369;W:569.378664828\n',
    'X:0.0;Y:0.0;SH/L:0.002\n',
])

defaultIparm_lbl.append('9m HIPD 151deg bank TOF data')
defaultIparms.append([
    '#GSAS-II instrument parameter file for 9m HIPD back scattering bank\n',
    'Type:PNT\n',
    'fltPath:10.32567;2-theta:151.0;Azimuth:0.0\n',
    'Zero:-0.773346536757;difC:5084.82763065;difA:-2.6304177486;difB:0.0\n',
    'alpha:5.0\n',
    'beta-0:0.0332763989665;beta-1:0.000964057827372;beta-q:0.0\n',
    'sig-0:0.0;sig-1:15.1402867268;sig-2:0.0;sig-q:0.0\n',
    'X:0.0;Y:0.0\n',
])

defaultIparm_lbl.append('10m TOF 90deg bank')
defaultIparms.append([
    '#GSAS-II instrument parameter file for 10m TOF 90deg bank\n',
    'Type:PNT\n',
    'fltPath:10;2-theta:90.0;Azimuth:0.0\n',
    'Zero:0.0;difC:3500.;difA:0.0;difB:0.0\n',
    'alpha:5.0\n',
    'beta-0:0.03;beta-1:0.004;beta-q:0.0\n',
    'sig-0:0.0;sig-1:80.0;sig-2:0.0;sig-q:0.0\n',
    'X:0.0;Y:0.0\n',
])

defaultIparm_lbl.append('63m POWGEN 90deg bank')
defaultIparms.append([
    '#GSAS-II instrument parameter file for POWGEN\n',
    'Type:PNT;fltPath:63.169;2-theta:90.0;Azimuth:0.0\n',
    'Zero:-4.96558487231;difC:22594.7440533;difA:-0.927945556608;difB:1.42511277627\n',
    'alpha:1.0\n',
    'beta-0:0.138077840635;beta-1:0.0029606795286;beta-q:0.0\n',
    'sig-0:24.8202075678;sig-1:-82.07196132;sig-2:269.925504862;sig-q:0.0\n',
    'X:-1.80259010604;Y:4.47209435997\n',
])
