"""
*ImageCalibrants: Calibration Standards*
----------------------------------------

GSASII powder calibrants as a dictionary ``ImageCalibrants.Calibrants`` with substances
commonly used for powder calibrations for image data.

Each entry in ``ImageCalibrants`` consists of::

  'key':([Bravais num,],[space group,],[(a,b,c,alpha,beta,gamma),],no. lines skipped,(dmin,pixLimit,cutOff))
  (The space group may be an empty string)

as an example::

  'LaB6  SRM660a':([2,],['',][(4.1569162,4.1569162,4.1569162,90,90,90),],0,(1.0,10,10)),

or where "Bravais num" and "(a,b,...)" are repeated in the case of mixtures:: 

  'LaB6 & CeO2':([2,0],['',''] [(4.1569,4.1569,4.1569,90,90,90),(5.4117,5.4117,5.4117,90,90,90)], 0, (1.0,2,1)),

To expand this list with locally needed additions, do not modify this file,
because you may lose these changes during a software update. Instead
duplicate the format of this file in a file named `UserCalibrants.py`
and there define the material(s) you want::

  Calibrants={
    'LaB6 skip 2 lines':([2,],['',],[(4.1569162,4.1569162,4.1569162,90,90,90),],2,(1.0,10,10)),
  }

New key values will be added to the list of options.
If a key is duplicated, the information in  `UserCalibrants.py` will
override the information in this file. 

Note, some useful Bravais numbers are: F-cubic=0, I-cubic=1, P-cubic=2, R3/m (hex)=3, P6=4, P4mmm=6
"""
Calibrants={
'':([0,],['',],[(0,0,0,0,0,0),],0,(1.0,10,10)),
'LaB6  SRM660b':([2,],[''],[(4.15689,4.15689,4.15689,90,90,90),],0,(1.0,10,10)),
'LaB6  SRM660a':([2,],[''],[(4.1569162,4.1569162,4.1569162,90,90,90),],0,(1.0,10,10)),
'LaB6  SRM660a skip 1':([2,],[''],[(4.1569162,4.1569162,4.1569162,90,90,90),],1,(1.0,10,10)),
'LaB6  SRM660': ([2,],[''],[(4.15695,4.15695,4.15695,90,90,90),],0,(1.0,10,10)),
'Si    SRM640c':([0,],[''],[(5.4311946,5.4311946,5.4311946,90,90,90),],0,(1.,10,10)),
'CeO2  SRM674b':([0,],[''],[(5.411651,5.411651,5.411651,90,90,90),],0,(1.0,2,1)),
'Al2O3 SRM676a':([3,],[''],[(4.759091,4.759091,12.991779,90,90,120),],0,(1.0,5,5)),
'Ni   @ 298K':([0,],[''],[(3.52475,3.52475,3.52475,90,90,90),],0,(1.0,10,10)),
'NaCl @ 298K':([0,],[''],[(5.6402,5.6402,5.6402,90,90,90),],0,(1.0,10,10)),
'NaCl even hkl only':([2,],[''],[(2.8201,2.8201,2.8201,90,90,90),],0,(1.0,10,10)),
'Ag behenate':([6,],[''],[(1.0,1.0,58.380,90,90,90),],0,(7.0,5,1)),
'Spun Si 3600 line/mm grating':([6,],[''],[(1.0,1.0,2777.78,90,90,90),],2,(200.,5,1)),
'Spun Si 7200 line/mm grating':([6,],[''],[(1.0,1.0,1388.89,90,90,90),],1,(200.,5,1)),
'Pt   @ 298K':([0,],[''],[(3.9231,3.9231,3.9231,90,90,90),],0,(1.0,5,1)),
'LaB6 & CeO2':([2,0],['','',],[(4.1569162,4.1569162,4.1569162,90,90,90),(5.411651,5.411651,5.411651,90,90,90)],0,(1.0,2,1)),
}
    
# this should not be duplicated in the UserCalibrants.py file:
try:
    import UserCalibrants as userFile
    Calibrants.update(userFile.Calibrants)
except:
    pass
