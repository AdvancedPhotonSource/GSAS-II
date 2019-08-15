"""
*Substances: Define Materials*
---------------------------------------------------------------------------------

Defines materials commonly found in small angle & reflectometry experiments.
GSASII substances as a dictionary ''Substances.Substances'' with these materials.

Each entry in ''Substances'' consists of::

     'key':{'Elements':{element:{'Num':float number in formula},...},'Density':value, 'Volume':,value}

Density & Volume are optional, if one missing it is calculated from the other; if both
are missing then Volume is estimated from composition & assuming 10A^3 for each atom,
Density is calculated from that Volume.
See examples below for what is needed.
"""
Substances = {
'Alumina':{'Elements':{'Al':{'Num':2.},'O':{'Num':3.}},'Density':3.986,},
'Water':{'Elements':{'O':{'Num':1.},'H':{'Num':2.}},'Density':1.0},
'Silicon':{'Elements':{'Si':{'Num':8.}},'Volume':160.209},
'a-Quartz':{'Elements':{'Si':{'Num':3.},'O':{'Num':6.}},'Volume':113.057},
'Ethanol':{'Elements':{'C':{'Num':2.},'O':{'Num':1},'H':{'Num':6.}},},
'Polyethylene':{'Elements':{'C':{'Num':1.},'H':{'Num':2.}},'Density':0.93,},
'Polystyrene':{'Elements':{'C':{'Num':1.},'H':{'Num':1.}},'Density':1.060,},
'Teflon':{'Elements':{'C':{'Num':1.},'F':{'Num':2.}},'Density':2.25,},
'Mylar':{'Elements':{'C':{'Num':5.},'H':{'Num':4.},'O':{'Num':2.}},'Density':1.38,},
'Iron':{'Elements':{'Fe':{'Num':4.}},'Density':7.87,},
'FeO-wustite':{'Elements':{'Fe':{'Num':4.},'O':{'Num':4.}},'Volume':79.285},
'Fe2O3-hematite':{'Elements':{'Fe':{'Num':12.},'O':{'Num':18.}},'Volume':301.689},
'Fe3O4-magnetite':{'Elements':{'Fe':{'Num':24.},'O':{'Num':32.}},'Volume':591.921},
'Zirconium':{'Elements':{'Zr':{'Num':2.}},'Density':6.51,},
'Carbon':{'Elements':{'C':{'Num':1.}},'Density':2.27,},
'Titanium':{'Elements':{'Ti':{'Num':1.}},'Density':4.51,},
'TiO2-rutile':{'Elements':{'Ti':{'Num':2.},'O':{'Num':4.}},'Volume':62.452},
'Chromium':{'Elements':{'Cr':{'Num':1.}},'Density':7.19,},
'Nickel':{'Elements':{'Ni':{'Num':4.}},'Density':8.90,},
'Copper':{'Elements':{'Cu':{'Num':4.}},'Density':8.96,},
'Hydroxyapatite':{'Elements':{'Ca':{'Num':5.},'P':{'Num':3.},'O':{'Num':13.},'H':{'Num':1.}},'Density':3.986,},
'Cr2O3':{'Elements':{'Cr':{'Num':2.},'O':{'Num':3.}},'Density':5.206,},
'ZrO2':{'Elements':{'Zr':{'Num':1.},'O':{'Num':3,}},'Density':6.134,},
'Y(0.16)Zr(0.84)O2':{'Elements':{'Y':{'Num':0.16},'Zr':{'Num':0.84},'O':{'Num':2.}},'Density':6.01,},
'Ag':{'Elements':{'Ag':{'Num':1}},'Volume':17.066},
'Al':{'Elements':{'Al':{'Num':1}},'Volume':16.582},
'Au':{'Elements':{'Au':{'Num':1}},'Volume':16.953},
'Co':{'Elements':{'Co':{'Num':1}},'Volume':11.0177},
'FeF2':{'Elements':{'Fe':{'Num':1},'F':{'Num':2}},'Volume':36.352},
'GaAs':{'Elements':{'Ga':{'Num':1},'As':{'Num':1}},'Volume':45.173},
'LaAlO3':{'Elements':{'La':{'Num':1},'Al':{'Num':1},'O':{'Num':3}},'Volume':54.503},
'LaFeO3':{'Elements':{'La':{'Num':1},'Al':{'Num':1},'O':{'Num':3}},'Volume':50.355},
'LaMnO3':{'Elements':{'La':{'Num':1},'Mn':{'Num':1},'o':{'Num':3}},'Volume':58.413},
'MgF2':{'Elements':{'Mg':{'Num':1},'F':{'Num':2}},'Volume':32.58},
'MgO':{'Elements':{'Mg':{'Num':1},'O':{'Num':1}},'Volume':17.977},
'MnF2':{'Elements':{'Mn':{'Num':1},'F':{'Num':2}},'Volume':38.56},
'NiO':{'Elements':{'Ni':{'Num':1},'O':{'Num':1}},'Volume':18.22},
'Pd':{'Elements':{'Pd':{'Num':1}},'Volume':14.738},
'Pt':{'Elements':{'Pt':{'Num':1}},'Volume':15.14},
'SrTiO3':{'Elements':{'Sr':{'Num':1},'Ti':{'Num':1},'O':{'Num':1}},'Volume':26.71},
'V':{'Elements':{'V':{'Num':1}},'Volume':19.26},
'protein':{'Elements':{'C':{'Num':9.25},'N':{'Num':2.34},'O':{'Num':2.77},'S':{'Num':0.22},'H':{'Num':14.3}},'Volume':288.8}, #insulin - typical?
}
# they should not be duplicated in the UserSubstances.py file:
try:
    import UserSubstances as userFile
    Substances.update(userFile.Substances)
except:
    pass
