"""
*Substances: Materials commonly found in small angle & reflectometry experiments*
---------------------------------------------------------------------------------

GSASII substances as a dictionary ''Substances.Substances'' with these materials

Each entry in ''Substances'' consists of:
     'key':{'Elements':{element:{'Num':number in formula},...},'Density':value,
        'Volume':,value}
Density & Volume are optional, if one missing it is calculated from the other; if both
     are missing then Volume is estimated from composition & assuming 10A^3 for each atom,
     Density is calculated from that Volume.
See examples below for what is needed.
"""
Substances = {
'Alumina':{'Elements':{'Al':{'Num':2},'O':{'Num':3}},'Density':3.986,},
'Water':{'Elements':{'O':{'Num':1},'H':{'Num':2}},'Density':1.0},
'Silicon':{'Elements':{'Si':{'Num':8}},'Volume':160.209},
'Ethanol':{'Elements':{'C':{'Num':2},'O':{'Num':1},'H':{'Num':6}},},
}
# they should not be duplicated in the UserSubstances.py file:
try:
    import UserSubstances as userFile
    Substances.update(userFile.Substances)
except:
    pass
