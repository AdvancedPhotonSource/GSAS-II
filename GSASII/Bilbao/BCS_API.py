# -*- coding: utf-8 -*-
"""
Created on 09/03/2026
A project to integrate automated running of the Bilbao Crystallographic Server
(https://cryst.ehu.es/)

This is an ongoing project, support for new tools will be added over time.

@author: sururi
"""
#%%
import os
import json
import requests
import re
import base64
import threading

# A single global used for the entire module
BCS_container = {'result': None,       # container for thread return value
                 'BCS_SERVER': None,   # server URL
                 'BCS_API_KEY': None,  # key to access server
                 'BCS_Sleep': None,    # callback sleep routine while waiting
                }

def initAPI(server=None,api_key=None,threadSleep=None):
    '''Check that the API is initialized, or initialize it, if called 
    with appropriate values. If the API is not initialized, attempt to 
    read from a .json file with this information in a few locations. 

    Throw an exception if called and the API is not initialized. 

    This routine is designed to be "low cost" so that it can be called
    from BCS_submit everytime it accesses the server. This means that 
    if the .json file is used, one can call the appropriate BCS_ routines
    without initializing, but one can override that by initializing 
    directly, as long as it is done before calling any of the BCS_
    routines.

    example use of threading for wxPython::

      def my_sleep():
          wx.GetApp().Yield()  # or your preferred event loop call
          time.sleep(0.1)
      initAPI(threadSleep=my_sleep)
    '''
    config_error_msg = """
The BCS API cannot be used without defining the BCS_API key and BCS server
address. This can be done in a call to initAPI or can be defined 
in a file named .BCSconfig.json that is located in your home directory 
or the location where the BCS_API.py module is located. Modify the
"config_template.json" file accordingly and save it as ".BCSconfig.json"
"""
    init = False
    if server:
        BCS_container['BCS_SERVER'] = server
        init = True
    if api_key:
        BCS_container['BCS_API_KEY'] = api_key
        init = True
    if threadSleep:
        init = True
        BCS_container['BCS_Sleep'] = threadSleep

    if init: return
    if (BCS_container['BCS_SERVER'] is not None
            and BCS_container['BCS_API_KEY'] is not None): return

    # If we have not initialized previously, set setup information from a file
    #   Look in the user's home directory or the location of this file
    #   TODO?: consider expanding the allowed locations? Environment vars?
    loc_list = (os.path.expanduser('~'), os.path.dirname(__file__))
    for d in loc_list:
        fil = os.path.join(d,'.BCSconfig.json')
        if os.path.exists(fil):
            try:
                with open(fil, 'r') as f:
                    config = json.load(f)
            except:
                print(f'read of {fil} failed')
                continue
        else:
            continue
        try:
            BCS_container['BCS_SERVER'] = config['BCS_SERVER']
            BCS_container['BCS_API_KEY'] = config['BCS_API_KEY']
        except:
            print(f"SERVER and/or API_KEY not defined in {fil}.")

    if (BCS_container['BCS_SERVER'] is None
            or BCS_container['BCS_API_KEY'] is None):
        print (config_error_msg)
        raise FileNotFoundError(".BCSconfig.json file not found.")

def validate_int(value, min_value='', max_value=''):
    if type(value)!=int:
        if(not re.match(r'^\d+$', value)):
            return False
        value = int(value)
    if min_value and value < min_value:
        return False
    if max_value and value > max_value:
        return False
    return True

def validate_number_fraction(value):
    value = str(value) # re is not applicable to numbers
    if not re.match('^[0123456789+-/.]+$',value):
        return False
    return True

def validate_yesno(value):
    if value not in ['yes','no']:
        return False
    return True

def validate_xyz(symop):
    """
    Ultra-simple check for XYZ format.
    """
    # Remove whitespace
    symop = symop.strip()

    # Check for exactly 2 commas
    if symop.count(',') != 2:
        return False

    # Split and check all parts exist
    parts = [p.strip() for p in symop.split(',')]
    if len(parts) != 3:
        return False

    # All parts must contain something
    if not all(parts):
        return False

    # Quick check for valid characters
    valid = 'xyzXYZ0123456789+-/., '
    for char in symop:
        if char not in valid:
            return False
    return True
def validate_msg2(msg):
    # Checks the BNS style format (###.###)
    if re.match(r'^\d+\.\d+$', msg):
        return True
    return False

#%% BCS functions
def BCS_submit(payload,write_to_file=''):
    initAPI()

    def BCS_request():
        '''Put in a Post request to the BCS and wait for the server response
        '''
        # Set up headers
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {BCS_container['BCS_API_KEY']}"
        }
        try:
            # Make the POST request
            response = requests.post(
                BCS_container['BCS_SERVER'],
                headers=headers,
                data=json.dumps(payload),
                timeout=300  # 30 second timeout
            )

            # Check if request was successful
            response.raise_for_status()

            # Parse JSON response
            result = response.json()
            if(write_to_file):
                with open(write_to_file, 'w') as f:
                    f.write(result['result']['output'])
            BCS_container['result'] = result

        except requests.exceptions.RequestException as e:
            print(f"Error making request: {e}")
            if hasattr(e, 'response') and e.response is not None:
                print(f"Response status: {e.response.status_code}")
                print(f"Response body: {e.response.text}")
            BCS_container['result'] = None
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON response: {e}")
            print(f"Raw response: {response.text}")
            BCS_container['result'] = None
    
    # Check if thread sleep function is defined
    if BCS_container['BCS_Sleep'] is None:
        BCS_request()         # Run synchronously (non-threaded)
    else:
        # Run the POST request in a background thread
        thread = threading.Thread(target=BCS_request)
        thread.start()        
        # Loop calling BCS_container['BCS_Sleep'] until thread finishes
        while thread.is_alive():
            BCS_container['BCS_Sleep']()
        thread.join()  # Ensure thread is fully joined
    
    return BCS_container['result']

def BCS_GENPOS(sgnum,write_to_file=''):
    if(not validate_int(sgnum,1,230)):
        print ("Space group number must be an integer between 1 and 230.")
        return False
    data = {
        "inputs": {
            "sgnum": str(sgnum)
        },
        "tool": 'GENPOS',
    }
    result = BCS_submit(data, write_to_file)
    if result:
        return result['result']['output']
    return None

def BCS_WYCKPOS(sgnum,write_to_file=''):
    if(not validate_int(sgnum,1,230)):
        print ("Space group number must be an integer between 1 and 230.")
        return False
    data = {
        "inputs": {
            "sgnum": str(sgnum)
        },
        "tool": 'WYCKPOS',
    }
    result = BCS_submit(data,write_to_file)
    if result:
        #print (json.dumps(result))
        return result['result']['output']
    return None

def BCS_HKLCOND(sgnum,write_to_file=''):
    if(not validate_int(sgnum,1,230)):
        print ("Space group number must be an integer between 1 and 230.")
        return False
    data = {
        "inputs": {
            "sgnum": str(sgnum)
        },
        "tool": 'HKLCOND',
    }
    result = BCS_submit(data,write_to_file)
    if result:
        #print (json.dumps(result))
        return result['result']['output']
    return None

def BCS_MAXSUB(sgnum,write_to_file=''):
    if(not validate_int(sgnum,1,230)):
        print ("Space group number must be an integer between 1 and 230.")
        return False
    data = {
        "inputs": {
            "sgnum": str(sgnum)
        },
        "tool": 'MAXSUB',
    }
    result = BCS_submit(data,write_to_file)
    if result:
        return result['result']['output']
    return None

def BCS_SERIES(sgnum,write_to_file=''):
    if(not validate_int(sgnum,1,230)):
        print ("Space group number must be an integer between 1 and 230.")
        return False
    data = {
        "inputs": {
            "sgnum": str(sgnum)
        },
        "tool": 'SERIES',
    }
    result = BCS_submit(data,write_to_file)
    if result:
        return result['result']['output']
    return None

def BCS_WYCKSETS(sgnum,write_to_file=''):
    if(not validate_int(sgnum,1,230)):
        print ("Space group number must be an integer between 1 and 230.")
        return False
    data = {
        "inputs": {
            "sgnum": str(sgnum)
        },
        "tool": 'WYCKSETS',
    }
    result = BCS_submit(data,write_to_file)
    if result:
        return result['result']['output']
    return None

def BCS_IDENTIFYGROUP(generators,write_to_file=''):
    # Generators: list of xyz ops.
    for generator in generators:
        if not validate_xyz(generator):
            print (f"{generator} is invalid")
            return False
        # else:
        #     print(f"{generator} is valid")
    data = {
        "inputs": {
            "generators": generators
        },
        "tool": 'IDENTIFY_GROUP',
    }
    result = BCS_submit(data,write_to_file)
    if result:
        return result['result']['output']
    return None

def BCS_kSUBGROUPSMAG(super='', # serial number of the space group of the parent paramagnetic phase
                      generators='', # generators defining the parent paramagnetic phase
                      km1x='0', # magnetic wave vectors
                      km1y='0',
                      km1z='0',
                      starmagnetica='no', # Choose the whole star of the propagation vector?
                      celtodas='no', #  Include the subgroups compatible with intermediate cells?
                      landau='no', # Show only subgroups that can be the result of a Landau-type         transition (single irrep order parameter)?
                      limite='spgroup', #  possible limitations of the subgroup list
                      # spgroup : Lowest space group to consider ('sub')
                      #  pgroup : Lowest point group to consider ('pointgroup': '0':all, ..., '122':m'-  3'm')
                      # crysclas: Lowest crystal system to consider ('crystalsystem': '0':all, ..., '7': cubic)
                      # maximal : Only maximal subgroups
                      sub='1.1', # Lowest space group to consider ('sub')
                      pointgroup='0', # Lowest point group to consider ('0':all, ..., '122':m'-3'm')
                      crystalsystem='0', # Lowest crystal system to consider ('0':all, ..., '7':cubic)
                      centrosymmetry='0', # Only centrosymmetric / non-centrosymmetric groups ('0': all, '1': centrosymmetric, '2': non-centrosymmetric)
                      polarity='0', # Only polar / non-polar groups ('0': all, '1': polar, '2': non-     polar)
                      listado='lista', # List / Graph of groups ('lista': list, 'graf': graph)
                      write_to_file='',
                      GSAS=False):
    inputs = {}
    # Generators: list of xyz ops.
    if len(generators):
        for generator in generators:
            if not validate_xyz(generator):
                print (f"{generator} is invalid")
                return False
        inputs['generators'] = generators
    elif super and validate_int(super,1,230):
        inputs['super'] = super
    else:
        print("Generators or parent space group not properly defined.")
        return False
    for k in [km1x, km1y, km1z]:
        if not validate_number_fraction(k):
            print("k-vectors are not properly defined.")
            return False
    inputs['km1x'] = km1x
    inputs['km1y'] = km1y
    inputs['km1z'] = km1z
    if not validate_yesno(starmagnetica):
        print("Improper star option.")
        return False
    inputs['starmagnetica'] = starmagnetica
    if not validate_yesno(celtodas):
        print("Improper intermediate cells option.")
        return False
    inputs['celtodas'] = celtodas
    if not validate_yesno(landau):
        print("Improper Landau-type transition option.")
        return False
    inputs['landau'] = landau
    if limite not in ['spgroup', 'pgroup', 'crysclas', 'maximal']:
        print("Improper limitation option.")
        return False
    inputs['limite'] = limite
    if not validate_msg2(sub):
        print("Improper lowest space group declaration.")
        return False
    inputs['sub'] = sub
    if not validate_int(pointgroup,0,122):
        print("Improper pointgroup specification.")
        return False
    inputs['pointgroup'] = pointgroup
    if not validate_int(crystalsystem,0,7):
        print("Improper crystal system specification.")
        return False
    inputs['crystalsystem'] = crystalsystem
    if not validate_int(centrosymmetry,0,2):
        print("Improper centrosymmetry specification.")
        return False
    inputs['centrosymmetry'] = centrosymmetry
    if not validate_int(polarity,0,2):
        print("Improper polarity specification.")
        return False
    inputs['polarity'] = polarity
    if listado not in ['lista','graf']:
        print("Improper listing/graph specification.")
        return False
    inputs['listado'] = listado
    if GSAS:
        tool = 'kSUBGROUPSMAG_GSAS'
    else:
        tool = 'kSUBGROUPSMAG'
    data = {
        "inputs": inputs,
        "tool": tool,
    }
    result = BCS_submit(data,write_to_file)
    if result:
        return(result['result']['output'])
    return None
def BCS_SUBGROUPS(super='', # serial number of the space group of the parent paramagnetic phase
                      generators='', # generators defining the parent paramagnetic phase
                      knm1x='0', # magnetic wave vectors
                      knm1y='0',
                      knm1z='0',
                      starmagnetica='no', # Choose the whole star of the propagation vector?
                      celtodas='no', #  Include the subgroups compatible with intermediate cells?
                      landau='no', # Show only subgroups that can be the result of a Landau-type         transition (single irrep order parameter)?
                      limite='spgroup', #  possible limitations of the subgroup list
                      # spgroup : Lowest space group to consider ('sub')
                      #  pgroup : Lowest point group to consider ('pointgroup': '0':all, ..., '122':m'-  3'm')
                      # crysclas: Lowest crystal system to consider ('crystalsystem': '0':all, ..., '7': cubic)
                      # maximal : Only maximal subgroups
                      sub='1', # Lowest space group to consider ('sub')
                      pointgroup='0', # Lowest point group to consider ('0':all, ..., '122':m'-3'm')
                      crystalsystem='0', # Lowest crystal system to consider ('0':all, ..., '7':cubic)
                      centrosymmetry='0', # Only centrosymmetric / non-centrosymmetric groups ('0': all, '1': centrosymmetric, '2': non-centrosymmetric)
                      polarity='0', # Only polar / non-polar groups ('0': all, '1': polar, '2': non-     polar)
                      strain='0', # Only proper ferroelastic phase transitions ('0': no, '1': yes)
                      listado='lista', # List / Graph of groups ('lista': list, 'graf': graph)
                      write_to_file='',
                      GSAS=False):
    inputs = {}
    # Generators: list of xyz ops.
    if len(generators):
        for generator in generators:
            if not validate_xyz(generator):
                print (f"{generator} is invalid")
                return False
        inputs['generators'] = generators
    elif super and validate_int(super,1,230):
        inputs['super'] = super
    else:
        print("Generators or parent space group not properly defined.")
        return False
    for k in [knm1x, knm1y, knm1z]:
        if not validate_number_fraction(k):
            print("k-vectors are not properly defined.")
            return False
    inputs['knm1x'] = knm1x
    inputs['knm1y'] = knm1y
    inputs['knm1z'] = knm1z
    if not validate_yesno(starmagnetica):
        print("Improper star option.")
        return False
    inputs['starmagnetica'] = starmagnetica
    if not validate_yesno(celtodas):
        print("Improper intermediate cells option.")
        return False
    inputs['celtodas'] = celtodas
    if not validate_yesno(landau):
        print("Improper Landau-type transition option.")
        return False
    inputs['landau'] = landau
    if limite not in ['spgroup', 'pgroup', 'crysclas', 'maximal']:
        print("Improper limitation option.")
        return False
    inputs['limite'] = limite
    if not validate_int(sub,1,230):
        print("Improper lowest space group declaration.")
        return False
    inputs['sub'] = sub
    if not validate_int(pointgroup,0,122):
        print("Improper pointgroup specification.")
        return False
    inputs['pointgroup'] = pointgroup
    if not validate_int(crystalsystem,0,7):
        print("Improper crystal system specification.")
        return False
    inputs['crystalsystem'] = crystalsystem
    if not validate_int(centrosymmetry,0,2):
        print("Improper centrosymmetry specification.")
        return False
    inputs['centrosymmetry'] = centrosymmetry
    if not validate_int(polarity,0,2):
        print("Improper polarity specification.")
        return False
    inputs['polarity'] = polarity
    if not validate_int(strain,0,1):
        print("Improper strain specification.")
        return False
    inputs['strain'] = strain
    if listado not in ['lista','graf']:
        print("Improper listing/graph specification.")
        return False
    inputs['listado'] = listado
    if GSAS:
        tool = 'SUBGROUPS_GSAS'
    else:
        tool = 'SUBGROUPS'
    data = {
        "inputs": inputs,
        "tool": tool,
    }
    result = BCS_submit(data,write_to_file)
    if result:
        return result['result']['output']
    return None

def BCS_PSEUDO(stru='', # Initial Structure (LS)
               maxdelta = '2', # Maximum allowed distance for pseudosymmetry search
               what = 'minsup', # Supergroups type for pseudosymmetry search
               # minsup : Minimal supergroups ('showindex')
               # celsup : Supergroups with specified k-index ('maxik')
               # mysuper : Specify supergroup transformation ('mysuper2', '[xyz][1-4]')
               # pseudocell : For monoclinic and triclinic structures: previous check of lattice pseudosymmetry ('angtol')
               showindex = 'no', # Show only indices in supergroups table?
               maxik = '1', # Supergroups with specified k-index ('1',...,'4')
               mysuper2 = '221', # Specified supergroup ('1',...,'230')
               # Transformation matrix from the specified supergroup
               x1='1',x2='0',x3='0',x4='0',
               y1='0',y2='1',y3='0',y4='0',
               z1='0',z2='0',z3='1',z4='0',
               trmat_abc = '', # alternative notation for the transformation matrix
                               # (overrides {xyz}{1234) if defined
               angtol = '5', # Angular tolerance for pseudolattice check,
               origin = 'yes', # Whether to search for origin shift in polar cases
               polares = '0.4', # Grid search length unit for polar cases
               write_to_file = ''):
    inputs = {}
    stru = stru.strip()
    if stru == '':
        print("Structure information is missing.")
        return False
    inputs['stru'] = stru
    if what not in ['minsup', 'celsup', 'mysuper', 'pseudocell']:
        print("Improper method definition.")
        return False
    inputs['what'] = what
    if not validate_number_fraction(maxdelta):
        print("Maximum delta is not properly defined.")
        return False
    inputs['maxdelta'] = maxdelta
    if not validate_yesno(showindex):
        print("Improper show index option.")
        return False
    inputs['showindex'] = showindex
    if not validate_int(maxik,1,4):
        print("Maximum ik option is not properly defined.")
        return False
    inputs['maxik'] = maxik
    if not validate_int(mysuper2,1,230):
        print("Maximum supergroup number is not properly defined.")
        return False
    inputs['mysuper2'] = mysuper2
    for xyz1234 in [x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4]:
        if not validate_number_fraction(xyz1234):
            print("Component(s) of the transformation matrix are not properly defined.")
            return False
    inputs['x1'] = x1
    inputs['x2'] = x2
    inputs['x3'] = x3
    inputs['x4'] = x4
    inputs['y1'] = y1
    inputs['y2'] = y2
    inputs['y3'] = y3
    inputs['y4'] = y4
    inputs['z1'] = z1
    inputs['z2'] = z2
    inputs['z3'] = z3
    inputs['z4'] = z4
    if(trmat_abc and not re.match('^[abc0123456789+-/.; rst]+$',trmat_abc)):
        print("Transformation matrix (abc format) is not properly defined.")
        return False
    trmat_abc = re.sub(r'[rst]', '0', trmat_abc)
    inputs['trmat_abc'] = trmat_abc
    if not validate_number_fraction(angtol):
        print("Angle tolerance is not properly defined.")
        return False
    inputs['angtol'] = angtol
    if not validate_yesno(origin):
        print("Origin shift option is not properly defined.")
        return False
    inputs['origin'] = origin
    if not validate_number_fraction(polares):
        print("Origin shift unit for polar structures is not properly defined.")
        return False
    inputs['polares'] = polares
    tool = 'PSEUDO'
    data = {
        "inputs": inputs,
        "tool": tool,
    }
    #print(json.dumps(data))

    result = BCS_submit(data,write_to_file)
    if result:
        return result['result']['output']
    return None

def BCS_CIF2STANDARD(input_CIF_file_path='', # Input CIF file path
                     output_CIF_file_path='', # If designated, the resulting CIF will be written here
                     write_to_file='',
                     response_format = 'json' # Returns with fields:
                                               # 'trmat': Transformation matrix relating the input setting to the standard
                                               # 'BCS': Structure definition in BCS format (standard setting)
                                               # 'CIF': Structure definition (standard setting)
                                               # 'SG_symbol': Space group symbol of the standard setting
                                               # 'SG_number': Space group ITA classification number
                    # If you want the direct html output, then set response_format to 'html'
                   ):
    with open(input_CIF_file_path, 'rb') as f:
        file_content = f.read()
        file_base64 = base64.b64encode(file_content).decode('utf-8')

    inputs = {
            "filename": input_CIF_file_path.split('/')[-1],
            "file_data": file_base64,
            "file_type": "application/octet-stream"
    }
    tool = 'CIF2STANDARD'
    data = {
        "inputs": inputs,
        "tool": tool,
        "response_format": response_format
    }
    #print(json.dumps(data))
    result = BCS_submit(data,write_to_file)
    if result:
        if(response_format == 'json'):
            r = json.loads(result['result']['output'])
            result['result']['output'] = r
            if(output_CIF_file_path != ''):
                with open(output_CIF_file_path, 'w') as f:
                    f.write(r['CIF'])
            return result['result']['output']
        return result['result']['output']
    return None

def parse_PSEUDO_results(pseudo_result):
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(pseudo_result, 'lxml')

    #region RemoveFullDetails
    # Remove the details to prevent double occurence == 0
    html_str = str(soup)
    stop_idx = html_str.find('Pseudosymmetry search full report')

    if stop_idx != -1:
        # Create a new soup from the HTML before the stop point
        html_before = html_str[:stop_idx]
        soup = BeautifulSoup(html_before, 'lxml')
    #endregion

    dict = {}

    search_results_table = soup.find('h3', string=re.compile(
        'Summary search results'))
    table = search_results_table.find_next('table')
    for row in table.find_all('tr'):
        row_list = []
        for col in row.find_all('td'):
            row_list.append(col.text)
        if (len(row_list)):
            i = int(row_list[0])
            symb, sgnum = row_list[1].split(" ")
            sgnum = re.sub(r'\(([0-9]+)\)', r'\1', sgnum)
            sgnum = int(sgnum)
            if (row_list[6] == '>tol'):
                row_list[6] = '-'
            dict_aux = {}
            dict_aux['symb'] = symb
            dict_aux['sgnum'] = sgnum
            dict_aux['index'] = row_list[2]
            dict_aux['i_k'] = row_list[3]
            dict_aux['trmat'] = row_list[4]
            dict_aux['Delta'] = row_list[6]
            dict_aux['Disp'] = row_list[7]
            dict_aux['subgroup_setting'] = None
            dict_aux['supergroup_setting'] = None

            dict[i] = dict_aux

    for supergroup_text in soup.find_all(
            lambda tag: tag.name == 'h4' and 'Supergroup' in tag.get_text()):
        stru_num = re.sub(r'^([\d]+)#.*', r"\1", supergroup_text.text)
        stru_num = int(stru_num)
        title_sub = supergroup_text.find_next('b', string=re.compile(
            r'\(subgroup setting\)'))
        sub_pre_block = title_sub.find_next('pre')
        dict[stru_num]['subgroup_setting'] = sub_pre_block.text
        title_super = supergroup_text.find_next('b', string=re.compile(
            r'\(supergroup setting\)'))
        super_pre_block = title_super.find_next('pre')
        dict[stru_num]['supergroup_setting'] = super_pre_block.text

    return dict

def parse_PSEUDO_list(pseudo_result):
    """
    Parses the result list of a PSEUDO minsup call.

    Returns a dictionary with the following keys per each list item:
    dict_keys(['symb', 'sgnum', 'index', 'i_k', 'trmat', 'latt_valid',
               'trd_cell', 'WP_valid'])
    """
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(pseudo_result, 'lxml')
    dict = {}
    search_results_table = soup.find(lambda tag: tag.name == 'h3' and 'Select minimal supergroups of space group' in tag.text)

    table = search_results_table.find_next('table')
    for row in table.find_all('tr'):
        row_list = []
        for col in row.find_all('td'):
            row_list.append(col.text)
        if (len(row_list)):
            i = int(row_list[0])
            dict_aux = {}
            dict_aux['symb'] = row_list[2].strip()
            dict_aux['sgnum'] = row_list[3].strip()
            dict_aux['index'] = row_list[4].strip()
            dict_aux['i_k'] = row_list[5].strip()
            dict_aux['trmat'] = row_list[6]
            aux_trd_cell = row_list[7]
            dict_aux['latt_valid'] = True
            if(re.match(r'^.*Lattice parameters don\'t comply.*',aux_trd_cell)):
                aux_trd_cell = re.sub(r'Lattice parameters don\'t comply.*','',aux_trd_cell)
                dict_aux['latt_valid'] = False
            dict_aux['trd_cell'] = aux_trd_cell
            dict_aux['WP_valid'] = True
            if(re.search(r'is invalid under',row_list[8])):
                dict_aux['WP_valid'] = False
            dict[i] = dict_aux
    return dict

def parse_WYCKPOS_results(wyckpos_results):
    """
    Parses the result of a BCS_WYCKPOS() call
    and returns a dictionary with the following:
    'centering' : A list of the centering operations (empty if primitive)
    'wps': A dictionary with the following keys:
        'multiplicity'
        'letter'
        'site_symmetry'
        'coordinates'
    """
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(wyckpos_results, 'lxml')
    table = soup.find('table', {'border': '5'})
    dict_res = {'centering': []}
    dict_wps = {}
    for row in table.find_all('tr')[1:]:
        #print(row.text)
        ls = []
        cols = row.find_all('td')
        if len(cols) == 1 and cols[0].text.find('+') >= 0:
            #print(f"Centering: {cols[0].text}")
            centering = cols[0].text
            centering = re.sub(r'\+', '', centering)
            centering = re.split(r'\s+', centering.strip())
            dict_res['centering'] = centering
            continue
        for col in row.find_all('td'):
            #print(col.text,end=' | ')
            ls.append(col.text)
        if re.match('^[0-9]', ls[0]):
            wp = "".join(ls[:2])
            dict_aux = {'multiplicity': int(ls[0]), 'letter': ls[1],
                        'site_symmetry': ls[2],
                        'coordinates': re.split(r'\s+', ls[3].strip())}
            dict_wps[wp] = dict_aux
    dict_res['wps'] = dict_wps
    return dict_res

def print_WYCKPOS_results(wyckpos_results_dic):
    """
    Returns the result of a parse_WYCKPOS_results(BCS_WYCKPOS()) call
    as a neat string
    :param:
        wyckpos_results_dic: parsed BCS_WYCKPOS(sgnum) result (i.e., parse_WYCKPOS_results(wyckpos_results) )
    :return:
        prettified string version of the
        parse_WYCKPOS_results(wyckpos_results) return dictionary

    Example:
    result_raw = BCS_WYCKPOS(221)
    result_dic = parse_WYCKPOS_results(result_raw)
    print (print_WYCKPOS_results(result_dic))
    """
    dic = wyckpos_results_dic
    out_string = []
    out_string.append(f"Centering: {dic['centering']}")
    wp_keys = list(dic['wps'].keys())
    wp_keys = sorted(wp_keys, key=lambda x: (int(re.match(r'(\d+)', x).group(1)),
                                           re.search(r'[a-z]+', x).group()))
    out_string.append("Wyckoff Positions")
    for key in wp_keys:
        out_string.append(f"*{key}*")
        for prop, val in dic['wps'][key].items():
            out_string.append(f"\t{prop}: {val}")
    return "\n".join(out_string)
