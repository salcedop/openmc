import argparse
from collections import namedtuple
from io import StringIO
import os
import shutil
from subprocess import Popen, PIPE, STDOUT, CalledProcessError
import sys
import tempfile

from . import endf


# For a given MAT number, give a name for the ACE table and a list of ZAID
# identifiers
ThermalTuple = namedtuple('ThermalTuple', ['name', 'zaids', 'nmix'])
_THERMAL_DATA = {
    1: ThermalTuple('hh2o', [1001], 1),
    2: ThermalTuple('parah', [1001], 1),
    3: ThermalTuple('orthoh', [1001], 1),
    5: ThermalTuple('hyh2', [1001], 1),
    7: ThermalTuple('hzrh', [1001], 1),
    8: ThermalTuple('hcah2', [1001], 1),
    10: ThermalTuple('hice', [1001], 1),
    11: ThermalTuple('dd2o', [1002], 1),
    12: ThermalTuple('parad', [1002], 1),
    13: ThermalTuple('orthod', [1002], 1),
    26: ThermalTuple('be', [4009], 1),
    27: ThermalTuple('bebeo', [4009], 1),
    31: ThermalTuple('graph', [6000, 6012, 6013], 1),
    33: ThermalTuple('lch4', [1001], 1),
    34: ThermalTuple('sch4', [1001], 1),
    37: ThermalTuple('hch2', [1001], 1),
    39: ThermalTuple('lucite', [1001], 1),
    40: ThermalTuple('benz', [1001, 6000, 6012], 2),
    41: ThermalTuple('od2o', [8016, 8017, 8018], 1),
    43: ThermalTuple('sisic', [14028, 14029, 14030], 1),
    44: ThermalTuple('csic', [6000, 6012, 6013], 1),
    46: ThermalTuple('obeo', [8016, 8017, 8018], 1),
    47: ThermalTuple('sio2-a', [8016, 8017, 8018, 14028, 14029, 14030], 3),
    48: ThermalTuple('uuo2', [92238], 1),
    49: ThermalTuple('sio2-b', [8016, 8017, 8018, 14028, 14029, 14030], 3),
    50: ThermalTuple('oice', [8016, 8017, 8018], 1),
    52: ThermalTuple('mg24', [12024], 1),
    53: ThermalTuple('al27', [13027], 1),
    55: ThermalTuple('yyh2', [39089], 1),
    56: ThermalTuple('fe56', [26056], 1),
    58: ThermalTuple('zrzrh', [40000, 40090, 40091, 40092, 40094, 40096], 1),
    59: ThermalTuple('cacah2', [20040, 20042, 20043, 20044, 20046, 20048], 1),
    75: ThermalTuple('ouo2', [8016, 8017, 8018], 1),
}

_TEMPLATE_MODER = """
moder / %%%%%%%%%%%%%%%%%%% Extract XS data%%%%%%%%%%%%%%%%%%%%%%%
1 {nmoder}
'{library} PENDF for {zsymam}'/
{nendf} {mat}
0 /
"""
_TEMPLATE_MODER_OUT = """
moder / %%%%%%%%%%%%%%%%%%% Extract XS data%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {noutmoder} /
0 /
"""
_TEMPLATE_RECONR = """
reconr / %%%%%%%%%%%%%%%%%%% Reconstruct XS for neutrons %%%%%%%%%%%%%%%%%%%%%%%
{nreconr_in} {nreconr}
'{library} PENDF for {zsymam}'/
{mat} 2/
{error}/ err
'{library}: {zsymam}'/
'Processed by NJOY'/
0/
"""

_TEMPLATE_BROADR = """
broadr / %%%%%%%%%%%%%%%%%%%%%%% Doppler broaden XS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {broadr_in} {nbroadr}
{mat} {num_temp} 0 0 0. /
{error} 2.e+6/ errthn
{temps}
0/
"""

_TEMPLATE_HEATR = """
heatr / %%%%%%%%%%%%%%%%%%%%%%%%% Add heating kerma %%%%%%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {nheatr_in} {nheatr} /
{mat} 3 /
302 318 402 /
"""

_TEMPLATE_UNRESR = """
unresr / %%%%%%%%%%%%%%%%%%%%%%%% Add probability tables %%%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {nunresr_in} {nunresr} /
{mat} {num_temp} {n_sigz} 1 /
{temps}
{sigz}
0/
"""

_TEMPLATE_PURR = """
purr / %%%%%%%%%%%%%%%%%%%%%%%% Add probability tables %%%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {npurr_in} {npurr} /
{mat} {num_temp} {n_sigz} 20 64 /
{temps}
{sigz}
0/
"""

_TEMPLATE_GROUPR = """
groupr/ %%%%%%%%%%%%%%%%%%%%%%%% Generate Multi-Group cross-section %%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {nunresr} 0 {ngroupr} /
{mat} 1 0 2 {lorde} {num_temp} {n_sigz} 1 /
'Running groupr for {mat}' /
{temps} /
{sigz} /
{num_groups} /
{group_struct} /
{end_lines}
"""

_TEMPLATE_GROUPR_BUILTIN_STRUCT = """
groupr/ %%%%%%%%%%%%%%%%%%%%%%%% Generate Multi-Group cross-section %%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {nunresr} 0 {ngroupr} /
{mat} 12 0 3 {lorde} {num_temp} {n_sigz} 1 /
'Running groupr for {mat}' /
{temps} /
{sigz} /
{end_lines}
"""

_TEMPLATE_GROUPR_ONE_OVER_E = """
groupr/ %%%%%%%%%%%%%%%%%%%%%%%% Generate Multi-Group cross-section %%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {nunresr} 0 {ngroupr} /
{mat} 5 0 4 {lorde} {num_temp} {n_sigz} 1 /
'Running groupr for {mat}' /
{temps} /
{sigz} /
0.1 0.025 0.82086e06 1.4e06 /
{end_lines}
"""

_TEMPLATE_ACER = """
acer / %%%%%%%%%%%%%%%%%%%%%%%% Write out in ACE format %%%%%%%%%%%%%%%%%%%%%%%%
{ninput} {nacer_in} 0 {nace} {ndir}
1 0 1 .{ext} /
'{library}: {zsymam} at {temperature}'/
{mat} {temperature}
1 1/
/
"""

_THERMAL_TEMPLATE_THERMR = """
thermr / %%%%%%%%%%%%%%%% Add thermal scattering data (free gas) %%%%%%%%%%%%%%%
0 {nthermr1_in} {nthermr1}
0 {mat} 12 {num_temp} 1 0 {iform} 1 221 1/
{temps}
{error} {energy_max}
thermr / %%%%%%%%%%%%%%%% Add thermal scattering data (bound) %%%%%%%%%%%%%%%%%%
{nthermal_endf} {nthermr2_in} {nthermr2}
{mat_thermal} {mat} 16 {num_temp} {inelastic} {elastic} {iform} {natom} 222 1/
{temps}
{error} {energy_max}
"""

_THERMAL_TEMPLATE_ACER = """
acer / %%%%%%%%%%%%%%%%%%%%%%%% Write out in ACE format %%%%%%%%%%%%%%%%%%%%%%%%
{nendf} {nthermal_acer_in} 0 {nace} {ndir}
2 0 1 .{ext}/
'{library}: {zsymam_thermal} processed by NJOY'/
{mat} {temperature} '{data.name}' /
{zaids} /
222 64 {mt_elastic} {elastic_type} {data.nmix} {energy_max} 2/
"""

def run(commands, tapein, tapeout, num_groups,input_filename=None,stdout=False,
        njoy_exec='njoy',groupr=False):
  """Run NJOY with given commands

    Parameters
    ----------
    commands : str
        Input commands for NJOY
    tapein : dict
        Dictionary mapping tape numbers to paths for any input files
    tapeout : dict
        Dictionary mapping tape numbers to paths for any output files
    input_filename : str, optional
        File name to write out NJOY input commands
    stdout : bool, optional
        Whether to display output when running NJOY
    njoy_exec : str, optional
        Path to NJOY executable

    Raises
    ------
    subprocess.CalledProcessError
        If the NJOY process returns with a non-zero status

  """

  if input_filename is not None:
      with open(input_filename, 'w') as f:
        f.write(commands)

    
  with tempfile.TemporaryDirectory() as tmpdir:
      
        
        # Copy evaluations to appropriates 'tapes'
    for tape_num, filename in tapein.items():
            tmpfilename = os.path.join(tmpdir, 'tape{}'.format(tape_num))
            shutil.copy(filename, tmpfilename)

        # Start up NJOY process
    njoy = Popen([njoy_exec], cwd=tmpdir, stdin=PIPE, stdout=PIPE,
                     stderr=STDOUT, universal_newlines=True)

    njoy.stdin.write(commands)
    njoy.stdin.flush()
    lines = []
    while True:
            # If process is finished, break loop
            line = njoy.stdout.readline()
            if not line and njoy.poll() is not None:
                break

            lines.append(line)
            if stdout:
                # If user requested output, print to screen
                
                print(line, end='')

        # Check for error
    if njoy.returncode != 0:
            raise CalledProcessError(njoy.returncode, njoy_exec,
                                     ''.join(lines))

        # Copy output files back to original directory
        
    for tape_num, filename in tapeout.items():
            tmpfilename = os.path.join(tmpdir, 'tape{}'.format(tape_num))
            if os.path.isfile(tmpfilename):
                print(filename)
                shutil.move(tmpfilename, filename)
    if (groupr):

       filename = 'njoy_output_groupr_{}'.format(num_groups)
    else:

       filename = 'njoy_output'
    outputfile = os.path.join(tmpdir,'output')
    shutil.move(outputfile,filename)
    
def make_pendf(filename, pendf='pendf', error=0.001, stdout=False):
    """Generate ACE file from an ENDF file

    Parameters
    ----------
    filename : str
        Path to ENDF file
    pendf : str, optional
        Path of pointwise ENDF file to write
    error : float, optional
        Fractional error tolerance for NJOY processing
    stdout : bool
        Whether to display NJOY standard output

    Raises
    ------
    subprocess.CalledProcessError
        If the NJOY process returns with a non-zero status

    """

    make_ace(filename, pendf=pendf, error=error, broadr=False,
             heatr=False, purr=False, acer=False, stdout=stdout)


def make_xs(filename,lorde=None, group_struct=None,num_groups=None,temperatures=None, background_xs=None,ace='ace', xsdir='xsdir', pendf=None,error=0.001, moder=True, broadr=True, heatr=True, unresr=False,purr=True, groupr=False,acer=True,one_over_e = False,
             stdout=False,**kwargs):
    """Generate incident neutron ACE file from an ENDF file

    Parameters
    ----------
    filename : str
        Path to ENDF file
    temperatures : iterable of float, optional
        Temperatures in Kelvin to produce ACE files at. If omitted, data is
        produced at room temperature (293.6 K).
    background_xs : iterable of float, optional
        Background cross-sections in barns, needed when taking into account
        self-shielding effects in unresolved resonance-range.
        if omitted, infinite dilution (i.e., 1E.10) is selected by default.
    group_struct : iterable of float, optional
        arbitrary group structure for GROUPR module. 
    num_groups : integer, bool
        number of groups in the fine-structure. 
        Important to note that this is minus one the size of your group_struct array. 
        defaults to 1.
    ace : str, optional
        Path of ACE file to write
    xsdir : str, optional
        Path of xsdir file to write
    pendf : str, optional
        Path of pendf file to write. If omitted, the pendf file is not saved.
    error : float, optional
        Fractional error tolerance for NJOY processing
    broadr : bool, optional
        Indicating whether to Doppler broaden XS when running NJOY
    unresr : bool, optional
        Indicating whether to take into account self-shielding effects in
        the unresolved resonance-range. 
    moder : bool, optional
        Allows for conversion of formatted mode to binary files and vice-versa.
    groupr : bool, optional
        Generates multi-group cross-sections. This module provides lots of flexi-
        bility to the user but here we limit ourselved to just using a couple of 
        weight-functions. The user is also forced to provide an arbitrary group 
        structure.
    lorde : integer, bool
        Legendre polynomial order of the angular flux. This is needed when computing the Bondarenko 
        model inside GROUPR. Possible values are in the range from 1-6. 
        optional, defaults to 1.
    one_over_e : bool, optional
        Thus far, when running GROUPR we have to provide an arbitrary group structure.
        Howerer, NJOY also needs a weight-function to compute the group cross-sections.
        If this variable is true then NJOY will use the following weight-function: 1/E + maxwellian +        fission spectrum. By default, a flat weight-function is selected. This option works well for         fine-group structure like the ones that will be needed when doing flux-tallies. 
    heatr : bool, optional
        Indicating whether to add heating kerma when running NJOY
    purr : bool, optional
        Indicating whether to add probability table when running NJOY
    acer : bool, optional
        Indicating whether to generate ACE file when running NJOY
  
    **kwargs
        Keyword arguments passed to :func:`openmc.data.njoy.run`

    Raises
    ------
    subprocess.CalledProcessError
        If the NJOY process returns with a non-zero status

    """
    ev = endf.Evaluation(filename)
    mat = ev.material
    zsymam = ev.target['zsymam']

    # Determine name of library
    library = '{}-{}.{}'.format(*ev.info['library'])

    if temperatures is None:
        temperatures = [293.6]
    num_temp = len(temperatures)
    temps = ' '.join(str(i) for i in temperatures)
    
    if ((background_xs is None) or (groupr)):
    #defaults to infinite dilution.
        background_xs = [1.e10]
    n_sigz = len(background_xs)
    sigz = ' '.join(str(i) for i in background_xs)
    
    if lorde is None:
        lorde = 1
 
    # Create njoy commands by modules
    commands = ""

    nendf, npendf = 20, 21
    tapein = {nendf: filename}
    tapeout = {}
    if pendf is not None:
        tapeout[npendf] = pendf
    ninput = nendf
    nlast = npendf - 1
    #mult can be either 1 or -1
    #depending on whether or not 
    #we call MODER. The thing is
    #MODER output files by convention have a minus
    #sign so all the modules call afterward need to
    #respond to this appropriately. MODER needs to be used
    #when planning to call GROUPR. In the case of ACER,
    #it is not necessary. 
    mult = 1
    # moder
    if moder:
       mult = -1
       ninput = mult * npendf
       nmoder = ninput
       commands += _TEMPLATE_MODER
       nlast = (nlast + 1) * mult
      
    # reconr

    nreconr_in = nlast
    nreconr = nreconr_in + mult
    commands += _TEMPLATE_RECONR
    nlast = nreconr

    # broadr
    if broadr:
        broadr_in = nlast
        nbroadr = broadr_in + mult
        commands += _TEMPLATE_BROADR
        nlast = nbroadr

    # heatr
    if heatr:
       nheatr_in = nlast
       nheatr = nheatr_in + mult
       commands += _TEMPLATE_HEATR
       nlast = nheatr

    #groupr

    if groupr:
       
       #unresr logic
       nunresr_in = nlast
       nunresr = nunresr_in + mult
       commands += _TEMPLATE_UNRESR
       nlast = nunresr
       #ngroupr_in = nunresr
       commands = commands.format(**locals())  
       ngroupr_in = nlast
       ngroupr = ngroupr_in + mult
       nlast = ngroupr
       fname = '{}' 
       #'card' #9 of the GROUPR module specifies whether or not
       #the cross-section should be in vector or matrix format.
       #In our case, the the multi-group cross-sections will
       #always be in vector format which is denoted by a '3 /'. 
       #the mt values of all the reactions of interests
       #need to be specified right after (see Pages 232 and 244 from the NJOY2016 manual.
       #Ideally, one would only need to specify the mt values correponding to
       #the 7 depletion reaction rates. However, not all nuclides
       #have all 7 and GROUPR will crash if it finds a nuclide does not
       #have one of the specified reactions. To get around this, we just tell
       #NJOY to compute multi-group cross-sections for all the available
       #reactions of each particular nuclide. This is done by specifying
       #a '0 /' after the '3 /'. However, if the user inputs multiple temperatures or materials,
       #then a '3 /\n0 /\n' pair must be added for each temperature..
       end_lines = ''
       for i in range(0,num_temp):
          end_lines += '3 /\n0 /\n'   
       end_lines += '0 /\n'
       if group_struct is None:
            #determine whether or not to use 1/E+fission-spectrum+maxwellian as
            #weight-function instead of a flat weight function.
            #the flat option is recommended by the NJOY manual 
            #when using a fine-group structure. 
            if (one_over_e):
              print("GROUPR module will use the following weight function: 1/E + maxwelliam + fission spectrum, with rrd 50-group structure")
              num_groups = 50
              commands += (_TEMPLATE_GROUPR_ONE_OVER_E).format(**locals())
            else:
              print("GROUPR module was specified without arbitrary structure, rr 50-group structure will be use instead")
              num_groups = 620
              commands += (_TEMPLATE_GROUPR_BUILTIN_STRUCT).format(**locals())
       else:
              commands += (_TEMPLATE_GROUPR).format(**locals())
       
       commands += 'stop\n'
       print(commands)
       run(commands,tapein, tapeout,num_groups=num_groups,stdout=stdout,groupr=groupr,**kwargs)

    #make sure to consider this logic only if GROUPR wasn't requested.
    #might be worth it to have a separate function called make_groupr
    else: 
     # purr
     if purr:
        npurr_in = nlast
        npurr = npurr_in + mult
        commands += _TEMPLATE_PURR
        nlast = npurr * mult


     commands = commands.format(**locals())

     # acer
     if acer:
        num_groups = None
        nacer_in = nlast * mult
        fname = '{}_{:.1f}'
        for i, temperature in enumerate(temperatures):
            # Extend input with an ACER run for each temperature
            nace = nacer_in * mult + 1 + 2*i
            ndir = nace + 1
            ext = '{:02}'.format(i + 1)
            commands += _TEMPLATE_ACER.format(**locals())

            # Indicate tapes to save for each ACER run
            tapeout[nace] = fname.format(ace, temperature)
            tapeout[ndir] = fname.format(xsdir, temperature)
 
     commands += 'stop\n'
     print(commands)
     run(commands,tapein, tapeout,num_groups=num_groups,stdout=stdout,groupr=groupr,**kwargs)

     if acer:
        with open(ace, 'w') as ace_file, open(xsdir, 'w') as xsdir_file:
            for temperature in temperatures:
                # Get contents of ACE file
                text = open(fname.format(ace, temperature), 'r').read()

                # If the target is metastable, make sure that ZAID in the ACE file reflects
                # this by adding 400
                if ev.target['isomeric_state'] > 0:
                    mass_first_digit = int(text[3])
                    if mass_first_digit <= 2:
                        text = text[:3] + str(mass_first_digit + 4) + text[4:]

                # Concatenate into destination ACE file
                ace_file.write(text)

                # Concatenate into destination xsdir file
                text = open(fname.format(xsdir, temperature), 'r').read()
                xsdir_file.write(text)

        # Remove ACE/xsdir files for each temperature
        
        for temperature in temperatures:
            os.remove(fname.format(ace, temperature))
            os.remove(fname.format(xsdir, temperature))

        
def make_ace_thermal(filename, filename_thermal, temperatures=None,
                     ace='ace', xsdir='xsdir', error=0.001, **kwargs):
    """Generate thermal scattering ACE file from ENDF files

    Parameters
    ----------
    filename : str
        Path to ENDF neutron sublibrary file
    filename_thermal : str
        Path to ENDF thermal scattering sublibrary file
    temperatures : iterable of float, optional
        Temperatures in Kelvin to produce data at. If omitted, data is produced
        at all temperatures given in the ENDF thermal scattering sublibrary.
    ace : str, optional
        Path of ACE file to write
    xsdir : str, optional
        Path of xsdir file to write
    error : float, optional
        Fractional error tolerance for NJOY processing
    **kwargs
        Keyword arguments passed to :func:`openmc.data.njoy.run`

    Raises
    ------
    subprocess.CalledProcessError
        If the NJOY process returns with a non-zero status

    """
    ev = endf.Evaluation(filename)
    mat = ev.material
    zsymam = ev.target['zsymam']

    ev_thermal = endf.Evaluation(filename_thermal)
    mat_thermal = ev_thermal.material
    zsymam_thermal = ev_thermal.target['zsymam']

    data = _THERMAL_DATA[mat_thermal]
    zaids = ' '.join(str(zaid) for zaid in data.zaids[:3])

    # Determine name of library
    library = '{}-{}.{}'.format(*ev_thermal.info['library'])

    # Determine if thermal elastic is present
    if (7, 2) in ev_thermal.section:
        elastic = 1
        mt_elastic = 223

        # Determine whether elastic is incoherent (0) or coherent (1)
        file_obj = StringIO(ev_thermal.section[7, 2])
        elastic_type = endf.get_head_record(file_obj)[2] - 1
    else:
        elastic = 0
        mt_elastic = 0
        elastic_type = 0

    # Determine number of principal atoms
    file_obj = StringIO(ev_thermal.section[7, 4])
    items = endf.get_head_record(file_obj)
    items, values = endf.get_list_record(file_obj)
    energy_max = values[3]
    natom = int(values[5])

    # Note that the 'iform' parameter is omitted in NJOY 99. We assume that the
    # user is using NJOY 2012 or later.
    iform = 0
    inelastic = 2

    # Determine temperatures from MF=7, MT=4 if none were specified
    if temperatures is None:
        file_obj = StringIO(ev_thermal.section[7, 4])
        endf.get_head_record(file_obj)
        endf.get_list_record(file_obj)
        endf.get_tab2_record(file_obj)
        params = endf.get_tab1_record(file_obj)[0]
        temperatures = [params[0]]
        for i in range(params[2]):
            temperatures.append(endf.get_list_record(file_obj)[0][0])

    num_temp = len(temperatures)
    temps = ' '.join(str(i) for i in temperatures)

    # Create njoy commands by modules
    commands = ""

    nendf, nthermal_endf, npendf = 20, 21, 22
    tapein = {nendf: filename, nthermal_endf:filename_thermal}
    tapeout = {}

    # reconr
    commands += _TEMPLATE_RECONR
    nlast = npendf

    # broadr
    nbroadr = nlast + 1
    commands += _TEMPLATE_BROADR
    nlast = nbroadr

    # thermr
    nthermr1_in = nlast
    nthermr1 = nthermr1_in + 1
    nthermr2_in = nthermr1
    nthermr2 = nthermr2_in + 1
    commands += _THERMAL_TEMPLATE_THERMR
    nlast = nthermr2

    commands = commands.format(**locals())

    # acer
    nthermal_acer_in = nlast
    fname = '{}_{:.1f}'
    for i, temperature in enumerate(temperatures):
        # Extend input with an ACER run for each temperature
        nace = nthermal_acer_in + 1 + 2*i
        ndir = nace + 1
        ext = '{:02}'.format(i + 1)
        commands += _THERMAL_TEMPLATE_ACER.format(**locals())

        # Indicate tapes to save for each ACER run
        tapeout[nace] = fname.format(ace, temperature)
        tapeout[ndir] = fname.format(xsdir, temperature)
    commands += 'stop\n'
    run(commands, tapein, tapeout, **kwargs)

    with open(ace, 'w') as ace_file, open(xsdir, 'w') as xsdir_file:
        # Concatenate ACE and xsdir files together
        for temperature in temperatures:
            text = open(fname.format(ace, temperature), 'r').read()
            ace_file.write(text)

            text = open(fname.format(xsdir, temperature), 'r').read()
            xsdir_file.write(text)

    # Remove ACE/xsdir files for each temperature
    for temperature in temperatures:
        os.remove(fname.format(ace, temperature))
        os.remove(fname.format(xsdir, temperature))
