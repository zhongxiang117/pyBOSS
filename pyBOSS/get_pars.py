"""Parse and get setting parameters

Priority:
    command line inputs > environmental variables > setting file > default configures

    specially, `PYBOSS_DATA_FILE_PATH` and `DATAFILE` in setting file can be
    both used for data files locating.
"""

from .__init__ import __version__
from .constants import PAR_SETTINGS, PARS_CMD, PYBOSS_DATA_FILE_PATH

import io
import os
import argparse
import configparser


# parameter groups for `argparse`
PAR_GROUPS = {
    'Calculation type setting': [
        'qmname', 'cmname', 'qmscale',
    ],
    'Factor to scale results': [
        'dielec', 'scllj', 'scl14c', 'scl14l',
    ],
    'External temperature and pressure': [
        't', 'p',
    ],
    'Solvent modes': [
        'svmod1', 'nmol', 'svmod2', 'nmol2', 'ibox', 'boxcut', 'rsolv', 
    ],
    'Solute setting': [
        'slfmt', 'ncent1', 'ncent2', 'izlong',
    ],
    'Polarization': [
        'polscx', 'polscs',
    ],
    'Cutoff settings': [
        'rcut', 'scut',
    ],
    'Calculation frequency': [
        'nvchg', 'nschg',
    ],
    'Perturbation settings': [
        'vdel', 'wkc', 'rdel', 'adel', 'rdels1', 'adels1', 'adels2', 'rdels2',
        'maxovl', 'irecnt', 
    ],
    'Charge settings': [
        'ich0', 'ich1', 'ich2',
    ],
    'Rdf settings': [
        'nrdf', 'nrdfs', 'nrdfs1', 'nrdfs2', 'nrdfa1', 'nrdfa2', 'nrdl',
    ],
    'Output settings': [
        'pltfmt', 
    ],
}

PAR_ENVS = [i.upper() for i in PARS_CMD.keys()]

# parameter types for `parse_pars`, all defined valid keys are listed in here
PAR_TYPES = {
    'int': [
        'ich0', 'ich1', 'ich2', 'icapat', 'nsaevl', 'nsatr', 'nmode', 'nsym', 'nconf', 'nmol',
        'ibox', 'nmol2', 'ncent1', 'ncent2', 'ncents', 'nrota1', 'nrota2', 'icut', 'irecnt',
        'indsol', 'izlong', 'maxovl', 'noxx', 'noss', 'nobndv', 'noangv', 'nonebn', 'nrdf',
        'nrdfs', 'nvchg', 'nschg', 'maxvar', 'nconsv', 'nbuse', 'nocut', 'nosmth', 'lhtsol',
        'iewald', 'isolec', 'ncons', 'icalc', 'inew', 'iprint', 'nrdl',
    ],
    'float': [
        'xstren', 'qmscale', 'caprad', 'capfk', 'ftol', 'satemp', 'sart', 'csetol', 'csrms',
        'csrtol', 'csreup', 'csaelo', 'csaeup', 'fscale', 'frqcut', 'boxcut', 'rdlmin', 'rdlinc',
        'eprmin', 'eprinc', 'edmin', 'edinc', 'essmin', 'essinc', 'ebsmin', 'ebsinc', 'vdel',
        'wkc', 'rdel', 'adel', 'rsolv', 'rdels1', 'adels1', 'rdels2', 'adels2', 'tlht', 'rcut',
        'scut', 'cutnb', 't', 'p', 'torcut', 'scllj', 'dielec', 'scl14c', 'scl14l', 'dielrf',
        'polscx', 'polscs', 'rc0', 'rc1', 'rc2'
    ],
    'iterint': ['icsopt', 'icutas', 'icutat', 'nrdfs1', 'nrdfs2', 'nrdfa1', 'nrdfa2'],
    'str': [
        'svmod1', 'qmname', 'cmname', 'slfmt', 'solor', 'slab', 'vxyz', 'optim', 'newzm',
        'svmod2', 'pltfmt', 'infile', 'upfile', 'save', 'average', 'slvzmat', 'bangpar',
        'waterbox', 'org1box', 'org2box', 'summary', 'output', 'pltfile', 'iflink', 'zmatrix',
        'pars', 'parameter'
    ],
}

def read_parameter_file(filename):
    """read parameter(PARAMETER) file
    
    Idea:
        To make life easier, line is right-padded to at least 80 chars.
        For each line, the read variable is appended a `_i` to indicate it should
        be an integer, `_f` for float number. And the last entry will always
        be the format and corresponding content.
        Array is specifically marked and dealt with.
    
    Rule:
        0) line are parsed to different variables in format: (var, value)
        1) check all parameters firstly, then exit if any error happens
        2) if value not exists:
            -> set to "" (empty string) for chars
            -> set to 0 (zero) for integer
            -> set to 0.0 (zero) for float
        3) if variable is required, add third item to True: (var, value, True)
        4) all items are in lower case
    """
    def _ls_process(ls,ev,dtype,errinfo):
        vs = []
        bo = True
        for i in ls:
            if i.strip() == '':
                vs.append(ev)
            else:
                try:
                    v = dtype(i)
                except ValueError:
                    print(f'Fatal: format: {errinfo[0]}:\n  -> {errinfo[1]}')
                    bo = False
                else:
                    vs.append(v)
        if bo:
            return vs
        return []

    if isinstance(filename,str):
        if os.path.isfile(filename):
            f = open(filename)
        else:
            print(f'Warning: not a file {filename}')
            return []
    elif hasattr(filename,'read'):      # for `io.StringIO()`
        f = filename
    else:
        print(f'Warning: unsupported streaming input: {filename}')
        return []

    mybool = True
    mypars = {}
    pars = []
    mypars['title'] = f.readline().strip()

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('svmod1', line[:5]),
        ('xstren_f', line[5:15]),
        ('A5,F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('qmname', line[:4]),
        ('cmname', line[4:8]),
        ('qmscale_f', line[8:13]),
        ('ich0_i', line[13:16]),
        ('ich1_i', line[16:19]),
        ('ich2_i', line[19:22]),
        ('2A4,F5.2,3I3', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('slfmt', line[:5]),
        ('A5', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('solor', line[:5]),
        ('slab', line[6:10]),
        ('vxyz', line[11:15]),
        ('A5,1X,A4,1X,A4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('icapat_i', line[:4]),
        ('caprad_f', line[4:14]),
        ('capfk_f', line[14:24]),
        ('I4,2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('optim', line[:5]),
        ('ftol_f', line[6:16]),
        ('newzm', line[16:21]),
        ('nsaevl_i', line[23:31]),
        ('nsatr_i', line[31:39]),
        ('satemp_f', line[39:49]),
        ('sart_f', line[49:59]),
        ('A5,1X,F10.6,A5,2X,2I8,2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        #('icsopt_i', line[:8]),        # 8 intergers
        ('csetol_f', line[8:18]),
        ('csrms_f', line[18:28]),
        ('csrtol_f', line[28:38]),
        ('csreup_f', line[38:48]),
        ('csaelo_f', line[48:58]),
        ('csaeup_f', line[58:68]),
        ('8I1,6F10.4', line.rstrip())
    ])
    ls = [i for i in line[:8]]
    vs = _ls_process(ls,0,int,pars[-1][-1])
    if vs:
        mypars['icsopt'] = vs
    else:
        mybool = False

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('nmode_i', line[:3]),
        ('nsym_i', line[3:6]),
        ('nconf_i', line[6:9]),
        ('fscale_f', line[9:19]),
        ('frqcut_f', line[19:29]),
        ('3I3,2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('nmol_i', line[:4]),
        ('ibox_i', line[4:8]),
        ('boxcut_f', line[8:14]),
        ('2I4,F6.2', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('svmod2', line[:5]),
        ('nmol2_i', line[5:9]),
        ('A5,I4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('ncent1_i', line[:3]),
        ('ncent2_i', line[3:6]),
        ('ncents_i', line[6:9]),
        #('icuta_i', line[9:12]),
        ('8I3', line.rstrip())
    ])
    ls = [line[i:i+3] for i in range(9,24,3)]
    vs = _ls_process(ls,0,int,pars[-1][-1])
    if vs:
        mypars['icutas'] = vs       # care: key change: icuta -> icutas
    else:
        mybool = False

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('nrota1_i', line[:3]),
        ('nrota2_i', line[3:6]),
        ('2I3', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('icut_i', line[:3]),
        #('icutat', line[3:48]),     # 15 integers
        ('16I3', line.rstrip())
    ])
    ls = [line[3+i:6+i] for i in range(0,45,3)]
    vs = _ls_process(ls,0,int,pars[-1][-1])
    if vs:
        mypars['icutat'] = vs
    else:
        mybool = False

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('irecnt_i', line[:3]),
        ('indsol_i', line[3:6]),
        ('izlong_i', line[6:9]),
        ('maxovl_i', line[9:12]),
        ('noxx_i', line[12:15]),
        ('noss_i', line[15:18]),
        ('nobndv_i', line[18:21]),
        ('noangv_i', line[21:24]),
        ('nonebn_i', line[24:27]),
        ('9I3', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('nrdf_i', line[:3]),
        ('nrdfs_i', line[3:6]),
        ('2I3', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    #pars.append([
    #    ('nrdfs1', line[:45]),
    #    ('15I3', line.rstrip())
    #])
    ls = [line[i:3+i] for i in range(0,45,3)]
    vs = _ls_process(ls,0,int,('15I3',line.rstrip()))
    if vs:
        mypars['nrdfs1'] = vs
    else:
        mybool = False

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    #pars.append([
    #    ('nrdfs2', line[:45]),
    #    ('15I3', line.rstrip())
    #])
    ls = [line[i:3+i] for i in range(0,45,3)]
    vs = _ls_process(ls,0,int,('15I3',line.rstrip()))
    if vs:
        mypars['nrdfs2'] = vs
    else:
        mybool = False

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    #pars.append([
    #    ('nrdfa1', line[:45]),
    #    ('15I3', line.rstrip())
    #])
    ls = [line[i:3+i] for i in range(0,45,3)]
    vs = _ls_process(ls,0,int,('15I3',line.rstrip()))
    if vs:
        mypars['nrdfa1'] = vs
    else:
        mybool = False

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    #pars.append([
    #    ('nrdfa2', line[:45]),
    #    ('15I3', line.rstrip())
    #])
    ls = [line[i:3+i] for i in range(0,45,3)]
    vs = _ls_process(ls,0,int,('15I3',line.rstrip()))
    if vs:
        mypars['nrdfa2'] = vs
    else:
        mybool = False

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('rdlmin_f', line[:10]),
        ('rdlinc_f', line[10:20]),
        ('2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('eprmin_f', line[:10]),
        ('eprinc_f', line[10:20]),
        ('2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('edmin_f', line[:10]),
        ('edinc_f', line[10:20]),
        ('2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('essmin_f', line[:10]),
        ('essinc_f', line[10:20]),
        ('2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('ebsmin_f', line[:10]),
        ('ebsinc_f', line[10:20]),
        ('2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('nvchg_i', line[:6]),
        ('nschg_i', line[6:12]),
        ('maxvar_i', line[12:18]),
        ('3I6', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('nconsv_i', line[:6]),
        ('nbuse_i', line[6:12]),
        ('nocut_i', line[12:18]),
        ('nosmth_i', line[18:24]),
        ('4I6', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('vdel_f', line[:10]),
        ('wkc_f', line[10:20]),
        ('2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('rdel_f', line[:10]),
        ('adel_f', line[10:20]),
        ('rsolv_f', line[20:30]),
        ('3F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('rdels1_f', line[:10]),
        ('adels1_f', line[10:20]),
        ('2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('rdels2_f', line[:10]),
        ('adels2_f', line[10:20]),
        ('2F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('lhtsol_i', line[:4]),
        ('tlht_f', line[4:14]),
        ('I4,F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('rcut_f', line[:10]),
        ('scut_f', line[10:20]),
        ('cutnb_f', line[20:30]),
        ('iewald_i', line[30:34]),
        ('3F10.4,I4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('t_f', line[:10]),
        ('p_f', line[10:20]),
        ('torcut_f', line[20:30]),
        ('scllj_f', line[30:40]),
        ('4F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('dielec_f', line[:10]),
        ('scl14c_f', line[10:20]),
        ('scl14l_f', line[20:30]),
        ('dielrf_f', line[30:40]),
        ('4F10.4', line.rstrip())
    ])

    dump = f.readline()
    line = f.readline().rstrip() + ' '*80
    pars.append([
        ('pltfmt', line[:5]),
        ('isolec_i', line[5:9]),
        ('polscx_f', line[9:20]),
        ('polscs_f', line[20:31]),
        ('A5,I4,2F11.4', line.rstrip())
    ])

    if isinstance(filename,str):
        f.close()
    for g in pars:
        for x in g[:-1]:      # always ignore the last
            k = x[0]
            v = x[1]
            bor = True if len(x)>2 else False
            if k.find('_i') != -1:
                nk = k.replace('_i','')
                nv = 0
                if v.strip() != '':
                    try:
                        nv = int(v)
                    except ValueError:
                        print(f'Fatal: format: {g[-1][0]}:\n  -> {g[-1][1]}')
                        mybool = False
                else:
                    if bor:
                        print(f'Fatal: format: {g[-1][0]}:\n  -> {g[-1][1]}')
                        mybool = False
            elif k.find('_f') != -1:
                nk = k.replace('_f','')
                nv = 0.0
                if v.strip() != '':
                    try:
                        nv = float(v)
                    except ValueError:
                        print(f'Fatal: format: {g[-1][0]}:\n  -> {g[-1][1]}')
                        mybool = False
                else:
                    if bor:
                        print(f'Fatal: format: {g[-1][0]}:\n  -> {g[-1][1]}')
                        mybool = False
            else:
                nk = k
                nv = v.strip().lower()      # important
                if bor and not nv:
                    print(f'Fatal: format: {g[-1][0]}:\n  -> {g[-1][1]}')
                    mybool = False
            mypars[nk] = nv
    mypars['bool'] = mybool
    return mypars


def write_settings_inner(envsettings={}, file='configs.txt'):
    cg = configparser.ConfigParser(allow_no_value=True)
    cg.add_section('root')

    c = 1
    for k,v in envsettings.items():
        cg.set('root', '# '+v[1], None)
        if isinstance(v[0],(tuple,list)):
            cg.set('root', k, ', '.join([str(i) for i in v[0]]))
        else:
            cg.set('root', k, str(v[0]))
        cg.set('root', '#@'+str(c), None)
        c += 1

    sio = io.StringIO()
    cg.write(sio)
    sio.seek(0)
    contents = sio.read()
    sio.close()
    contents = contents.replace('[root]\n','')
    for t in range(c-1,0,-1):
        contents = contents.replace('#@'+str(t),'')
    with open(file, 'wt') as f: f.write(contents)


def parse_pars(dpars):
    if not dpars: return {}
    pars = {}
    for k,v in dpars.items():
        if k in PAR_TYPES['int']:
            try:
                pars[k] = int(v)
            except ValueError:
                print(f'Fatal: wrong value type: (int): {v}')
                return {}
        elif k in PAR_TYPES['float']:
            try:
                pars[k] = float(v)
            except ValueError:
                print(f'Fatal: wrong value type: (float): {v}')
                return {}
        elif k in PAR_TYPES['iterint']:
            vl = v.replace(',', ' ').replace(';', ' ').replace('[', ' ').replace(']', ' ').split()
            ls = []
            for v in vl:
                try:
                    ls.append(int(v))
                except ValueError:
                    print(f'Fatal: wrong value type: (list-int, separated by comma/space): {v}')
                    return {}
            pars[k] = ls
        elif k in PAR_TYPES['str']:
            pars[k] = v
        else:
            print(f'Warning: unknown key pairs: {k} <-> {v}')
    return pars


def read_config_file(filename):
    if isinstance(filename,str):
        if os.path.isfile(filename):
            f = open(filename)
        else:
            print(f'Warning: not a file {filename}')
            return {}
    elif hasattr(filename,'read'):      # for `io.StringIO()`
        f = filename
    else:
        print(f'Warning: unsupported streaming input: {filename}')
        return []
    contents = '[root]\n' + f.read()
    f.close()

    cg = configparser.ConfigParser()
    cg.read_string(contents)
    return parse_pars({k:v for k,v in cg['root'].items()})


def get_env_settings():
    envs = {}
    for k in PARS_CMD:
        f = os.getenv(k.upper())
        if f:
            envs[k] = f
    return envs


def myargs():
    envs = ', '.join([i.upper() for i in PARS_CMD])
    parser = argparse.ArgumentParser(
        description='PyBOSS: Python Version of BOSS(from William L. Jorgensen) program',
        allow_abbrev=False,
        usage=argparse.SUPPRESS,
        epilog=f"""
        Input type: `f` for file, `s` for string, `i` for integer, `j` for float.
        Environmental Variables will also be used: {envs}. Priority: command line
        arguments > environmental variables > setting file > default configures
    """
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=__version__
    )

    usage = parser.add_argument_group('Multiple configurations')
    usage.add_argument(
        '-bc','--boss-parfile',
        metavar='f',
        type=str,
        dest='boss_parfile',
        help='BOSS `pmfpar` file, options inside it will be overwritten by command line inputs'
    )
    usage.add_argument(
        '-pyc','--py-config',
        metavar='f',
        type=str,
        dest='py_config',
        help='Python type configure file, lower priority than command line inputs and `--boss-config`'
    )
    usage.add_argument(
        '-path',
        metavar='dir',
        type=str,
        dest='path',
        help='command line option for `PYBOSS_DATA_FILE_PATH`'
    )
    usage.add_argument(
        '--output-full-settings',
        dest='bo_full',
        action='store_true',
        help='output full settings, if no inputs, a template file will be generated'
    )
    usage.add_argument(
        '--output-boss-parfile',
        dest='bo_bosspar',
        action='store_true',
        help='output collected parameters to BOSS type parfile'
    )
    usage.add_argument(
        '--no-output-input-settings',
        dest='no_bo_input',
        action='store_true',
        help='no output collected settings, from both command line and environmental variables'
    )
    usage.add_argument(
        '--debug',
        dest='debug',
        action='store_true',
        help='debug mode, more information will be printed out'
    )

    files = parser.add_argument_group('Options for inputs and outputs')
    for k,v in PARS_CMD.items():
        if isinstance(v[0],str):
            t = 'f'
            p = str
        elif isinstance(v[0],int):
            t = 'i'
            p = int
        else:
            t = 'j',
            p = float
        files.add_argument(
            '-'+k,
            metavar=t,
            type=p,
            help=v[1]
        )

    gsets = {i:j for i,j in PAR_SETTINGS.items()}
    for i,ks in PAR_GROUPS.items():
        g = parser.add_argument_group(i)
        for k in ks:
            v = gsets.pop(k)
            if isinstance(v[0],(tuple,list)):
                if isinstance(v[0][0],str):
                    t = 's'
                    p = str
                elif isinstance(v[0][0],int):
                    t = 'i'
                    p = int
                else:
                    t = 'j',
                    p = float
                g.add_argument(
                    '-'+k,
                    nargs=len(v[0]),
                    metavar=t,
                    type=p,
                    help=v[1]
                )
            else:
                if isinstance(v[0],str):
                    t = 's'
                    p = str
                elif isinstance(v[0],int):
                    t = 'i'
                    p = int
                else:
                    t = 'j',
                    p = float
                g.add_argument(
                    '-'+k,
                    metavar=t,
                    type=p,
                    help=v[1]
                )

    m = parser.add_argument_group('Miscellaneous, currently may not be supported')
    for k,v in gsets.items():
        if isinstance(v[0],(tuple,list)):
            if isinstance(v[0][0],str):
                t = 's'
                p = str
            elif isinstance(v[0][0],int):
                t = 'i'
                p = int
            else:
                t = 'j',
                p = float
            m.add_argument(
                '-'+k,
                nargs=len(v[0]),
                metavar=t,
                type=p,
                help=v[1]
            )
        else:
            if isinstance(v[0],str):
                t = 's'
                p = str
            elif isinstance(v[0],int):
                t = 'i'
                p = int
            else:
                t = 'j',
                p = float
            m.add_argument(
                '-'+k,
                metavar=t,
                type=p,
                help=v[1]
            )

    args = parser.parse_args()
    cpars = {k:v for k,v in vars(args).items() if v}

    # default configures
    paths = [p for p in PYBOSS_DATA_FILE_PATH if os.path.isdir(p)]
    ps = os.getenv('PYBOSS_DATA_FILE_PATH')
    if ps:
        ls = []
        for p in ps.split(':'):
            p = p.strip()
            if os.path.isdir(p):
                ls.append(p)
        paths.extend(ls[-1:0:-1])   # care, EVs priority are in reversed order
    paths.append('.')               # current path is appended
    if args.path and os.path.isdir(args.path):      # command line setting
        paths.append(args.path)

    if args.debug:
        print('Debug: PyBOSS searching data folder path: high -> low')
        for p in paths[-1:0:-1]:
            print(f'  -> {p}')

    # Priority: cmd-BOSS-pars > cmd-PyBOSS-pars > settingfile-BOSS-pars > settingfile-PyBOSS-pars
    bossp = None
    boboss = False
    if args.parameter:
        for p in paths[-1:0:-1]:
            np = os.path.join(p,args.parameter)
            if os.path.isfile(np):
                bossp = np
                boboss = True
                break
    elif args.pars:
        for p in paths[-1:0:-1]:
            np = os.path.join(p,args.pars)
            if os.path.isfile(np):
                bossp = np
                boboss = True
                break
    else:
        np = PARS_CMD['parameter'][0]
        for p in paths[-1:0:-1]:
            np = os.path.join(p,args.pars)
            if os.path.isfile(np):
                bossp = np
                break
        if not bossp:
            np = PARS_CMD['pars'][0]
            for p in paths[-1:0:-1]:
                np = os.path.join(p,args.pars)
                if os.path.isfile(np):
                    bossp = np
                    break
    if boboss: args.bo_full = True
    vpars = {}
    if bossp:
        print(f'Note: parameter file used (lowest priority): {bossp}')
        vpars = read_parameter_file(bossp)
        vpars['pars'] = bossp

    # PyBOSS configures
    ppars = {}
    if args.py_config:
        if os.path.isfile(args.py_config):
            print(f'Note: reading parameters from PyBOSS configures: {args.py_config}')
            ppars = read_config_file(args.py_config)
            ppars['pars'] = args.py_config
        else:
            print(f'Warning: PyBOSS configure: not a file: {args.py_config}')

    # BOSS configures
    bpars = {}
    if args.boss_parfile:
        args.bo_full = True
        print(f'Note: reading parameters from BOSS parfile: {args.boss_parfile}')
        bpars = read_parameter_file(args.boss_parfile)
        bpars['pars'] = args.boss_parfile

    # environmental variables
    epars = get_env_settings()

    # default configures
    dpars = {}
    for k,v in PAR_SETTINGS.items():
        dpars[k] = v[0]
    for k,v in PARS_CMD.items():
        dpars[k] = v[0]

    # now, do collection based on their priorities
    dpars.update(vpars)
    dpars.update(ppars)
    dpars.update(bpars)
    dpars.update(epars)
    #dpars.update(cpars)

    # `upars` used for user-setting-outputs, do not printout all settings
    if args.bo_full:
        upars = vpars
    else:
        upars = {}
    if args.py_config:
        upars.update(ppars)
    if args.boss_parfile:
        upars.update(bpars)
    upars.update(epars)
    #upars.update(cpars)

    # `cpars` are from command line, there are some keys we do not need
    allkeys = [*PAR_SETTINGS.keys(), *PARS_CMD.keys()]
    for k,v in cpars.items():
        if k in allkeys:
            dpars[k] = v
            upars[k] = v
    
    # add `PYBOSS_DATA_FILE_PATH`
    dpars['pyboss_data_file_path'] = paths

    if args.bo_full:
        # prepare for writing
        pdpars = {}
        for k,v in dpars.items():
            if k in PAR_SETTINGS:
                pdpars[k] = [v,PAR_SETTINGS[k][1]]
            elif k in PARS_CMD:
                pdpars[k] = [v,PARS_CMD[k][1]]
        write_settings_inner(pdpars)
    elif not args.no_bo_input:
        # prepare for writing
        pupars = {}
        for k,v in upars.items():
            if k in PAR_SETTINGS:
                pupars[k] = [v,PAR_SETTINGS[k][1]]
            elif k in PARS_CMD:
                pupars[k] = [v,PARS_CMD[k][1]]
        write_settings_inner(pupars)

    if args.bo_bosspar:
        write_boss_parfile(dpars)

    return dpars


def write_boss_parfile(dpars):
    """opposite operation with `read_parameter_file`"""
    pass


