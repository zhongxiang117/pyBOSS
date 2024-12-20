import os


def get_ff_bond_pars(data,msg=True):
    """
    Format: (I4,1X,I2,1X,A2,1X,3F10.6)
    Allow errors exist
    """
    pars = []
    for line in data:
        line = line.rstrip()
        if not line or line[0] == '#': continue
        line += ' '*80
        try:
            t = int(line[:4])
            s = int(line[5:7])
            a = line[8:10].strip()
            chg = float(line[11:21])
            sig = float(line[21:31])
            eps = float(line[31:41])
        except ValueError:
            if msg:
                print(f'Warning: bond parameters: {line.rstrip()}')
        else:
            pars.append((t,s,a,chg,sig,eps))
    return pars


def get_ff_dihedral_pars(data,msg=True):
    """
    Format: (I3,2X,4(F8.4,1X,I1),2X,4(A2,1X))
    Allow errors exist
    """
    pars = []
    for line in data:
        line = line.rstrip()
        if not line or line[0] == '#': continue
        line += ' '*80
        try:
            t = int(line[:3])
            v1 = float(line[5:13])
            v1i = int(line[14]) if line[14] != ' ' else 0
            v2 = float(line[15:23])
            v2i = int(line[24]) if line[24] != ' ' else 0
            v3 = float(line[25:33])
            v3i = int(line[34]) if line[34] != ' ' else 0
            v4 = float(line[35:43])
            #v4i = int(line[44]) if line[44] != ' ' else 0
            v4i = 0         # place-holder
            s = line[47:58]
            p = [s[0:2], s[3:5], s[6:8], s[9:11]]
            for i in p:
                if i[0] == ' ' or i == '  ':
                    raise ValueError
        except ValueError:
            if msg:
                print(f'Warning: dihedral parameters: {line.rstrip()}')
        else:
            pars.append((t,v1,v2,v3,v4,v1i,v2i,v3i,v4i,*p))
    return pars


def _separator(prolines):
    before = after = len(prolines)
    bo = True
    for i,l in enumerate(prolines):
        if bo:
            if len(l) < 4 or l[:4] == '    ':
                bo = False
                before = i
        else:
            if not (len(l) < 4 or l[:4] == '    '):
                after = i
                break
    return (before,after)


def read_ff_file(filename,skip=None,msg=True):
    """read force field parameter file

    Return:
        bond_pars: 
            [(index, atom_type_num, atom_type_name, charge, sigma, epsilon),   ]
        
        dihedral_pars:
            [(index, v1, v2, v3, v4, I1, I2, I3, I4, type),   ]
    """
    if not os.path.isfile(filename):
        print(f'Warning: not a file {filename}')
        return (tuple(),tuple())
    prolines = []
    with open(filename,'rt') as f:
        if skip:
            for i in range(skip):
                f.readline()
        prolines = f.readlines()        # newline mark '\n' is in end
    before,after = _separator(prolines)
    bond_pars = get_ff_bond_pars(prolines[:before],msg)
    dihedral_pars = get_ff_dihedral_pars(prolines[after:],msg)
    return (bond_pars, dihedral_pars)


def get_oplsaa_bond_pars(data,msg=True):
    """
    Format: (A2,1X,A2,2F10.4)
    Allow errors exist
    """
    pars = []
    for line in data:
        line = line.rstrip()
        if not line or line[0] == '#': continue
        line += ' '*80
        try:
            v = float(line[5:15])
            t = float(line[15:25])
        except ValueError:
            if msg:
                print(f'Warning: oplsaa bond parameters: {line.rstrip()}')
        else:
            a = line[0:2]
            b = line[3:5]
            pars.append((a,b,v,t))
    return pars


def get_oplsaa_angle_pars(data,msg=True):
    """
    Format: (A2,1X,A2,1X,A2,2X,2F10.3)
    Allow errors exist
    """
    pars = []
    for line in data:
        line = line.rstrip()
        if not line or line[0] == '#': continue
        line += ' '*80
        try:
            v = float(line[10:20])
            t = float(line[20:30])
        except ValueError:
            if msg:
                print(f'Warning: oplsaa angle parameters: {line.rstrip()}')
        else:
            a = line[0:2]
            b = line[3:5]
            c = line[6:8]
            pars.append((a,b,c,v,t))
    return pars


def read_oplsaa_file(filename,skip=None,msg=True):
    """read oplsaa parameter file

    Return:
        oplsaa_bond_pars: 
            [(index, atom_type_num, atom_type_name, charge, sigma, epsilon),   ]
        
        oplsaa_angle_pars:
            [(index, v1, v2, v3, v4, I1, I2, I3, I4, type),   ]
    """
    if not os.path.isfile(filename):
        print(f'Warning: not a file {filename}')
        return (tuple(),tuple())
    prolines = []
    with open(filename,'rt') as f:
        if skip:
            for i in range(skip):
                f.readline()
        prolines = f.readlines()        # newline mark '\n' is in end
    before,after = _separator(prolines)
    bond_pars = get_oplsaa_bond_pars(prolines[:before],msg)
    angle_pars = get_oplsaa_angle_pars(prolines[after:],msg)
    return (bond_pars, angle_pars)


