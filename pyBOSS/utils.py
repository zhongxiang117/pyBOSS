import os
import math
import numpy as np


def read_zmat(filename):
    """read zmatrix
    
    keys:
        title:
        number_of_entries:
        entry_(i)    : until to number of entries(included), start from 1
            entry format: (atom_type,initial,final,b,bv,a,av,d,dv,resid,resname)
    """
    if not os.path.isfile(filename):
        print(f'Warning: not a file {filename}')
        return []
    zmat = {'number_of_entries':0, 'title':'',}
    alldata = []
    with open(filename,'rt') as f:
        zmat['title'] = f.readline().strip()
        cnt = 1
        data = []
        for line in f:
            line = line.rstrip() + ' '*80
            if line[:4] == '    ':          # means end reading
                break
            if line[:3].lower() == 'ter':   # reset
                zmat['number_of_atoms_for_entry_'+str(cnt)] = len(data)
                cnt += 1
                alldata.extend(data)
                data = []
            atype = line[5:8].strip()
            try:
                beg = int(line[9:13])
                end = int(line[14:18])
                bn = int(line[19:23])
                bv = float(line[23:35])
                an = int(line[35:39])
                av = float(line[39:51])
                dn = int(line[51:55])
                dv = float(line[55:67])
                sid = line[72:76].strip()
                if sid:
                    resid = int(line[72:76])
                else:
                    resid = 0
            except ValueError:
                print(f'Fatal: wrong defined line:\n  -> {line.rstrip()}')
                return {}
            resname = line[68:71].strip()
            data.append((atype,beg,end,bn,bv,an,av,dn,dv,resid,resname))
        if data:
            zmat['number_of_atoms_for_entry_'+str(cnt)] = len(data)
            alldata.extend(data)
        else:
            cnt -= 1
        zmat['number_of_entries'] = cnt
        zmat['data'] = alldata
        num_atoms = len(zmat['data'])

        # now read geometry bonds
        bonds_geometry = []
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip() + ' '*80
                if line[:4] == '    ': break
                try:
                    b = int(line[0:4])
                    e = int(line[4:8])
                    v = float(line[8:20])
                    if e not in [1,2,3]:
                        raise ValueError
                except ValueError:
                    print(f'Fatal: wrong line: geometry bonds\n  -> {line.rstrip()}')
                    return {}
                # check
                if b > num_atoms or b <= 0:
                    print(f'Fatal: geometry bonds: wrong value:\n  -> {line.rstrip()}')
                    return {}
                bonds_geometry.append((b-1,e,v))
        zmat['bonds_geometry'] = tuple(set(bonds_geometry))

        # now read variable bonds
        bonds_variable = []
        bonds_equal = []
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip() + ' '*80
                if line[:4] == '    ': break
                if line[4] not in [' ','-','=']:
                    print(f'Fatal: wrong line: variable bonds\n  -> {line.rstrip()}')
                    return {}
                try:
                    if line[5:9] == '    ':
                        b = 1
                        e = int(line[0:4])
                    else:
                        b = int(line[0:4])
                        e = int(line[5:9])
                except ValueError:
                    print(f'Fatal: wrong line: variable bonds\n  -> {line.rstrip()}')
                    return {}
                # check
                if b > num_atoms or e > num_atoms or b <= 0 or e <= 0:
                    print(f'Fatal: additional bonds: wrong value:\n  -> {line.rstrip()}')
                    return {}
                if line[4] == '=':
                    bonds_equal.append((b-1,e-1))
                else:
                    if b > e: b, e = e, b
                    bonds_variable.extend(list(range(b-1,e)))
        zmat['bonds_variable'] = tuple(set(bonds_variable))
        zmat['bonds_equal'] = tuple(set(bonds_equal))

        # now read additional bonds
        bonds_additional = []
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip() + ' '*80
                if line[:4] == '    ': break
                try:
                    b = int(line[0:4])
                    e = int(line[4:8])
                except ValueError:
                    print(f'Fatal: wrong line: additional bonds\n  -> {line.rstrip()}')
                    return {}
                # check
                if b > num_atoms or e > num_atoms or b <= 0 or e <= 0:
                    print(f'Fatal: additional bonds: wrong value:\n  -> {line.rstrip()}')
                    return {}
                bonds_additional.append((b-1,e-1))
        zmat['bonds_additional'] = tuple(set(bonds_additional))
    
        # now read harmonic bonds
        bonds_harmonic = []
        zmat['centroid_on_atom'] = ()
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip() + ' '*80
                if line[:4] == '    ': break
                try:
                    b = int(line[0:4])
                    e = int(line[4:8])
                    v1 = float(line[8:18])
                    v2 = float(line[18:28])
                    v3 = float(line[28:38])
                    v4 = float(line[38:48])
                except ValueError:
                    print(f'Fatal: wrong line: harmonic bonds\n  -> {line.rstrip()}')
                    return {}
                # check
                if b > num_atoms or e > num_atoms or b <= 0 or e <= 0:
                    print(f'Fatal: harmonic bonds: wrong value:\n  -> {line.rstrip()}')
                    return {}
                if e == 9999:
                    zmat['centroid_on_atom'] = (b-1,v2,v3,v4)
                else:
                    bonds_harmonic.append((b-1,e-1,v1,v2,v3,v4))
        zmat['bonds_harmonic'] = tuple(set(bonds_harmonic))

        # now read variable angles
        angles_variable = []
        angles_equal = []
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip() + ' '*80
                if line[:4] == '    ': break
                if line[4] not in [' ','-','=']:
                    print(f'Fatal: wrong line: variable angles\n  -> {line.rstrip()}')
                    return {}
                try:
                    if line[5:9] == '    ':
                        b = 1
                        e = int(line[0:4])
                    else:
                        b = int(line[0:4])
                        e = int(line[5:9])
                except ValueError:
                    print(f'Fatal: wrong line: variable angles\n  -> {line.rstrip()}')
                    return {}
                if b > num_atoms or e > num_atoms or b <= 0 or e <= 0:
                    print(f'Fatal: harmonic bonds: wrong value:\n  -> {line.rstrip()}')
                    return {}
                if line[4] == '=':
                    angles_equal.append((b-1,e-1))
                else:
                    if b > e: b, e = e, b
                    angles_variable.extend(list(range(b-1,e)))
        zmat['angles_variable'] = tuple(set(angles_variable))
        zmat['angles_equal'] = tuple(set(angles_equal))

        # now read additional angles
        angles_additional = []
        zmat['angles_auto'] = False
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip() + ' '*80
                if line[:4] == '    ': break
                if line[:4].lower() == 'auto':
                    zmat['angles_auto'] = True
                    continue
                try:
                    b = int(line[0:4])
                    e = int(line[4:8])
                    u = int(line[8:12])
                except ValueError:
                    print(f'Fatal: wrong line: additional angles\n  -> {line.rstrip()}')
                    return {}
                # check
                if b > num_atoms or e > num_atoms or b <= 0 or e <= 0 or u > num_atoms or u <= 0:
                    print(f'Fatal: additional angles: wrong value:\n  -> {line.rstrip()}')
                    return {}
                angles_additional.append((b-1,e-1,u-1))
        zmat['angles_additional'] = tuple(set(angles_additional))
    
        # now read dihedrals
        dihedrals_variable = []
        dihedrals_equal = []
        dihedrals_flip = []
        dihedrals_geometry = []
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip() + ' '*80
                if line[:4] == '    ': break
                if line[:4].lower() == 'flip':
                    try:
                        v = float(line[4:8])
                        e = int(line[8:12])
                    except ValueError:
                        print(f'Fatal: wrong line: variable dihedrals\n  -> {line.rstrip()}')
                        return {}
                    # check
                    if e > num_atoms or e <= 0:
                        print(f'Fatal: variable dihedrals: wrong value:\n  -> {line.rstrip()}')
                        return {}
                    dihedrals_flip.append((e-1,v))
                    continue
                if line[4] in [' ','-','=']:
                    try:
                        if line[5:9] == '    ':
                            b = 1
                            e = int(line[0:4])
                        else:
                            b = int(line[0:4])
                            e = int(line[5:9])
                    except ValueError:
                        print(f'Fatal: wrong line: variable dihedrals\n  -> {line.rstrip()}')
                        return {}
                    # check
                    if b > num_atoms or e > num_atoms or b <= 0 or e <= 0:
                        print(f'Fatal: variable dihedrals: wrong value:\n  -> {line.rstrip()}')
                        return {}
                    if line[4] == '=':
                        dihedrals_equal.append((b-1,e-1))
                    else:
                        if b > e: b, e = e, b
                        dihedrals_variable.extend(list(range(b-1,e)))
                    continue
                try:
                    b = int(line[0:4])
                    e = int(line[4:8])
                    u = int(line[8:12])
                    v = float(line[12:24])
                except ValueError:
                    print(f'Fatal: wrong line: geometry dihedrals\n  -> {line.rstrip()}')
                    return {}
                # check
                if b > num_atoms or b <= 0:
                    print(f'Fatal: geometry dihedrals: wrong value:\n  -> {line.rstrip()}')
                    return {}
                dihedrals_geometry.append((b,e,u,v))
        zmat['dihedrals_variable'] = tuple(set(dihedrals_variable))
        zmat['dihedrals_flip'] = tuple(set(dihedrals_flip))
        zmat['dihedrals_equal'] = tuple(set(dihedrals_equal))
        zmat['dihedrals_geometry'] = tuple(set(dihedrals_geometry))

        # now read additional dihedrals
        dihedrals_additional = []
        zmat['dihedrals_auto'] = False
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip() + ' '*80
                if line[:4] == '    ': break
                if line[:4].lower() == 'auto':
                    zmat['dihedrals_auto'] = True
                    continue
                try:
                    b = int(line[0:4])
                    e = int(line[4:8])
                    u = int(line[8:12])
                    f = int(line[12:16])
                    v1 = int(line[16:20])
                    v2 = int(line[20:24])
                except ValueError:
                    print(f'Fatal: wrong line: additional dihedrals\n  -> {line.rstrip()}')
                    return {}
                # check
                if b > num_atoms or e > num_atoms or b <= 0 or e <= 0 or u > num_atoms \
                    or f > num_atoms or u <= 0 or f <= 0:
                    print(f'Fatal: additional dihedrals: wrong value:\n  -> {line.rstrip()}')
                    return {}
                dihedrals_additional.append((b-1,e-1,u-1,f-1,v1,v2))
        zmat['dihedrals_additional'] = tuple(set(dihedrals_additional))
    
        # now read custom bond_pars
        bond_pars = []
        bo = False
        if line and line[:4] == '    ':      # double check
            for line in f:
                line = line.rstrip()
                if len(line) >= 4 and line[:4].lower() == ' fin':   # mark: ` Final`
                    bo = True
                    break
        if bo:
            for line in f:
                line = line.rstrip()
                if not line: continue
                line = line + ' '*80
                try:
                    t = int(line[:4])
                    s = int(line[5:7])
                    a = line[8:10].strip()
                    chg = float(line[11:21])
                    sig = float(line[21:31])
                    eps = float(line[31:41])
                except ValueError:
                    print(f'Warning: zmatrix bond parameters: {line.rstrip()}')
                else:
                    bond_pars.append((t,s,a,chg,sig,eps))
        zmat['bond_pars'] = tuple(bond_pars)

    #TODO
    assert len(bonds_geometry) > 0
    assert len(bonds_equal) == 0
    assert len(bonds_additional) == 0
    assert len(bonds_harmonic) == 0
    assert len(angles_equal) == 0
    assert zmat['angles_auto'] == False
    assert len(angles_additional) == 0
    assert len(dihedrals_additional) == 0
    assert len(dihedrals_flip) == 0
    assert len(dihedrals_geometry) == 0

    return zmat


def norcross(a,b):
    """calculate normalization for cross product a and b, 3x1"""
    x = a[1]*b[2] - b[1]*a[2]
    y = b[0]*a[2] - a[0]*b[2]
    z = a[0]*b[1] - b[0]*a[1]
    p = x*x + y*y + z*z
    s = 1.0/math.sqrt(p) if p > 0.0 else 0.0
    return [x*s, y*s, z*s]


def hessian(a,b,c):
    """calculate internal Hessian matrix, 3x3"""
    vca = [a[i]-c[i] for i in range(3)]
    vcb = [b[i]-c[i] for i in range(3)]
    y = norcross(vca,vcb)
    x = norcross(vcb,y)
    z = norcross(x,y)
    return [x,y,z]


def genatom(a,b,c,r,angle,dihedral):
    """a is the reference, b and c are anchors, angle/dihedral should be radian"""
    st = math.sin(angle)
    x = r * st * math.cos(dihedral)
    y = r * st * math.sin(dihedral)
    z = -r * math.cos(angle)
    h = hessian(a,b,c)
    xx = x*h[0][0] + y*h[1][0] + z*h[2][0] + c[0]
    yy = x*h[0][1] + y*h[1][1] + z*h[2][1] + c[1]
    zz = x*h[0][2] + y*h[1][2] + z*h[2][2] + c[2]
    return [xx,yy,zz]


def zmat2cor(zmat,radian=None):
    """convert zmatrix to coordinates
    
    Rule:
        firat atom will be (0.0, 0.0, 0.0)
        second atom will at positive x-axis
        third atom will in positive xy plane
    
    Return:
        cor: List[(x,y,z),   ]
    """
    if len(zmat) == 0: return []
    if len(zmat) == 1: return [0.0, 0.0, 0.0]

    if radian is not True:
        # avoid change the initial values
        nzmat = [[j for j in i] for i in zmat]
        for v in nzmat[1:]:
            v[3] = v[3] * math.pi / 180.0
            v[5] = v[5] * math.pi / 180.0
        zmat = nzmat

    cor = [[0.0,0.0,0.0], [zmat[1][1],0.0,0.0], ]
    if len(zmat) == 2: return cor

    y = zmat[2][1] * math.sin(zmat[2][3])
    x = zmat[2][1] * math.cos(zmat[2][3])
    if zmat[2][0] != 1:
        x = cor[1][0] - x
    cor.append([x,y,0.0])

    for v in zmat[3:]:
        xyz = genatom(cor[v[4]-1],cor[v[2]-1],cor[v[0]-1],v[1],v[3],v[5])
        cor.append(xyz)
    return cor


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


def find_all_rings_dfs(nbor):
    visited = [False for i in range(len(nbor))]
    dfs = []
    rings = []
    while False in visited:
        t = visited.index(False)
        visited[t] = True
        for v in nbor[t]:
            dfs.append([t,v])
        while len(dfs):
            g = dfs.pop()
            t = g[-1]
            if visited[t]:
                for v in nbor[t]:
                    if v == g[0]:
                        if len(g) > 2:
                            rings.append(g)
            else:
                visited[t] = True
                for v in nbor[t]:
                    if v == g[0]:
                        if len(g) > 2:
                            rings.append(g)
                    else:
                        for i in range(len(dfs)):
                            if dfs[i][-1] == t:
                                dfs[i].append(v)
                        for i in range(len(g)):
                            dfs.append([*g[i:],v])
    return rings


class ClosedForm:
    """Calculate RMSD for obmols by using closedForm algorithm

    Args:
        vl : 2D n*3f : List[[float, float, float]] : Left input matrix
        vr : 2D n*3f : List[[float, float, float]] : Right input matrix
        Left matrix is used for reference
            -> new_vr_fit_on_vl = BestFit * vr

    Methods:
        centroid
        calc_N
        calc_left_rotM
        calc_bestfit

    Reference:
        Closed-form solution of absolute orientation using unit quaternions
        Berthold K. P. Horn
        J. Opt. Soc. Am. A 4, 629-642 (1987)
    """
    def __init__(self,vl=None,vr=None,*args,**kws):
        self.vl = vl
        self.vr = vr

    def calc_bestfit(self,vl=None,vr=None,centroid=None,centroid_vl=True,centroid_vr=True):
        """calc best fit structure
        both vl(reference) and vr(target) can be reused
        centroid: whether centroid results: None(Yes), False(on vl), True(on vr)
        """
        if vl is None: vl = self.vl
        if vr is None: vr = self.vr
        if centroid_vl:
            cvl,tl = self.centroid(vl)
        else:
            cvl = vl
            tl = [0.0, 0.0, 0.0]
        if centroid_vr:
            cvr,tr = self.centroid(vr)
        else:
            cvr = vr
            tr = [0.0, 0.0, 0.0]
        M = self.calc_left_rotM(cvl,cvr)
        fit = []
        for v in cvr:
            rx = v[0]*M[0][0] + v[1]*M[1][0] + v[2]*M[2][0]
            ry = v[0]*M[0][1] + v[1]*M[1][1] + v[2]*M[2][1]
            rz = v[0]*M[0][2] + v[1]*M[1][2] + v[2]*M[2][2]
            fit.append([rx,ry,rz])
        if centroid is True:
            tk = tr
        elif centroid is False:
            tk = tl
        else:
            tk = [0.0, 0.0, 0.0]
        for i in range(len(fit)):
            fit[i][0] += tk[0]
            fit[i][1] += tk[1]
            fit[i][2] += tk[2]
        return fit

    def centroid(self,v):
        """calc centroid vector and translation vector"""
        x = [i[0] for i in v]
        y = [i[1] for i in v]
        z = [i[2] for i in v]
        t = len(v)
        ax = sum(x) / t
        ay = sum(y) / t
        az = sum(z) / t
        # care, do not use in-place operation
        cv = [[0.0 for i in range(3)] for j in range(len(v))]
        for i in range(t):
            cv[i][0] = v[i][0] - ax
            cv[i][1] = v[i][1] - ay
            cv[i][2] = v[i][2] - az
        return cv,(ax,ay,az)

    def calc_N(self,vl,vr):
        """calc 4x4 real symmetric N matrix"""
        XxYx = 0.0
        XxYy = 0.0
        XxYz = 0.0
        XyYx = 0.0
        XyYy = 0.0
        XyYz = 0.0
        XzYx = 0.0
        XzYy = 0.0
        XzYz = 0.0
        # careful of the sequence: X-r, Y-l
        # for rotation l = Rr
        for i,p in enumerate(vl):
            XxYx += p[0] * vr[i][0]
            XxYy += p[0] * vr[i][1]
            XxYz += p[0] * vr[i][2]

            XyYx += p[1] * vr[i][0]
            XyYy += p[1] * vr[i][1]
            XyYz += p[1] * vr[i][2]

            XzYx += p[2] * vr[i][0]
            XzYy += p[2] * vr[i][1]
            XzYz += p[2] * vr[i][2]

        N = [[0.0, 0.0, 0.0, 0.0] for i in range(4)]

        N[0][0] = XxYx + XyYy + XzYz
        N[0][1] = XyYz - XzYy
        N[0][2] = XzYx - XxYz
        N[0][3] = XxYy - XyYx

        N[1][0] = N[0][1]
        N[1][1] = XxYx - XyYy - XzYz
        N[1][2] = XxYy + XyYx
        N[1][3] = XzYx + XxYz

        N[2][0] = N[0][2]
        N[2][1] = N[1][2]
        N[2][2] = -XxYx + XyYy - XzYz
        N[2][3] = XyYz + XzYy

        N[3][0] = N[0][3]
        N[3][1] = N[1][3]
        N[3][2] = N[2][3]
        N[3][3] = -XxYx - XyYy + XzYz

        return N

    def calc_left_rotM(self,vl,vr):
        """calc rotation matrix for vr*M,
        quaternion is got from the vector which is corresponding to
        largest positive eigenvalue

        M  : 2D 3*3f : rotation matrix for vr, vr*M, in element-wise operation
        """
        N = self.calc_N(vl,vr)
        values,vectors = np.linalg.eig(N)
        ndx = np.where(values == max(values))
        ndx = ndx[0][0]
        # For numpy, eigenvectors are correspondingly put in column
        # note, this vector has already been normalized
        V = vectors[:,ndx]
        M = [[0.0, 0.0, 0.0] for i in range(3)]

        M[0][0] = 1 - 2 * (V[2]*V[2] + V[3]*V[3])
        M[0][1] = 2 * (V[1]*V[2] - V[3]*V[0])
        M[0][2] = 2 * (V[1]*V[3] + V[2]*V[0])

        M[1][0] = 2 * (V[1]*V[2] + V[3]*V[0])
        M[1][1] = 1 - 2 * (V[1]*V[1] + V[3]*V[3])
        M[1][2] = 2 * (V[2]*V[3] - V[1]*V[0])

        M[2][0] = 2 * (V[1]*V[3] - V[2]*V[0])
        M[2][1] = 2 * (V[2]*V[3] + V[1]*V[0])
        M[2][2] = 1 - 2 * (V[1]*V[1] + V[2]*V[2])

        return M


def align_onto_z_axis(cor):
    n = 0
    rmax = 0.0
    for i,l in enumerate(cor):
        r = l[0]*l[0] + l[1]*l[1] + l[2]*l[2]
        if r > rmax:
            rmax = r
            n = i
    v = cor[n]

    # first, rotate y axis to yz-plane
    t = pow(v[0]*v[0]+v[2]*v[2],0.5)
    ct = v[2] / t
    a = math.acos(ct)
    st = math.sin(a)
    if v[0] < 0.0: st = -st
    ry = [                  # counter clockwise
        [ct,  0.0,   st],
        [0.0, 1.0,  0.0],
        [-st, 0.0,   ct]
    ]

    # then, rotate x axis to xz-plane
    t = pow(rmax,0.5)
    st = v[1] / t
    a = math.asin(st)
    ct = math.cos(a)
    if v[1] < 0.0:
        st = -st
        ct = -ct
    rx = [                  # counter clockwise
        [1.0,   0.0, 0.0],
        [0.0,   ct,  st],
        [0.0,   -st,   ct]
    ]

    rmat = [
        [sum([ry[i][k]*rx[k][j] for k in range(3)]) for j in range(3)] for i in range(3)
    ]

    fin = []
    for c in cor:
        t = [sum([c[j]*rmat[j][i] for j in range(3)]) for i in range(3)]
        fin.append(t)

    return fin


def randu(x):
    imod = 1048573
    x = x * 1173 + 458753759
    x = x % imod
    y = x / imod
    return x,y


def randint(n1, n2, x=None):
    """number of `n2` replaced in `n1`, thus `n1>=n2`"""
    if n1 < n2: n1, n2 = n2, n1
    if not x: x = 786123
    nmr = [0 for i in range(n2)]
    for v in range(n2):
        while True:
            x,y = randu(x)
            t = int(n1*y) + 1
            bo = True
            for i in range(v):
                if nmr[i] == t:
                    bo = False
                    break
            if bo:
                nmr[v] = t
                break
    return x,[i-1 for i in sorted(nmr)]










