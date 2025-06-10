from PyBOSS.utils import randint, genatom
from .constants import SOLVENT_MODE, SOLVENT_PARS
from .inter import wxpot

import os
import math


def get_solvents_data(svmod1,svmod2,bond_pars):
    g_c_esq = 332.06
    solvent1_pars = SOLVENT_PARS[svmod1.lower()]
    solvent2_pars = SOLVENT_PARS[svmod2.lower()]
    s_mod1 = solvent1_pars[1]
    s_mod2 = solvent2_pars[1]

    spars = {}
    if s_mod1 == 1 or s_mod2 == 1:
        spars = SOLVENT_MODE[1]
    elif s_mod1 == 2 or s_mod2 == 2:
        spars = SOLVENT_MODE[2]

    solventpars = {
        '1': {
            'name': solvent1_pars[0],
            'modenum': solvent1_pars[1],
            'atom_names': solvent1_pars[2],
            'atom_indexes': solvent1_pars[3],
            'natoms': len(solvent1_pars[2])
        },
        '2': {
            'name': solvent2_pars[0],
            'modenum': solvent2_pars[1],
            'atom_names': solvent2_pars[2],
            'atom_indexes': solvent2_pars[3],
            'natoms': len(solvent2_pars[2])
        },
    }
    solventpars.update(spars)
    solventpars['natmx'] = max(len(solvent1_pars[2]),len(solvent2_pars[2]))

    sqrtesq = pow(g_c_esq,0.5)
    for s in ['1','2']:
        p = solventpars[s]
        bp = []
        for i in p['atom_indexes']:
            for t in bond_pars:
                if t[0] == i:
                    bp.append(t)
                    break
        p['bond_pars'] = bp
        p['QW'] = [l[3]*sqrtesq for l in bp]
        p['AW'] = [pow(4.0*l[5]*(l[4]**12),0.5) for l in bp]
        p['BW'] = [pow(4.0*l[5]*(l[4]**6), 0.5) for l in bp]

    inters = {}
    n = len(solventpars['1']['QW'])
    for i in range(n):
        q = solventpars['1']['QW'][i]
        a = solventpars['1']['AW'][i]
        b = solventpars['1']['BW'][i]
        for j in range(n):
            qq = solventpars['1']['QW'][j] * q
            aa = solventpars['1']['AW'][j] * a
            bb = solventpars['1']['BW'][j] * b
            key = f'{i}-{j}'
            inters[key] = [qq,aa,bb]
    solventpars['1-1'] = inters

    return solventpars


def read_waterboxes_tip4p(filename):
    """
        format:
            header: (I4,3F12.8)
            data:   (6F13.8)
    """
    if not os.path.isfile(filename):
        print(f'Warning: not a file: {filename}')
        return {}
    waterboxes = {}
    with open(filename,'rt') as f:
        while True:
            line = f.readline()
            if not line: break
            try:
                nmol = int(line[:4])
                bx = float(line[4:16])
                by = float(line[16:28])
                bz = float(line[28:40])
            except (ValueError,TypeError):      # TypeError when list empty
                print(f'Fatal: TIP4P waterbox: wrong input: {line}')
                return waterboxes
            data = []
            for l in range(nmol*4*3//6):
                line = f.readline()
                try:
                    a = float(line[:13])
                    b = float(line[13:26])
                    c = float(line[26:39])
                    d = float(line[39:52])
                    e = float(line[52:65])
                    t = float(line[65:78])
                except (ValueError,TypeError):
                    print(f'Fatal: TIP4P waterbox: wrong input: {line}')
                    return {}
                else:
                    data.extend([a,b,c,d,e,t])
            box = []
            s = 0
            for i in range(4):
                rl = []
                s += 3 * nmol
                for j in range(3):
                    u = s - (3-j) * nmol
                    rl.append([data[u+k] for k in range(nmol)])
                box.append(rl)
            mols = []
            for i in range(nmol):
                mols.append([
                    [box[k][0][i], box[k][1][i], box[k][2][i]] for k in range(4)
                ])
            waterboxes[nmol] = {
                'edge' : [bx,by,bz],
                'xyz' : mols
            }
    return waterboxes


def wxpot_init(
    nm, key, icut=None, waterbox=None, solventsdata=None, scut=None,
    solutezmat=None, solutesdata=None, movetype=None,
    **kws
):
    """function only used for initializing solvent input"""
    assert movetype != 2

    edge = waterbox['edge']
    wxyz = waterbox['xyz'][nm]

    wik = 1000000000.0
    rmins = [1000000000.0 for i in range(solutezmat['number_of_entries'])]
    for i,d in enumerate(solutesdata[key]['xyz']):
        if solutesdata[key]['atomtypes'][i] == 0: continue
        x = abs(d[0] - wxyz[0][0])
        if x > edge[0]: x -= edge[0]*2.0        # minimum image conversion
        y = abs(d[1] - wxyz[0][1])
        if y > edge[1]: y -= edge[1]*2.0
        z = abs(d[2] - wxyz[0][2])
        if z > edge[2]: z -= edge[2]*2.0

        for j,offset in enumerate(solutesdata['atoms_number_offset']):
            if i >= offset[0] and i < offset[1]:
                break
        u = x*x + y*y + z*z
        if u < rmins[j]:
            rmins[j] = u
        wik = min(u,wik)

    su2 = scut*scut
    if wik > su2:
        return [0.0, 0.0]

    sl2 = (scut-0.5)**2
    sul = 1.0 / (su2-sl2)
    iyes = [1 for i in range(len(solutezmat['data']))]

    if wik >= sl2:
        scale = sul * (su2-wik)
    else:
        scale = 1.0
    if icut != 0:
        for i in range(solutezmat['number_of_entries']):
            if rmins[i] > su2:
                for j in range(*solutesdata['atoms_number_offset'][i]):
                    iyes[j] = 0
    aw = solventsdata['1']['AW']
    bw = solventsdata['1']['BW']
    qw = solventsdata['1']['QW']
    anew = solutesdata[key]['xyz']
    a = solutesdata[key]['A']
    b = solutesdata[key]['B']
    q = solutesdata[key]['Q']
    elj = ecoul = 0.0
    for i,p in enumerate(solutezmat['data']):
        if p[1] <= 0 or iyes[i] == 0: continue
        x = abs(anew[i][0] - wxyz[0][0])         # LJ only for oxygen atom
        y = abs(anew[i][1] - wxyz[0][1])
        z = abs(anew[i][2] - wxyz[0][2])
        xim = 0.0 if x <= edge[0] else edge[0]*2.0
        yim = 0.0 if y <= edge[1] else edge[1]*2.0
        zim = 0.0 if z <= edge[2] else edge[2]*2.0
        rr = (x-xim)**2 + (y-yim)**2 + (z-zim)**2
        r1 = pow(rr,0.5)
        r6 = 1.0 / (rr*rr*rr)
        elj += r6 * (a[i]*aw[0]*r6-b[i]*bw[0]) * scale
        ecoul += q[i] * qw[0] * scale / r1
        for j in range(1,len(wxyz)):
            x = abs(anew[i][0] - wxyz[j][0])
            y = abs(anew[i][1] - wxyz[j][1])
            z = abs(anew[i][2] - wxyz[j][2])
            rr = (x-xim)**2 + (y-yim)**2 + (z-zim)**2
            r1 = pow(rr,0.5)
            ecoul += q[i] * qw[j] * scale / r1
    return [elj,ecoul]


class InitSolvents:
    def __init__(
        self, svmod1=None, svmod2=None, bond_pars=None,
        waterboxfile=None, ibox=None, irn=None, nmol=None, nmol2=None,
        *args, **kws
    ):
        self.nice = True
        self.info = ''

        waterboxes = read_waterboxes_tip4p(waterboxfile)

        self.ibox = ibox if ibox else 216
        self.waterbox = waterboxes[self.ibox]
        self.kws['ac'] = self.waterbox['xyz']
        self.kws['edge'] = e = self.waterbox['edge']
        self.kws['edg2'] = e / 2.0

        self.irn = irn if irn else 786123
        self.nmol = nmol
        self.nmol2 = nmol2
        self.solventsdata = get_solvents_data(svmod1, svmod2, bond_pars)

        self.kws = kws
        self.kws['aw'] = []
        self.kws['bw'] = []
        self.kws['qw'] = []
        for k in ['1','2']:
            self.kws['aw'].append(self.solventsdata[k]['aw'])
            self.kws['bw'].append(self.solventsdata[k]['bw'])
            self.kws['qw'].append(self.solventsdata[k]['qw'])

    def run(self,key=None):
        self.buildup()
        self.boss_type_filter()
        self.replace_solvent_atoms()
        if self.solventsdata['1']['modenum'] <= 2:
            self.fix_solvent_waterbox()
        self.solventsdata.update(self.waterbox)

    def buildup(self,key=None):
        if key is None: key = 'reference'
        if key != 'reference': assert False
        if key:
            self.energies = []

            self.kws['movtyp'] = 1
            self.kws['ncutat'] = None
            self.kws['nmov'] = None

            for nm in range(self.ibox):
                self.kws['nm'] = nm
                energy = wxpot(**self.kws)
                self.energies.append(energy)
        else:
            assert False

    def boss_type_filter(self):
        # to maximumly keep consistence
        esw = [sum(self.energies[i]) for i in range(self.ibox)]
        sortedndx = list(range(self.ibox))
        for i in range(self.ibox):
            for j in range(i):
                if esw[j] >= esw[i]:
                    e = esw.pop(i)
                    esw.insert(j,e)
                    n = sortedndx.pop(i)
                    sortedndx.insert(j,n)
        # to keep consistent
        newndx = sorted([sortedndx[i] for i in range(self.nmol)])
        # update
        self.waterbox['xyz'] = [self.waterbox['xyz'][i] for i in newndx]
        return sortedndx

    def replace_solvent_atoms(self,irn=None):
        if not irn: irn = self.irn
        natoms2 = len(self.solventsdata['2']['atom_names'])
        irn, nmr = randint(self.nmol,self.nmol2,irn)
        anew = [[0.0, 0.0, 0.0] for i in range(50)]
        anew[1][0] = 0.945
        anew[2][0] = -0.45374566
        anew[2][1] = -1.35610283
        for i in nmr:
            xyz = self.waterbox['xyz'][i][0]
            for j in range(1,natoms2):
                for k in range(3):
                    self.waterbox['xyz'][i][j][k] = anew[j][k] + xyz[k]
        self.irn = irn
        self.nmr = nmr
        self.waterbox['nmr'] = nmr

    def fix_solvent_waterbox(self):
        rcd = 0.9572
        thbcd = 104.52 / 180.0 * math.pi
        rcm = 0.15
        th2 = thbcd / 2.0
        wbox = self.waterbox['xyz']
        for i in range(self.nmol):
            if i in self.nmr: continue
            a = wbox[i][1]
            b = wbox[i][0]
            c = wbox[i][2]
            d = genatom(a,b,c,rcd,thbcd,0.0)
            a = genatom(d,c,b,rcd,thbcd,0.0)
            c = genatom(d,a,b,rcd,thbcd,0.0)
            for j in range(3):
                wbox[i][1][j] = a[j]
                wbox[i][2][j] = c[j]
                if len(wbox[i]) == 4:
                    e = genatom(d,a,b,rcm,th2,0.0)
                    wbox[i][3][j] = e[j]
        self.waterbox['xyz'] = wbox






