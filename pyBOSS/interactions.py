

class Interactions:
    def __init__(
        self, solventsdata=None, solutezmat=None, solutesdata=None, movetype=None,
        islab=None, rcut=None, irfon=None, dielrf_factor=None, ncutas=None,
        icutas=None, icutat=None, nmol=None, nmol2=None, scut=None, icut=None,
        ncent1=None, ncent2=None, rdfs_solute_solvent=None, isolec=None,
        wkc=None, essmin=None, essinc=None, nocut=None, noss=None,
        *args, **kws
    ):
        self.solutezmat = solutezmat
        self.solutesdata = solutesdata

        self.solventsdata = solventsdata
        self.edge = solventsdata['edge']
        self.aoo = solventsdata['AOO']
        self.coo = solventsdata['COO']
        self.qq = solventsdata['QQ']
        self.kq = solventsdata['KQ']
        self.nmr = solventsdata['nmr']
        self.edge2 = [i*2.0 for i in self.edge]
        self.vnew = self.edge2[0] * self.edge2[1] * self.edge2[2]

        self.islab = islab
        self.rcut = rcut
        self.rcutsq = self.rcut * self.rcut
        self.rl2 = (self.rcut-0.5) * (self.rcut-0.5)
        self.rul = 1.0 / (self.rcutsq - self.rl2)
        self.scut = scut
        self.sl2 = (self.scut-0.5) * (self.scut-0.5)
        self.scutsq = self.scut * self.scut
        self.sul = 1.0 / (self.scutsq-self.sl2)
        self.rc3 = 1.0 / (self.rcutsq*self.rcut)
        self.rc9 = self.rc3 * self.rc3 * self.rc3

        self.irfon = irfon
        self.rfact = dielrf_factor
        self.ncutas = ncutas
        self.icutas = icutas
        self.icutat = icutat
        self.icut = icut
        self.nmol = nmol if nmol else len(self.solventsdata['xyz'])
        self.nmol2 = nmol2 if nmol2 else len(self.nmr)
        self.myncent1 = ncent1 - 1        # input index starts from 1
        self.myncent2 = ncent2 - 1
        self.rdfs_solute_solvent = rdfs_solute_solvent
        self.isolec = isolec
        self.movetype = movetype
        self.mymaxwik = 100000000.0
        self.wkc = wkc
        self.essmin = essmin
        self.essinc = essinc
        self.nocut = nocut
        self.noss = noss
        self.slvatoms1 = self.solventsdata['1']['natoms']
        self.slvatoms2 = self.solventsdata['2']['natoms']
        self.slvnatoms = max(self.slvatoms1, self.slvatoms2)
        self.slvmode1 = self.solventsdata['1']['modenum']
        self.slvmode2 = self.solventsdata['2']['modenum']
        self.twopi = 3.141592653589793 * 2.0


    def sspot(self, nm, nmov, anew, keyi, keyj):
        wxyz = self.solventsdata['xyz'][nm]

        nt1 = 1 if nmov in self.nmr else 0
        nt2 = 1 if nm in self.nmr else 0

        ai = self.solventsdata[keyi]['AW']
        bi = self.solventsdata[keyi]['BW']
        qi = self.solventsdata[keyi]['QW']
        bpi = self.solventsdata[keyi]['bond_pars']
        aj = self.solventsdata[keyj]['AW']
        bj = self.solventsdata[keyj]['BW']
        qj = self.solventsdata[keyj]['QW']
        bpj = self.solventsdata[keyj]['bond_pars']

        wik = self.mymaxwik
        for i in range(self.ncutas[nt1]):
            ni = self.icutas[nt1][i]
            for j in range(self.ncutas[nt2]):
                nj = self.icutas[nt2][j]
                zz = abs(anew[ni][2]-wxyz[nj][2])
                izz = self.edge2[2] if zz > self.edge[2] else 0.0
                mwik = (zz-izz) * (zz-izz)
                if mwik <= self.rcutsq:
                    xx = abs(anew[ni][0]-wxyz[nj][0])
                    ixx = self.edge2[0] if xx > self.edge[0] else 0.0
                    mwik += (xx-ixx) * (xx-ixx)
                    if mwik <= self.rcutsq:
                        yy = abs(anew[ni][1]-wxyz[nj][1])
                        iyy = self.edge2[1] if yy > self.edge[1] else 0.0
                        mwik += (yy-iyy) * (yy-iyy)
                        wik = min(mwik, wik)
        if wik > self.rcutsq:
            return [wik, 0.0, 0.0]

        sspot = sslj = 0.0
        if self.nmol2 > 0:
            for i in range(len(ai)):
                if bpi[i][0] <= 0: continue
                aii = ai[i]
                bii = bi[i]
                qii = qi[i]
                for j in range(len(aj)):
                    if bpj[j][0] <= 0: continue
                    xx = abs(anew[i][0]-wxyz[j][0])
                    ixx = self.edge2[0] if xx > self.edge[0] else 0.0
                    yy = abs(anew[i][1]-wxyz[j][1])
                    iyy = self.edge2[1] if yy > self.edge[1] else 0.0
                    zz = abs(anew[i][2]-wxyz[j][2])
                    izz = self.edge2[2] if zz > self.edge[2] else 0.0
                    rr = (xx-ixx)*(xx-ixx) + (yy-iyy)*(yy-iyy) + (zz-izz)*(zz-izz)
                    if aii != 0.0:
                        r6 = 1.0 / (rr*rr*rr)
                        aaij = aii * aj[j]
                        bbij = bii * bj[j]
                        xlj = (aaij*r6-bbij) * r6
                        sspot += xlj
                        sslj += xlj
                    qqij = qii * qj[j]
                    if qqij != 0.0:
                        r1 = 1.0 / pow(rr,0.5)
                        if self.irfon == 0:
                            sspot += qqij * r1
                        else:
                            sspot += qqij * (r1 + rr*self.rfact)
        else:
            assert False
        if wik > self.rl2:
            factor = self.rul * (self.rcutsq-wik)
            sspot = sspot * factor
            sslj = sslj * factor
        
        return [wik, sspot, sslj]

    def wwpot(self, nm, anew):
        wxyz = self.solventsdata['xyz'][nm]

        zz = abs(anew[0][2] - wxyz[0][2])
        if self.islab == 1:
            izz = self.edge2[2] if zz > self.edge[2] else 0.0
            wik = (zz-izz) * (zz-izz)
        else:
            wik = zz * zz

        bo = True
        if wik > self.rcutsq:
            bo = False
        else:
            xx = abs(anew[0][0] - wxyz[0][0])
            ixx = self.edge2[2] if xx > self.edge[2] else 0.0
            wik += (xx-ixx) * (xx-ixx)
            if wik > self.rcutsq:
                bo = False
            else:
                yy = abs(anew[0][1] - wxyz[0][1])
                iyy = self.edge2[2] if yy > self.edge[2] else 0.0
                wik += (yy-iyy) * (yy-iyy)
                if wik > self.rcutsq:
                    bo = False
        if not bo:
            return [wik, 0.0, 0.0]

        if bo:
            r6 = 1.0 / (wik*wik*wik)
            wwpot = r6 * (self.aoo*r6-self.coo)
            sslj = wwpot
            n = len(anew)
            knt = 0
            for i in range(self.kq,n):
                for j in range(self.kq,n):
                    xx = abs(anew[i][0]-wxyz[j][0])
                    yy = abs(anew[i][1]-wxyz[j][1])
                    zz = abs(anew[i][2]-wxyz[j][2])
                    rr = (xx-ixx)**2 + (yy-iyy)**2 + (zz-izz)**2
                    assert rr > 0.0
                    r1 = 1.0 / pow(rr,0.5)
                    if self.irfon == 0:
                        wwpot += self.qq[knt] * r1
                    else:
                        wwpot += self.qq[knt] * (r1 + rr*self.rfact)
                    knt += 1
            if wik >= self.rl2:
                wwpot = wwpot * self.rul * (self.rcutsq-wik)
        else:
            wwpot = 0.0
            sslj = 0.0
        return [wik, wwpot,sslj]

    def run_solvent_solvent_interactions(self):
        self.energies = {}
        self.wiks = []
        for i in range(self.nmol-1):
            if i in self.solventsdata['nmr']:
                keyi = '2'
                svmodi = self.solventsdata['2']['modenum']
            else:
                keyi = '1'
                svmodi = self.solventsdata['1']['modenum']
            anew = self.solventsdata['xyz'][i]
            for j in range(i+1,self.nmol):
                if len(self.solventsdata['2']['atom_names']) == 0:
                    assert False
                else:
                    if j in self.solventsdata['nmr']:
                        keyj = '2'
                        svmodj = self.solventsdata['2']['modenum']
                    else:
                        keyj = '1'
                        svmodj = self.solventsdata['1']['modenum']
                    if svmodi <= 2 and svmodj <= 2:
                        e = self.wwpot(j, anew)
                    else:
                        e = self.sspot(j, i, anew, keyi, keyj)
                    self.energies[f'{i}-{j}'] = [*e[1:], keyi, keyj]
                    self.wiks.append(e[0])


    def sxpot(self, nm, nmov=None, movetype=None, icut=None):
        if movetype is None: movetype = self.movetype
        if icut is None: icut = self.icut

        centerasol = self.solutesdata['reference']['xyz']
        n1 = 1 if nm in self.nmr else 0

        if movetype == 2:
            assert False
            if nmov in self.nmr:
                skey = '2'
                n1 = 1
            else:
                skey = '1'
                n1 = 0
            natoms = self.solventsdata[skey]['natoms']
        else:
            anm = self.solventsdata['xyz'][nm]
            if nm in self.nmr:
                skey = '2'
            else:
                skey = '1'
            aw = self.solventsdata[skey]['AW']
            bw = self.solventsdata[skey]['BW']
            qw = self.solventsdata[skey]['QW']

        if self.icut < 2:
            assert False
            wik = self.mymaxwik
            for j in range(self.ncutas[n1]):
                nc = self.icutas[n1][j]
                twik = 0.0
                for i in range(3):
                    if movetype == 1 or movetype == 2:
                        cnew = 0.5 * (centerasol[self.myncent1][i]+centerasol[self.myncent2][i])
                        cnm = anm[nc][i]
                    else:
                        cnew = centerasol[nc][i]
                        cnm  = anm[nc][i]
                    x = abs(cnew-cnm)
                    if self.islab == 1 and i == 3: break
                    if x > self.edge[i]: x = x - self.edge2[i]
                    twik += x*x
                wik = min(wik, twik)
        else:
            rmins = [self.mymaxwik for i in range(self.solutezmat['number_of_entries'])]
            wik = self.mymaxwik
            for i,d in enumerate(self.solutesdata['initial_bond_pars']):
                if d[1] == 0: continue

                if self.icut >= 3:
                    if i+1 not in self.icutat: continue
                
                for m,offset in enumerate(self.solutesdata['atoms_number_offset']):
                    if i >= offset[0] and i < offset[1]:
                        break
                
                for nc in self.icutas[n1]:
                    x = 0.0
                    for j in range(3):
                        y = abs(centerasol[i][j]-anm[nc][j])
                        if self.islab == 1 and j == 3: break
                        x += y*y
                    rmins[m] = min(rmins[m], x)
                    wik = min(wik,x)

        elj = [0.0, 0.0, 0.0]
        ecoul = [0.0, 0.0, 0.0]
        ecsx = elsx = 0.0
        if wik > self.scutsq:
            return [wik, elj, ecoul, ecsx, elsx]

        assert False
        iyes = [1 for i in range(len(self.solutezmat['data']))]
        if wik >= self.sl2:
            scale = self.sul * (self.scutsq-wik)
        else:
            scale = 1.0

        if self.icut != 0:
            for i in range(self.solutezmat['number_of_entries']):
                if rmins[i] > self.scutsq:
                    for j in range(*self.solutesdata['atoms_number_offset'][i]):
                        iyes[j] = 0

        elj = [0.0, 0.0, 0.0]
        ecoul = [0.0, 0.0, 0.0]
        ecsx = elsx = 0.0
        rnew = []
        for n,key in enumerate(['reference','first','second']):
            a = self.solutesdata[key]['A']
            b = self.solutesdata[key]['B']
            q = self.solutesdata[key]['Q']
            asol = self.solutesdata[key]
            for i,p in enumerate(self.solutezmat['data']):
                if p[1] <= 0 or iyes[i] == 0: continue

                for m,offset in enumerate(self.solutesdata['atoms_number_offset']):
                    if i >= offset[0] and i < offset[1]:
                        break
                
                x = abs(asol[i][0][0]-anm[0][0])
                y = abs(asol[i][0][1]-anm[0][1])
                z = abs(asol[i][0][2]-anm[0][2])
                xim = self.edge2[0] if x > self.edge[0] else 0.0
                yim = self.edge2[1] if y > self.edge[1] else 0.0
                zim = self.edge2[2] if z > self.edge[2] else 0.0
                for j in range(natoms):
                    x = abs(asol[i][j][0]-anm[j][0])
                    y = abs(asol[i][j][1]-anm[j][1])
                    z = abs(asol[i][j][2]-anm[j][2])
                    rr = (x-xim)**2 + (y-yim)**2 + (z-zim)**2
                    r1 = pow(rr,0.5)
                    r6 = 1.0 / (rr*rr*rr)
                    el = scale * r6 * (a[i]*aw[j]) * r6 - b[i]*bw[j]
                    elj[n] += el
                    ec = scale * q[i] * qw[j]
                    if self.irfon == 0:
                        ec = ec / r1
                    else:
                        ec = ec * (1.0/r1 + rr*self.rfact)
                    ecoul[n] += ec
                    if n == 1:
                        if m == self.isolec:
                            ecsx += ec
                            elsx += el
                        if [i+1,j+1] in self.rdfs_solute_solvent:
                            rnew.append(r1)
        return [wik, elj, ecoul, ecsx, elsx]

    def wxpot(self, nm, key=None, movetype=None, icut=None):
        if movetype is None: movetype = self.movetype
        if icut is None: icut = self.icut
        if key is None: key = 'reference'

        if movetype == 1:
            centerasol = self.solventsdata['xyz'][0]
        elif movetype == 2:
            assert False
            centerasol = self.solutesdata['reference']['xyz']

        anm = self.solventsdata['xyz'][nm]
        wik = self.mymaxwik
        rmins = [self.mymaxwik for i in range(self.solutezmat['number_of_entries'])]

        if self.icut < 2:
            assert False
            wik = 0.0
            for i in range(3):
                if movetype == 1 or movetype == 2:
                    cnew = 0.5 * (centerasol[self.myncent1][i]+centerasol[self.myncent2][i])
                else:
                    assert False
                x = abs(cnew-anm[0][i])
                if self.islab == 1 and i == 3: break
                if x > self.edge[i]: x = x - self.edge2[i]
                wik += x*x
        else:
            for i,d in enumerate(self.solutesdata['reference']['xyz']):
                if self.solutesdata['initial_bond_pars'][i][1] == 0: continue

                if self.icut >= 3:
                    if i+1 not in self.icutat: continue

                x = abs(d[0] - anm[0][0])
                if x > self.edge[0]: x -= self.edge2[0]     # minimum image conversion
                y = abs(d[1] - anm[0][1])
                if y > self.edge[1]: y -= self.edge2[1]
                z = abs(d[2] - anm[0][2])
                if z > self.edge[2]: z -= self.edge2[2]

                for j,offset in enumerate(self.solutesdata['atoms_number_offset']):
                    if i >= offset[0] and i < offset[1]:
                        break
                u = x*x + y*y + z*z
                if u < rmins[j]:
                    rmins[j] = u
                wik = min(u,wik)

        if wik > self.scutsq:
            return [wik, 0.0, 0.0, 0.0, 0.0]

        assert False
        iyes = [1 for i in range(len(self.solutezmat['data']))]
        if wik >= self.sl2:
            scale = self.sul * (self.scutsq-wik)
        else:
            scale = 1.0
        if self.icut != 0:
            for i in range(self.solutezmat['number_of_entries']):
                if rmins[i] > self.scutsq:
                    for j in range(*self.solutesdata['atoms_number_offset'][i]):
                        iyes[j] = 0
        aw = self.solventsdata['1']['AW']
        bw = self.solventsdata['1']['BW']
        qw = self.solventsdata['1']['QW']
        if movetype != 2:
            anew = self.solutesdata[key]['xyz']
            a = self.solutesdata[key]['A']
            b = self.solutesdata[key]['B']
            q = self.solutesdata[key]['Q']
            elj = ecoul = 0.0
            for i,p in enumerate(self.solutezmat['data']):
                if p[1] <= 0 or iyes[i] == 0: continue
                x = abs(anew[i][0] - anm[0][0])         # LJ only for oxygen atom
                y = abs(anew[i][1] - anm[0][1])
                z = abs(anew[i][2] - anm[0][2])
                xim = 0.0 if x <= self.edge[0] else self.edge2[0]
                yim = 0.0 if y <= self.edge[1] else self.edge2[1]
                zim = 0.0 if z <= self.edge[2] else self.edge2[2]
                rr = (x-xim)**2 + (y-yim)**2 + (z-zim)**2
                r1 = pow(rr,0.5)
                r6 = 1.0 / (rr*rr*rr)
                elj += r6 * (a[i]*aw[0]*r6-b[i]*bw[0]) * scale
                ecoul += q[i] * qw[0] * scale / r1
                for j in range(1,len(anm)):
                    x = abs(anew[i][0] - anm[j][0])
                    y = abs(anew[i][1] - anm[j][1])
                    z = abs(anew[i][2] - anm[j][2])
                    rr = (x-xim)**2 + (y-yim)**2 + (z-zim)**2
                    r1 = pow(rr,0.5)
                    ecoul += q[i] * qw[j] * scale / r1
        else:
            assert False

        return [wik, elj, ecoul, 0.0, 0.0]

    def run_solute_solvent_interactions(self,movetype=None,icut=None):
        self.energies = []
        self.wiks = []
        for i in range(self.nmol):
            if i in self.nmr:
                modsv = self.solventsdata['2']['modenum']
            else:
                modsv = self.solventsdata['1']['modenum']
            if modsv > 2:
                e = self.sxpot(i,movetype=movetype,icut=icut)
            else:
                e = self.wxpot(i,movetype=movetype,icut=icut)
            self.wiks.append(e[0])
            self.energies.append(e[1:])

    def normlize_wiks(self):
        self.wiks = [1.0/(i+self.wkc) for i in self.wiks]
        total = sum(self.wiks)
        self.wiks = [i/total for i in self.wiks]
        self.maxwik = max(self.wiks)
    
    def calc_overflow_energies(self):
        self.ioess = [0 for i in range(50)]
        for e in [sum(i) for i in self.energies]:
            j = int((e-self.essmin)/self.essinc)
            if j > 0 and j <= 50:
                self.ioess[j-1] += 1

    def ecut(self):
        if self.slvnatoms <= 0 or self.nocut == 1 or self.noss == 1:
            return 0.0

        if self.slvatoms2 == 0:
            if self.slvmode1 <= 2 or self.slvmode1 == 13:
                return 0.0
            xn = self.nmol
            xm = 0.0
        else:
            xn = self.nmol - self.nmol2
            xm = self.nmol2
            if self.slvmode1 <= 2: xn = 0.0
            if self.slvmode2 <= 2: xm = 0.0

        xnn = xn * xn
        xmm = xm * xm
        xnm = xn * xm
        ecut = 0.0
        for i in range(self.slvatoms1):
            ai = self.solventsdata['1']['AW'][i]
            bi = self.solventsdata['1']['BW'][i]
            for j in range(self.slvatoms1):
                aj = self.solventsdata['1']['AW'][j]
                bj = self.solventsdata['1']['BW'][j]
                ecut += self.twopi*xnn/(3.0*self.vnew) * (ai*aj*self.rc9/3.0 - bi*bj*self.rc3)

        for i in range(self.slvatoms2):
            ai = self.solventsdata['2']['AW'][i]
            bi = self.solventsdata['2']['BW'][i]
            for j in range(self.slvatoms2):
                aj = self.solventsdata['2']['AW'][j]
                bj = self.solventsdata['2']['BW'][j]
                ecut += self.twopi*xmm/(3.0*self.vnew) * (ai*aj*self.rc9/3.0 - bi*bj*self.rc3)

        for i in range(self.slvatoms1):
            ai = self.solventsdata['1']['AW'][i]
            bi = self.solventsdata['1']['BW'][i]
            for j in range(self.slvatoms2):
                aj = self.solventsdata['2']['AW'][j]
                bj = self.solventsdata['2']['BW'][j]
                ecut += self.twopi*xnm/(3.0*self.vnew) * (ai*aj*self.rc9/3.0 - bi*bj*self.rc3)

        return ecut









