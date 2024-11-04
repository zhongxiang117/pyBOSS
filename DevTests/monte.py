from PyBOSS.utils import randu
from PyBOSS.boss_calc import PyBOSSCalc
from PyBOSS.interactions import Interactions
from PyBOSS.constants import PAR_ENERGIES

import math


class PyMonte:
    def __init__(
        self,solventsdata=None,solutesdata=None,solutezmat=None,ncons=None,
        innerpars=None,

        irn=None,nmol=None,ncutas=None,icutas=None,islab=None,nmr=None,
        rcut=None,eprinc=None,edinc=None,essinc=None,ebsinc=None,beta=None,
        nschg=None,nvchg=None,nopref=None,pvcon=None,nofep=None,
        nosolu=None,nsolut=None,igbsa=None,lhtsol=None,betlht=None,
        *args,**kws
    ):
        assert nosolu is None or nosolu == 0
        self.nosolu = 0
        assert nsolut is None or nsolut == 1
        self.nsolut = 1
        assert igbsa is None or igbsa == 0
        self.igbsa = 0

        self.solventsdata=solventsdata
        self.solutesdata= solutesdata
        self.solutezmat = solutezmat
        self.innerpars = innerpars
        self.irn = irn
        self.nmol = nmol
        self.ncutas = ncutas
        self.icutas = icutas
        self.islab = islab
        self.nmr = nmr if nmr else self.solventsdata['nmr']
        self.rcut = rcut
        self.nschg = nschg
        self.nvchg = nvchg
        self.ncons = ncons
        self.beta = beta
        self.nopref = nopref
        self.nfl = 0
        self.nsolmv = 0
        self.eprinc = eprinc
        self.edinc = edinc
        self.essinc = essinc
        self.ebsinc = ebsinc
        self.ntdih = len(self.solutesdata['pert_delta_dihedrals'])
        self.ntbnd = len(self.solutesdata['pert_delta_bonds'])
        self.ntang = len(self.solutesdata['pert_delta_angles'])
        self.natmx = self.solventsdata['natmx']
        self.pvcon = pvcon
        self.nofep = nofep
        self.lhtsol = lhtsol
        self.betlht = betlht

        self.xmlbet = (self.nmol+self.nsolut) / self.beta
        self.epin = 1.0 / self.eprinc
        self.dince = 1.0 / self.edinc
        self.epins = 1.0 / self.essinc
        self.dinces = 1.0 / self.ebsinc
        self.mhist = self.ncons//(120*self.nschg) + 1
        self.ihist = 0
        self.ndihmx = 8 if self.ntdih > 8 else self.ntdih
        self.rncsq = (self.rcut+1.5)**2
        self.tk20 = 20.0 / self.beta
        self.tk25 = 25.0 / self.beta
        self.nschg5 = 5 * self.nschg
        self.nsolup = 20 * self.nschg
        self.nfrqsa = min(1+self.ncons//10,self.nsolup)
        self.ncon = 0
        self.movvol = 0
        self.movsol = 0
        self.movtyp = 0
        self.histdi = [[0.0 for i in range(8)] for j in range(120)]
        self.nflptr = 0
        self.nflprj = 0
        self.nsolmv = 0
        self.E = {k:v[0] for k,v in PAR_ENERGIES.items()}
        for k,v in self.E.items():
            if v is None:
                self.E[k] = []


    def setup(self):
        self.calc_hydrogen_da()

        self.myinter = Interactions(
            solventsdata=self.solventsdata, solutezmat=self.solutezmat,
            solutesdata=self.solutesdata,
            **self.innerpars
        )
        self.myinter.run_solute_solvent_interactions(movetype=1)
        self.myinter.normlize_wiks()

        ecut = self.myinter.ecut()

        radsa = self.solutesdata['radius_sa']
        asol = self.solutesdata['reference']['xyz']

        self.mypbc = PyBOSSCalc(radsa=radsa, asol=asol, rsolv=self.innerpars['rsolv'])

        sasa = self.mypbc.calc_amber_sasa(self.solutesdata['amber_sasa_atomtypes'])

        self.solutesdata['sasa'] = sasa






    def run(self):
        nsatm = len(self.solutesdata['data'])
        self.ncon += 1
        nrepet = 1

        etot = 0.0
        hinter = 0.0
        v = 0.0
        vsq = 0.0
        h = 0.0
        hinsq = 0.0
        vh = 0.0

        edelg = 0.0
        udel1 = 0.0
        edelg2 = 0.0
        udel2 = 0.0
        sos1 = sos2 = 0.0

        lhtsol = self.lhtsol
        self.isol = 0

        enew = 0.0
        eold = 0.0
        vnew = 0.0
        vold = 0.0

        xr = [[0 for i in range(15)] for j in range(100)]
        xrs = [[0 for i in range(15)] for j in range(100)]
        ior = [[0 for i in range(15)] for j in range(100)]
        iors = [[0 for i in range(15)] for j in range(100)]
        xedist = [0.0 for i in range(50)]
        xepar = [0.0 for i in range(50)]
        xess = [0.0 for i in range(50)]
        ioepar = [0 for i in range(50)]
        ioess = [0 for i in range(50)]
        xbesol = [0.0 for i in range(50)]
        dtorad = math.pi / 180.0

        if self.ncon % self.nschg5 == 0:
            self.calc_solutes_dipole()
            if self.nmol != 0:
                self.calc_hydrogen_da()

        if self.ncon % self.nfrqsa == 0:
            if self.nosolu != 0:
                # use PyBOSSCalc
                pass
            self.ssljco()
        
        bo70 = False
        if self.ncon % self.nvchg == 0 or self.nvchg == 999999:
            ##goto 40
            movvol = 1
            self.irn, y = randu(self.irn)
            nmov = int(self.nmol*y+1)
            #goto 70
            bo70 = True
        elif self.ncon % self.nschg == 0:
            ##goto 50
            movtyp = 1
            nmov = 0
            self.movsol = 1
            if self.ntdih != 0:
                if (self.ncon//self.nschg) % self.mhist == 0:
                    for i in range(self.ndihmx):
                        tmp = self.solutezmat['data'][i][8]
                        self.histdi[self.ihist][i] = 1 + int(48.0*tmp/math.pi/2.0)
            bo70 = True
        else:
            while True:
                self.irn, y = randu(self.irn)
                nmov = int(y+1.0)
                if self.nopref == 0:
                    self.irn, y = randu(self.irn)
                    wr = self.myinter.wiks[nmov] / self.myinter.maxwik
                    if wr > y:
                        break
                else:
                    break
            self.movsvn()
            #goto 80
        if bo70:
            self.movmol()
            if self.movvol == 1:
                pass
            if self.nflptr < 0:
                self.nflptr = -self.nflptr
        #80
        delh = enew - eold
        fac = 1.0
        if self.movtyp != 0:
            if self.movsol != 1:
                delh += self.pvcon*(vnew-vold) - self.xmlbet*math.log(vnew/vold)
            if delh <= 0.0:
                pass
                #goto 90
            else:
                if self.nopref == 0:
                    pass
                    #fac = wnew(nmov) / self.myinter.wiks[nmov]
        if delh > self.tk20:
            #goto 550
            pass
        elif delh < -self.tk20:
            #goto 90
            pass
        else:
            self.irn, y = randu(self.irn)
            if self.movsol == 1:
                if self.isol == lhtsol or lhtsol == 9999:
                    xbeta = self.betlht
            else:
                xbeta = self.beta
            emet = fac * math.exp(-xbeta*delh)
            if emet < y:
                #goto 550
                pass
        #90
        xrep = nrepet + 0.0
        x = xrep * eold
        etot += x
        x = xrep * vold
        v += x
        vsq += x * vold
        y = eold + self.pvcon*vold
        x = xrep * y
        h += x
        vh += x * vold
        if self.igbsa == 1:
            assert False
        else:
            self.E['esx'] += xrep * self.E['esonol']
            self.E['esxco'] += xrep * self.E['esonco']
            self.E['esxlj'] += xrep * self.E['esonlo']
            self.E['esx1'] += xrep * self.E['esol1']
            self.E['esx2'] += xrep * self.E['esol2']
        xztmp = 0.0
        for k in ['edihol','enbol','ebndol','eangol','esbnol','esanol','esnbol','esdiol','ebcold']:
            xztmp += self.E[k]
        einter = self.E['old'] - xztmp
        x = einter + self.pvcon*vold
        y1 = xrep * x
        hinter += y1
        hinsq += y1 * x
        if self.nofep == 0:
            x = self.E['exxold']+self.E['esonol']+self.E['edihol']+self.E['enbol']+self.E['ebndol']+self.E['eangol']+epoold
            del1 = self.E['exxol1']+self.E['esol1']+self.E['ediol1']+self.E['enbol1']+self.E['ebnol1']+self.E['eanol1']+self.E['epool1'] - x
            del2 = self.E['exxol2']+self.E['esol2']+self.E['ediol2']+self.E['enbol2']+self.E['ebnol2']+self.E['eanol2']+self.E['epool2'] - x
            x = math.exp(-self.beta*del1)*xrep
            edelg = edelg + x
            udel1 = udel1 + (y+del1)*x
            x = math.exp(-self.beta*del2)*xrep
            edelg2 = edelg2 + x
            udel2 = udel2 + (y+del2)*x
            x = math.exp(-self.beta*del1*0.5)*xrep
            sos1  = sos1 + x
            x = math.exp(-self.beta*del2*0.5)*xrep
            sos2  = sos2 + x

        if nsatm > 1:
            self.E['exx'] += xrep * self.E['exxold']
            self.E['exx1'] += xrep * self.E['exxol1']
            self.E['exx2'] += xrep * self.E['exxol2']
            self.E['epol'] += xrep * self.E['epoold']
            self.E['epol1'] += xrep * self.E['epool1']
            self.E['epol2'] += xrep * self.E['epool2']
            self.E['exxco'] += xrep * self.E['exxolc']
            self.E['exxlj'] += xrep * self.E['exxoll']
            self.E['ebc'] += xrep * self.E['ebcold']
            self.E['ebc1'] += xrep * self.E['ebcol1']
            self.E['ebc2'] += xrep * self.E['ebcol2']

        if self.ntdih != 0:
            self.E['edih'] += xrep * self.E['edihol']
            self.E['edih1'] += xrep * self.E['ediol1']
            self.E['edih2'] += xrep * self.E['ediol2']

        self.E['enb'] += xrep * self.E['enbol']
        self.E['enb1'] += xrep * self.E['enbol1']
        self.E['enb2'] += xrep * self.E['enbol2']

        if self.ntbnd != 0:
            self.E['ebnd'] += xrep * self.E['ebndol']
            self.E['ebnd1'] += xrep * self.E['ebnol1']
            self.E['ebnd2'] += xrep * self.E['ebnol2']

        if self.ntang != 0:
            self.E['eang'] += xrep * self.E['eangol']
            self.E['eang1'] += xrep * self.E['eanol1']
            self.E['eang2'] += xrep * self.E['eanol2']

        if self.natmx == 0:
            pass
            #goto 140
        
        for i in range(self.innerpars['nrdl']):
            for j in range(self.innerpars['nrdf']):
                xr[i][j] += xrep * ior[i][j]

            for j in range(self.innerpars['nrdfs']):
                xrs[i][j] += xrep * iors[i][j]
        
        i = int(self.dince*(eonold-self.innerpars['edmin'])) + 1
        if i >= 1 and i <= 50:
            xedist[i] += xrep
        
        for i in range(50):
            xepar[i] += xrep * ioepar[i]
            xess[i] += xrep * ioess[i]
        
        i = int(self.dinces*(self.E['esonol']-self.innerpars['ebsmin'])) + 1
        if i >= 1 and i <= 50:
            xbesol[i] += xrep
        
        for i in range(len(phi)):
            j = int(phi[i]/60.0/dtorad) + 1
            if j > 60: j -= 60
            iphid[i][j] += nrepet
        
        naccept += nrepet
        if last == 0:
            pass
            #goto 690
        eold = enew
        vold = vnew
        epoold = eponew
        epool1 = epone1
        epool2 = epone2
        nrepet = 1

        if movsol != 1:
            if movvol != 1:
                self.irn, y = randu(self.irn)
                nmr = int(nmol+1.0)
                movacp[nmov] = movacp[nmov] + 1
                esonol = esone
                esonco = esonc
                esonlo = esonl
                esol1 = eson1
                ess1[nmov] = es1
                esol2 = eson2
                ess2[nmov] = es2

                j = int(epins*(ess(nmov)-essmin)) + 1
                k = int(epins*(esol-essmin)) + 1

                ess[nmov] = esol
                essc[nmov] = ecsx
                essl[nmov] = elsx
                if j > 1 and j <= 50:
                    ioess[j] = ioess[j] - 1
                if k > 1 and k <= 50:
                    ioess[k] = ioess[k] + 1

                ioepar = [0 for i in range(50)]
                eonold = 0.0

                for i in range(nmov):
                    n = nmol*(i-1) + nmov - i - (i*(i-1))//2
                    eij[n] = emov[i]

                j = nmol*(nmov-1) - nmov - (nmov*(nmov-1))//2
                for i in range(nmov+1, nmol):
                    n = j+i
                    eij[n] = emov[i]

                for i in range(1, NMR-1):
                    n = nmol*(i-1) + nmr - i - (i*(i-1))//2
                    e = eij[n]
                    eonold = eonold+e
                    j = int(epin*(e-eprmin)) + 1
                    if j > 1 and j <= 50:
                        ioepar[j] += 1

                k = nmol*(nmr-1) - nmr - (nmr*(nmr-1))//2
                for i in range(nmr+1, nmol):
                    e = eij[i+k]
                    eonold += e
                    j = int(epin*(e-eprmin)) + 1
                    if j > 1 and j <= 50:
                        ioepar[j] += 1
                    
                ntyp = nstyp[nmov]
                natoms = nsvat[ntyp]
                
                """
                do 240 j = 1, 3
                do 230 k = 1, natoms
                    x = anew(k,j)
                    anew(k,j) = ac(nmov,k,j)
                    ac(nmov,k,j) = x
    230          continue
    240       continue
                """

                if ntyp == 2:
                    pass
                    #go to 260
                if nrdfs == 0:
                    pass
                    #go to 260

                movtyp = 2
                if modsv1 <= 2:
                    x = wxpot(0)
                if modsv1 > 2:
                    x = sxpot(0)

                movtyp = 0
                for i in range(1, nrdfs):
                    m = knew(i)
                    if m <= nrdl and m >= 1:
                        iors[m][i] += 1
                    m = int(rinc*(rnew[i]-rdlmin)) + 1
                    if m <= nrdl and m >= 1:
                        iors[m][i] -= 1
                
                #260
                if nbuse == 1:
                    nbmov[nmov] += 1
                    if nbmov[nmov] == nblim:
                        pass
                        # call neibor
                
                if nrdf == 0:
                    pass
                    #goto 530
                
                #270
                while True:
                    ntyp = nstyp[nmr]
                    if ntyp == 2:
                        self.irn, y = randu(self.irn)
                        nmr = int(xmol*y+1.0)
                    else:
                        break

                """
                    NATOMS = NSVAT(1)
                    NCEN = 1
                    IF (MODCUS.EQ.1) NCEN = NCENTS
                    DO 290 I = 1, 3
                    DO 280 J = 1, NATOMS
                        ANEW(J,I) = AC(NMR,J,I)
        280          CONTINUE
        290       CONTINUE
                    DO 310 I = 1, NRDF
                    DO 300 J = 1, NRDL
                        IOR(J,I) = 0
        300          CONTINUE
        310       CONTINUE
                """












    def movmol(self):
        pass



    def ssljco(self):
        pass
    
    def movsvn(self):
        pass

    def calc_solutes_dipole(self):
        pass








    def calc_neighbor(self):
        edge = self.solventsdata['edge']
        wxyz = self.solventsdata['xyz']
        nbmov = [0 for i in range(self.nmol)]
        nbor = {}
        for i in range(self.nmol-1):
            n1 = 1 if i in self.nmr else 0
            for j in range(i+1,self.nmol):
                n2 = 1 if j in self.nmr else 0
                bo = False
                for k in range(self.ncutas[n1]):
                    nci = self.icutas[n1][k] - 1
                    for l in range(self.ncutas[n2]):
                        ncj = self.icutas[n2][l] - 1
                        dx = abs(wxyz[i][nci][0] - wxyz[j][ncj][0])
                        if dx > edge[0]: dx -= edge[0]*2.0
                        dy = abs(wxyz[i][nci][1] - wxyz[j][ncj][1])
                        if dy > edge[1]: dy -= edge[1]*2.0
                        dz = abs(wxyz[i][nci][2] - wxyz[j][ncj][2])
                        if self.islab != 1:
                            if dz > edge[1]: dz -= edge[2]*2.0
                        if dx*dx+dy*dy+dz*dz > self.rncsq:
                            continue
                        bo = True
                        nbmov[i] += 1
                        nbmov[j] += 1
                        nbor[f'{i}-{nbmov[i]}'] = j
                        nbor[f'{j}-{nbmov[j]}'] = i
                    if bo: break
                if bo: break
        self.nbor = nbor

    def calc_hydrogen_da(self,key=None):
        # need more debug -- hydrogen donor & acceptor
        if not key: key = 'reference'
        xyz = self.solutesdata[key]['xyz']
        atomtypes = self.solutesdata[key]['atomtypes']
        a = self.solutesdata[key]['A']
        self.numhdonor = 0
        self.numhacceptor = 0
        for rxyz,st,a in zip(xyz,atomtypes,a):
            #                                 heteroatom only if hydrogen
            if st in [7,8,16] or (st == 1 and a <= 0.0):
                for i,wxyz in enumerate(self.solventsdata['xyz']):
                    if i in self.nmr:
                        bp = self.solventsdata['2']['bond_pars']
                    else:
                        bp = self.solventsdata['1']['bond_pars']
                    sn = len(bp)
                    for j in range(len(bp)):
                        sp = bp[j][2]
                        mxyz = xyz[j]
                        if a == 1:
                            if sp in [7,8,16]:
                                d = sum([rxyz[t]-mxyz][t] for t in range(3))
                                if d < 6.25:
                                    self.numhdonor += 1
                        elif sp == 1:
                            d = sum([rxyz[t]-mxyz][t] for t in range(3))
                            if d < 6.25:
                                self.numhacceptor += 1






