import math


class PyBOSSCalc:
    """BOSS Calc for Solvent Accessible Surface Area calculation
    
    Note:
        I do not know how does those codes work, but they work, by literally
        transforming them to Python.
    """
    def __init__(
        self, radsa, asol, rsolv, amber_sasa_atomtypes=None,
        debug=None, *args, **kws
    ):
        self.amber_sasa_atomtypes = amber_sasa_atomtypes
        if debug:
            self.mydebuglist = []
        radsa = [i for i in radsa]      # deep copy
        n = len(radsa)
        for i in range(n):
            if radsa[i] > 0.0:
                radsa[i] += rsolv

        occlud = [False for i in range(n)]
        for i in range(n):
            if radsa[i] <= 0.0:
                occlud[i] = True

        pi = math.pi
        pii = pi * 2.0
        piv = pi * 4.0
        degr = 180.0 / pi
        tol = 0.0001
        tolrt = pi * pow(tol,0.5)
        tolsq = tol * tol

        sected = [False for i in range(n*n//2)]
        m = 0
        for i in range(n-1):
            ri = asol[i]
            for j in range(i+1,n):
                rj = asol[j]
                t = (ri[0]-rj[0])**2 + (ri[1]-rj[1])**2 + (ri[2]-rj[2])**2
                rij = pow(t,0.5)
                if abs(radsa[i]-radsa[j]) > rij-tol:
                    if radsa[i] >= radsa[j]:
                        occlud[j] = True
                    else:
                        occlud[i] = True
                elif radsa[i]+radsa[j] > rij+tol:
                    sected[m] = True
                m += 1

        m = 0
        for i in range(n-1):
            for j in range(i+1,n):
                if occlud[i] or occlud[j]:
                    sected[m] = False
                m += 1

        # make index starts from 1
        sected.insert(0,None)
        occlud.insert(0,None)
        # `radsa` and `asol` are outside reference
        radsa.insert(0,None)
        asol.insert(0,None)
        s = [0.0 for i in range(n+1)]
        v = [0.0 for i in range(n+1)]

        kmax = 31
        ik = [0.0 for i in range(kmax)]
        cs = [0.0 for i in range(kmax)]
        th = [[0.0 for i in range(kmax)] for j in range(2)]
        cut = [False for i in range(kmax)]


        # note: `##number` means done
        sa = vl = 0.0
        ntot = n
        ninx = 0
        for i in range(1,n+1):          # 20
            #tmp = input(f'\ncontinue? i = {i}')
            if occlud[i]:
                ##goto 19
                if occlud[i]:
                    si = 0.0
                    vi = 0.0
                #print(' si, vi, radsa[i] = ', si,vi,radsa[i])
                s[i] = si * radsa[i] * radsa[i]
                v[i] = (si+vi) * radsa[i]**3 / 3.0
                ri = radsa[i]
                if ri > 0.0:
                    ri = ri - rsolv
                continue

            si = piv
            vi = 0.0
            xi = asol[i][0] / radsa[i]
            yi = asol[i][1] / radsa[i]
            zi = asol[i][2] / radsa[i]
            nint = 0
            m = i - ntot
            boj = True
            for j in range(1,n+1):      # 18
                #tmp = input(f'\ncontinue? i = {i}, j = {j}')
                if j <= i:
                    m = m + ntot - j
                else:
                    m += 1
                #print(' -> i, j, m = ', i,j,m)
                if j == i or not sected[m]:
                    ##goto 18
                    continue

                rj = radsa[j] / radsa[i]
                t = (asol[i][0]-asol[j][0])**2 + (asol[i][1]-asol[j][1])**2 + (asol[i][2]-asol[j][2])**2
                rij = pow(t,0.5) / radsa[i]
                csi = 0.5 * (1.0+rij*rij-rj*rj) / rij
                rmi = 1.0 - csi*csi
                #print('> rj, rij, csi, rmi = ',rj,rij,csi,rmi)
                if rmi < 0.0:
                    ##goto 18
                    continue

                rmi = pow(rmi,0.5)
                csj = 0.5 * (rj*rj+rij*rij-1.0) / rj / rij
                if abs(csj) > 1.0:
                    ##goto 18
                    continue

                csj = math.cos(0.5*(pi+math.acos(csi)-math.acos(csj)))
                rmj = pow(1.0-csj*csj, 0.5)
                xm = xi + csi * (asol[j][0]/radsa[i]-xi) / rij
                ym = yi + csi * (asol[j][1]/radsa[i]-yi) / rij
                zm = zi + csi * (asol[j][2]/radsa[i]-zi) / rij
                
                a = self.orient(xm, ym, zm, xi, yi, zi)
                
                # True: normal    False: break    None: continue
                boj = True
            #18
                eclipsd = True if csi < 0.0 else False
                ninx = 0
                init = 0
                th0 = 0.0
                mk = i - ntot
                #print('ninx, init, th0, mk = ',ninx,init,th0,mk)
                for k in range(1,ntot+1):       # 9
                    #tmp = input(f'\n now for k = {k}')
                    if k <= i:
                        mk = mk + ntot - k
                    else:
                        mk += 1
                    #print(' mk = ',mk)
                    if k == i or k == j or not sected[mk]: continue

                    xk = asol[k][0] / radsa[i]
                    yk = asol[k][1] / radsa[i]
                    zk = asol[k][2] / radsa[i]

                    # update
                    xk,yk,zk = self.xformx(-1, a, xm, ym, zm, xk, yk, zk)

                    rk = (radsa[k]/radsa[i])**2 - zk*zk
                    rmk = pow(xk*xk+yk*yk, 0.5)
                    #print('after: xk, yk, zk = ', xk, yk, zk)
                    #print('rk, rmk = ',rk, rmk)
                    if debug:
                        self.mydebuglist.extend([rk,rmk])

                    if rk < tol:
                        if zk > 0.0: sected[mk] = False
                        continue
                    elif rmk > tol:
                        cmm = (rmi*rmi+rmk*rmk-rk) / (2.0*rmi*rmk)
                    elif rk > rmi*rmi:
                        cmm = -1.0
                    else:
                        cmm = 1.0
                    #print('k, rmk, zk, rk, cmm = ', k, rmk, zk, rk, cmm)
                    if cmm > 1.0-tol:
                        if zk > 0.0: sected[mk] = False
                        ##goto 9
                        continue
                    elif cmm < tol-1.0:
                        sected[m] = False
                        if csi*rmk > rmi*(csi+zk) or csi+zk > 0.0:
                            boj = None
                            ##goto 18
                            break
                        occlud[i] = True
                        ##goto 19
                        boj = False
                        if occlud[i]:
                            si = 0.0
                            vi = 0.0
                        #print(' si, vi, radsa[i] = ', si,vi,radsa[i])
                        s[i] = si * radsa[i] * radsa[i]
                        v[i] = (si+vi) * radsa[i]**3 / 3.0
                        ri = radsa[i]
                        if ri > 0.0:
                            ri = ri - rsolv
                        break

                    ninx += 1
                    if ninx > kmax:
                        print('Warning: ')      # error
                    
                    ik[ninx] = k
                    cs[ninx] = cmm
                    th[0][ninx] = math.acos(xk/rmk)
                    if yk < 0.0:
                        th[0][ninx] = pii - th[0][ninx]
                    th[1][ninx] = math.acos(cs[ninx])
                    #print(' -> ninx th[1] = ',ninx, th[0][ninx], th[1][ninx])
                    boinx = True
                    for inx in range(ninx-1,0,-1):     #  8
                        dnx = abs(th[0][inx]-th[0][ninx])
                        if dnx > pi: dnx = pii - dnx
                        if (th[1][inx]+th[1][ninx]) > (pii-dnx-tolrt):
                            sected[m] = False
                            boj = None
                            ##goto 18
                            break
                        elif th[1][inx] > th[1][ninx]+dnx-tolrt:
                            ninx = ninx - 1
                            boinx = False
                            ##goto 9
                            break
                        elif th[1][ninx] > th[1][inx]+dnx-tolrt:
                            for ix in range(inx,ninx):
                                th[0][ix] = th[0][ix+1]
                                th[1][ix] = th[1][ix+1]
                                cs[ix] = cs[ix+1]
                                ik[ix] = ik[ix+1]
                                if ix+1 == init: init = ix
                            ninx = ninx - 1
                    if boj is None: break
                    if not boinx: continue
                    #print(' -> th0 = ',th0)

                    if th[1][ninx] > th0:
                        th0 = th[1][ninx]
                        init = ninx
                    cut[ninx] = True
                #9
            #18
        #20
                if boj is None:
                    continue
                elif boj is False:
                    break

                vj = pi
                sj = pii * (1.0-csi)

                #print(' vj, sj, init = ', vj, sj, init)
                if debug:
                    self.mydebuglist.extend([vj,sj])
                if ninx == 0:
                    next = init
                    ##goto 17
                    onext = next
                    #1001
                    if vj*rmi*rmi < tol:
                        sected[m] = False
                        ##goto 18
                        continue
                    else:
                        nint += 1
                    si = si - sj
                    vi = vi + vj*rmi*rmi*(csi-csj*rmi/rmj)
                    #print(' -> si, vi = ', si, vi)
                    continue   ##18

                overlap = False
                for ix in range(1,ninx+1):      # 11
                    if not cut[ix]:
                        ##goto 11
                        continue
                    yn0 = 1.0
                    yn1 = -1.0
                    for inx in range(1,ninx+1): # 10
                        if inx == ix or not cut[inx]:
                            continue        # goto 10
                        thx = th[0][inx] - th[0][ix]
                        #print(' -> thx = ',thx)
                        if thx < -pi:
                            thx = thx + pii
                        elif thx > pi:
                            thx = thx - pii
                        if th[1][inx]+th[1][ix] > abs(thx)-tol:
                            if abs(thx) < pi-tolsq:
                                ynx = (cs[inx]-math.cos(thx)*cs[ix]) / math.sin(thx)
                                if thx > 0.0 and ynx < yn0:
                                    yn0 = ynx
                                elif thx < 0.0 and ynx > yn1:
                                    yn1 = ynx
                            else:
                                print('Error')
                    if yn1 >= yn0-tol:
                        cut[ix] = False
                    if ix == init and yn1 > -1.0:
                        overlap = True
                #11
                    #print('...ix, ik, cs', ix,ik[ix],cs[ix],th[0][ix],th[1][ix],yn0,yn1,cut[ix],degr)

                if not cut[init]: sected[m] = False
                if not sected[m]:
                    ##goto 18
                    continue
                
                sn0 = pow(1.0-cs[init]*cs[init],0.5)
                next = init
                onext = -1
            #18
                while True:
                    #13
                    next,nint,si,vi,sj,vj,sn0 = self.boss_13(
                        next,nint,si,vi,sj,vj,sn0,    onext,cut,
                        eclipsd,cs,th,csi,ninx,tol,tolrt,tolsq,pi,pii
                    )

                    if onext == next:
                        ##goto 1001
                        break

                    #17
                    onext = next
                    if next != init:
                        ##goto 13
                        pass
                    else:
                        break
            #18
        #20
                #1001
                if vj*rmi*rmi < tol:
                    sected[m] = False
                    ##goto 18
                    continue
                else:
                    nint += 1
                si = si - sj
                vi = vi + vj*rmi*rmi*(csi-csj*rmi/rmj)
                #print(' -> si, vi = ', si, vi)
            ##END 18
        #20
            if boj is False: continue
            if si < 0.0: si = 0.0
            #19
            if occlud[i]:
                si = 0.0
                vi = 0.0
            #print(' si, vi, radsa[i] = ', si,vi,radsa[i])
            s[i] = si * radsa[i] * radsa[i]
            v[i] = (si+vi) * radsa[i]**3 / 3.0
            ri = radsa[i]
            if ri > 0.0:
                ri = ri - rsolv
        ##END 20

        # later process, important!
        asol.pop(0)
        if debug:
            self.asol = asol
        self.s = s[1:]
        self.v = v[1:]

    def boss_13(self,
        next,nint,si,vi,sj,vj,sn0,  onext,cut,
        eclipsd,cs,th,csi,ninx,tol,tolrt,tolsq,pi,pii
    ):
        #13
        last = next
        next = last
        vk = sk = 0.0
        thx = pii
        for inx in range(1,ninx+1):     # 15
            if inx == last or not cut[inx]:
                continue
            dtx = th[0][inx] - th[0][last]
            #print(' --> dtx = ',dtx)
            if dtx < 0.0:
                dtx = dtx + pii
            if dtx < thx:
                thx = dtx
                next = inx
        if onext == next:
            ##goto 1001
            return next,nint,si,vi,sj,vj,sn0

        sn1 = pow(1.0-cs[next]*cs[next], 0.5)
        vk = th[1][next] - cs[next]*sn1
        sk = 0.0
        if abs(csi) > tol and abs(cs[next]) > tol:
            tmp = (1.0/cs[next]/cs[next]-1.0) / (1.0/cs[next]/cs[next]+1.0/csi/csi-1.0)
            sk = 2.0 * math.acos(pow(tmp,0.5)) - pi
        elif abs(csi) > tol:
            sk = -pi
        #print(' -> 1 -- sk = ',sk)
        if cs[next] < 0.0:
            sk = -pii - sk
        if eclipsd:
            sk = -sk
        sk = -2.0 * th[1][next] * csi - sk
        #print(' -> 2 -- sk = ',sk)
        if eclipsd and cs[next] < 0.0:
            sk = -sk - pii
        dtx = th[1][next] + th[1][last] - thx
        if dtx > tolrt:
            snx = math.sin(thx)
            cnx = math.cos(thx)
            ynx = (cs[next]-cnx*cs[last]) / snx
            xkn = cnx*cs[last] + snx*sn0
            ykn = snx*cs[last] - cnx*sn0
            xko = 0.5 * (cs[next]+xkn)
            yko = 0.5 * (sn1+ykn)
            css = xko*xko + yko*yko

            #print(' -> xkn, ykn, xko, yko, css = ', xkn, ykn, xko, yko, css)
            if css < tolsq or abs(csi) < tol:
                pass
            elif css < 1.0 - tol:
                if pi - dtx >= 0:
                    cmm = pow(css,0.5)
                else:
                    cmm = -pow(css,0.5)
                tmp = pow((1.0/css-1.0)/(1.0/css+1.0/csi/csi-1.0), 0.5)
                skm = 2.0*math.acos(tmp) - pi

                #print(' -> 1 -- skm = ',skm)
                if dtx > pi:
                    skm = -pii - skm
                if eclipsd:
                    skm = -skm
                skm = -dtx*csi - skm
                if eclipsd and dtx > pi:
                    skm = -skm - pii
                #print(' -> 2 -- skm = ',skm)
                if cs[next] >= 0:
                    xkm = 1.0
                else:
                    xkm = -1.0
                if cs[next]*cs[next] > tol:
                    xkm = 1.0/cs[next]/pow(1.0/cs[next]/cs[next]+1.0/csi/csi-1.0, 0.5)
                zkm = pow(1.0-xkm*xkm, 0.5)
                xkn = 1.0 if cs[last] >= 0 else -1.0
                #print(' -> 3 -- skm, zkn, xkn = ',skm, xkn)

                if cs[last]*cs[last] > tol:
                    xkn = 1.0/cs[last]/pow(1.0/cs[last]/cs[last]+1.0/csi/csi-1.0, 0.5)
                zkn = pow(1.0-xkn*xkn, 0.5)
                ykn = xkn * snx
                xkn = xkn * cnx
                xpo = 1.0 if cmm >= 0.0 else -1.0
                if cmm*cmm > tol:
                    xpo = 1.0/cmm/pow(1.0/cmm/cmm+1.0/csi/csi-1.0, 0.5)
                zpo = pow(1.0-xpo*xpo, 0.5)
                ypo = xpo * yko / cmm
                xpo = xpo * xko / cmm

                #print(' -> xpo, zpo, ypo = ', xpo, zpo, ypo)
                asa = math.acos(1.0-0.5*((xkm-xkn)**2+ykn*ykn+(zkm-zkn)**2))
                asb = math.acos(1.0-0.5*((xkm-xpo)**2+ypo*ypo+(zkm-zpo)**2))
                asc = math.acos(1.0-0.5*((xkn-xpo)**2+(ykn-ypo)**2+(zkn-zpo)**2))
                #print(' -> asa, asb, asc = ', asa, asb, asc)
                asa = pi - asa

                if cs[next] < 0.0:
                    asa = pi - asa
                    asb = pi - asb
                if cs[last] < 0.0:
                    asa = pi - asa
                    asc = pi - asc
                if cmm < 0.0:
                    asb = pi - asb
                    asc = pi - asc
                skk = asa + asb + asc - pi
                if eclipsd: skk = -skk
                sk = sk - skm - skk

            #print(' cnx, cs(next), sn1 = ', cnx, cs[next], sn1, snx, cnx)
            xkm = cnx*cs[next] + snx*sn1
            ykm = snx*cs[next] - cnx*sn1
            vk = vk - 0.5*(dtx+cs[last]*(ynx-sn0+ykm)-ynx*xkm)
            #print(' -> skm, skk, sk = ',skm, skk, sk)
            #print(' -> xkm, ykm, vk = ',xkm, ykm, vk)
        
        vj = vj - vk
        sj = sj - sk
        sn0 = sn1
        cut[next] = False
        #print(' vj, sj, sn0 = ', vj, sj, sn0)
        return next,nint,si,vi,sj,vj,sn0

    def xformx(self, iopt, a, xm, ym, zm, x, y, z):
        if iopt == -1:
            nx = (x-xm)*a[0][0] + (y-ym)*a[1][0] + (z-zm)*a[2][0]
            ny = (x-xm)*a[0][1] + (y-ym)*a[1][1] + (z-zm)*a[2][1]
            nz = (x-xm)*a[0][2] + (y-ym)*a[1][2] + (z-zm)*a[2][2]
        else:
            nx = x*a[0][0] + y*a[0][1] + z*a[0][2] + xm
            ny = x*a[1][0] + y*a[1][1] + z*a[1][2] + ym
            nz = x*a[2][0] + y*a[2][1] + z*a[2][2] + zm
        return [nx,ny,nz]

    def orient(self, x0, y0, z0, x, y, z):
        a0 = 0.0
        a1 = 1.0
        tol = 0.000001
        snx = a0
        csx = a1
        sny = a0
        csy = a1
        snz = a0
        csz = a1
        dx = x - x0
        dy = y - y0
        dz = z - z0
        d2 = pow(dx*dx+dz*dz, 0.5)
        if d2 > tol:
            sny = dx / d2
            csy = -dz / d2
        d3 = pow(dx*dx+dy*dy+dz*dz,0.5)
        if d3 > tol:
            snx = -dy / d3
            csx = d2 / d3
        a = [[0.0 for i in range(3)] for j in range(3)]
        a[0][0] = csy*csz - snx*sny*snz
        a[1][0] = -csx*snz
        a[2][0] = sny*csz + snx*csy*snz
        a[0][1] = csy*snz + snx*sny*csz
        a[1][1] = csx*csz
        a[2][1] = sny*snz - snx*csy*csz
        a[0][2] = -csx*sny
        a[1][2] = snx
        a[2][2] = csx*csy
        return a

    def calc_amber_sasa(self,amber_sasa_atomtypes=None):
        if not amber_sasa_atomtypes:
            amber_sasa_atomtypes = self.amber_sasa_atomtypes
        mv = [0.0, 0.0, 0,0, 0.0]
        ms = [0.0, 0.0, 0,0, 0.0]
        for t,tv,ts in zip(amber_sasa_atomtypes,self.v,self.s):
            mv[t-1] += tv
            ms[t-1] += ts
        sasa = {}
        for i,key in enumerate(['Hydrophobic','Hydrophilic','Aromatic','Weakly Polar']):
            sasa[key] = [ms[i],mv[i]]
        return sasa




