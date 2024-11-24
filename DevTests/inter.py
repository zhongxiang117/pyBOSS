def hbond(nsatm, iztyp, asol, iatno, nmol, nsvat, nstyp, isatno, ac, aw, a, ityp):
    rhb2 = 6.25
    naccp = 0
    ndon = 0
    for i in range(nsatm):
        if iztyp[i] <= 0:
            continue
        j = iatno[i]
        x = asol[i][0]
        y = asol[i][1]
        z = asol[i][2]
        bo = False
        if j == 1:
            if a[ityp[i]] == 0:
                bo = True
        if j in [7, 8]:
            bo = True
        if not (bo or j==16):
            continue
        for l in range(nmol):
            ltyp = nstyp[l]
            for k in range(nsvat[ltyp]):
                m = isatno[k][ltyp]
                if j == 1:
                    if m in [7, 8, 16]:
                        r = (x - ac[l][k][0])**2 + (y - ac[l][k][1])**2 + (z - ac[l][k][2])**2
                        if r <= rhb2:
                            ndon += 1
                elif m == 1:
                    if aw[k][ltyp] != 0:
                        continue
                    r = (x - ac[l][k][0])**2 + (y - ac[l][k][1])**2 + (z - ac[l][k][2])**2
                    if r <= rhb2:
                        naccp += 1
    return naccp, ndon


def ssljco(nmol, noss, nstyp, modsv1, modsv2, anew, nsvat, ac):
    esljol = 0.0
    escool = 0.0
    sslj = 0.0
    if nmol == 0:
        return esljol, sslj, escool
    if noss == 1:
        return esljol, sslj, escool
    movtyp = 0
    n = nmol - 1
    for i in range(n):
        k = i + 1
        nmov = i
        jtyp = nstyp[i]
        natoms = nsvat[jtyp]
        ntyp = modsv1
        if jtyp == 2:
            ntyp = modsv2
        for l in range(natoms):
            for m in range(3):
                anew[l][m] = ac[i][l][m]
        e = 0.0
        for j in range(k,nmol):
            if nsvat[1] == 0:
                if modsv1 > 2:
                    e = sspot(j)
                else:
                    e = wwpot(j)
            else:
                jtyp = modsv1
                if nstyp[j] == 2:
                    jtyp = modsv2
                if ntyp <= 2 and jtyp <= 2:
                    e = wwpot(j)
                else:
                    e = sspot(j)
            esljol += sslj
            escool += e - sslj
    esljol += ecut()
    return esljol, sslj, escool


def sspot(
    nm, edge, islab, nstyp, nsvat, ac, nmov, ncutas, icutas, anew, edg2, rcutsq, nmol2,
    isvaty, qw, aw, bw, iewald, irfon, rffac, natmx, aas, bbs, qqs, rl2, rul, ru2
):
    sspot = 0.0
    sslj = 0.0
    xedg = edge[0]
    yedg = edge[1]
    zedg = edge[2]
    
    if islab == 1:
        zedg = 0.0

    n1 = nstyp[nmov]
    n2 = nstyp[nm]
    nat2 = nsvat[n2]
    
    c = [[0.0,0.0,0.0] for i in range(nat2)]
    for i in range(3):
        for j in range(nat2):
            c[j][i] = ac[nm][j][i]

    wik = 1.0e20
    for i in range(ncutas[n1]):
        nw = icutas[n1][i]
        for j in range(ncutas[n2]):
            nc = icutas[n2][j]
            xiz = 0.0
            zz = abs(anew[nw][2] - c[nc][2])
            if zz > edg2[2]:
                xiz = zedg
            twik = (zz - xiz) ** 2
            if twik > rcutsq:
                continue
            xix = 0.0
            xx = abs(anew[nw][0] - c[nc][0])
            if xx > edg2[0]:
                xix = xedg
            twik += (xx - xix) ** 2
            if twik > rcutsq:
                continue
            xiy = 0.0
            yy = abs(anew[nw][1] - c[nc][1])
            if yy > edg2[1]:
                xiy = yedg
            twik += (yy - xiy) ** 2
            wik = min(wik, twik)

    if wik > rcutsq:
        return sspot

    if nmol2 != 0:
        for i in range(nsvat[n1]):
            if isvaty[i][n1] <= 0:
                continue

            x = anew[i][0]
            y = anew[i][1]
            z = anew[i][2]
            qi = qw[i][n1]
            ai = aw[i][n1]
            iflag = 1 if ai == 0.0 else 0
            bi = bw[i][n1]

            for j in range(nat2):
                if isvaty[j][n2] <= 0:
                    continue
                
                xx = abs(x - c[j][0])
                yy = abs(y - c[j][1])
                zz = abs(z - c[j][2])
                xix = xedg if xx > edg2[0] else 0.0
                xiy = yedg if yy > edg2[1] else 0.0
                xiz = zedg if zz > edg2[2] else 0.0
                
                rr = (xx - xix) ** 2 + (yy - xiy) ** 2 + (zz - xiz) ** 2
                if iflag != 1:
                    r6 = 1.0 / (rr ** 3)
                    aaij = ai * aw[j][n2]
                    bbij = bi * bw[j][n2]
                    xlj = (aaij*r6 - bbij) * r6
                    sspot += xlj
                    sslj += xlj
                
                qqij = qi * qw[j][n2]
                if qqij != 0.0:
                    if iewald == 1:
                        raise NotImplementedError('Ewald not supported')
                    else:
                        r1 = 1.0 / pow(rr,0.5)
                        if irfon == 0:
                            sspot += qqij * r1
                        else:
                            sspot += qqij * (r1 + rr * rffac)
    else:
        k = 0
        for i in range(natmx):
            if isvaty[i][0] <= 0:
                continue

            x = anew[i][0]
            y = anew[i][1]
            z = anew[i][2]

            for j in range(natmx):
                if isvaty[j][0] <= 0:
                    continue
                
                xx = abs(x - c[j][0])
                yy = abs(y - c[j][1])
                zz = abs(z - c[j][2])
                xix = xedg if xx > edg2[0] else 0.0
                xiy = yedg if yy > edg2[1] else 0.0
                xiz = zedg if zz > edg2[2] else 0.0
                
                rr = (xx - xix) ** 2 + (yy - xiy) ** 2 + (zz - xiz) ** 2
                if aas[k] != 0.0:
                    r6 = 1.0 / (rr ** 3)
                    xlj = (aas[k] * r6 - bbs[k]) * r6
                    sspot += xlj
                    sslj += xlj

                if qqs[k] != 0.0:
                    if iewald == 1:
                        raise NotImplementedError('Ewald not supported')
                    else:
                        r1 = 1.0 / pow(rr,0.5)
                        if irfon == 0:
                            sspot += qqs[k] * r1
                        else:
                            sspot += qqs[k] * (r1 + rr*rffac)
                k += 1

    if wik < rl2:
        return sspot
    
    feat = rul * (ru2 - wik)
    sspot *= feat
    sslj *= feat

    return sspot


def wwpot(
    nm, anew, ac, islab, edg2, edge, rcutsq, aoo, coo, natoms, irfon, rffac, qq, rl2, ru2, rul,
    kq, iewald
):
    xiz = 0.0
    zz = abs(anew[0][2] - ac[nm][0][2])
    if islab != 1:
        if zz > edg2[2]:
            xiz = edge[2]
    wik = (zz - xiz) ** 2

    if wik > rcutsq:
        return 0.0

    xix = 0.0
    xx = abs(anew[0][0] - ac[nm][0][0])
    if xx > edg2[0]:
        xix = edge[0]
    wik += (xx - xix) ** 2

    if wik > rcutsq:
        return 0.0

    xiy = 0.0
    yy = abs(anew[0, 1] - ac[nm, 0, 1])
    if yy > edg2[1]:
        xiy = edge[1]
    wik += (yy - xiy) ** 2

    if wik > rcutsq:
        return 0.0

    r6 = 1.0 / (wik * wik * wik)
    val = sslj = r6 * (aoo * r6 - coo)

    anm = [[0.0,0.0,0.0] for i in range(natoms)]
    for i in range(3):
        for j in range(natoms):
            anm[j][i] = ac[nm][j][i]

    knt = 0
    for i in range(kq, natoms):
        for j in range(kq, natoms):
            xx = abs(anew[i][0] - anm[j][0])
            yy = abs(anew[i][1] - anm[j][1])
            zz = abs(anew[i][2] - anm[j][2])
            rr = (xx - xix) ** 2 + (yy - xiy) ** 2 + (zz - xiz) ** 2
            if iewald == 1:
                raise NotImplementedError('Ewald not supported')
            else:
                r1 = 1.0 / pow(rr,0.5)
                if irfon == 0:
                    val += qq[knt] * r1
                else:
                    val += qq[knt] * (r1 + rr * rffac)
            knt += 1
    if wik < rl2:
        return val

    val *= rul * (ru2 - wik)
    return val


def ecut(natmx, nocut, noss, rcut, nmol, nsvat, modsv1, nmol2, modsv2, aw, bw, vnew, TWOPI):
    val = 0.0

    if natmx == 0 or nocut == 1 or noss == 1:
        return val

    rc3 = 1.0 / (rcut ** 3)
    rc9 = rc3 ** 3
    xn = float(nmol)
    xm = 0.0

    if nsvat[1] == 0:
        if modsv1 <= 2 or modsv1 == 13:
            return val
    else:
        xn = float(nmol - nmol2)
        xm = float(nmol2)
        if modsv1 <= 2:
            xn = 0.0
        if modsv2 <= 2:
            xm = 0.0

    for i in range(natmx):
        ai = aw[i][0]
        bi = bw[i][0]
        aj = aw[i][1]
        bj = bw[i][1]

        for j in range(natmx):
            val += (TWOPI * xn * xn / (3.0 * vnew)) * (ai * aw[j][0] * rc9 / 3.0 - bi * bw[j][0] * rc3)
            val += (TWOPI * xm * xm / (3.0 * vnew)) * (aj * aw[j][1] * rc9 / 3.0 - bj * bw[j][1] * rc3)
            val += (TWOPI * xn * xm / (3.0 * vnew)) * (ai * aw[j][1] * rc9 / 3.0 - bi * bw[j][1] * rc3) * 2.0

    return val


def isolut(iatom, nsolut, nlatm):
    for i in range(nsolut):
        if iatom <= nlatm[i]: break
    return i


def cappot(r2, capfk, capsq5, caprad):
    if r2 <= capsq5: return 0.0
    val = capfk * (pow(r2,0.5) - caprad - 5) ** 2
    return val


def sxpot(
    nm, movtyp, nstyp, nsvat, ac, nmov, nsatm, icut, ncutas, icutas, anew, ncent1, ncent2,
    asol, islab, edg2, edge, nsolut, icutat, iztyp, nrdfs, scutsq, sl2, su2, sul, nsatno, nofep,
    a, b, q, ityp, ityp1, ityp2,  anew1, anew2, aw, bw, qw, irfon, iewald, rffac, isolec,
    nrdfa1, nrdfa2, asol1, asol2, capfk, icapat
):
    if movtyp != 2:
        n1 = nstyp[nm]
        natoms = nsvat[n1]
        anm = [[0.0,0.0,0.0] for i in range(natoms)]
        for i in range(3):
            for j in range(natoms):
                anm[i][j] = ac[nm][j][i]
    else:
        n1 = nstyp[nmov]
        natoms = nsvat[n1]

    iyes = [1 for i in range(nsatm)]
    if icut < 2:
        wik = 1e20
        for j in range(ncutas[n1]):
            twik = 0.0
            nc = icutas[n1][j]
            for i in range(3):
                cnew = anew[nc][i]
                cnm = anm[i][nc]
                if movtyp == 1:
                    cnew = 0.5 * (anew[ncent1][i] + anew[ncent2][i])
                elif movtyp == 2:
                    cnm = 0.5 * (asol[ncent1][i] + asol[ncent2][i])
                x = abs(cnew - cnm)
                if not(islab == 1 and i == 2):
                    if x > edg2[i]:
                        x -= edge[i]
                twik += x**2
            wik = min(wik, twik)
    else:
        rmin = [1e20 for i in range(nsolut)]
        wik = 1e20
        for i in range(nsatm):
            bo = True
            if icut >= 3:
                if i in icutat:
                    bo = False
            else:
                if iztyp[i] == -1:
                    continue
            if bo: continue

            m = isolut(i)
            for l in range(ncutas[n1]):
                nc = icutas[n1][l]
                x = 0.0
                for j in range(3):
                    if movtyp != 2:
                        y = abs(anew[i][j] - anm[j][nc])
                    else:
                        y = abs(anew[nc][j] - asol[i][j])
                    if not(islab == 1 and j == 2):
                        if y > edg2[j]:
                            y -= edge[j]
                    x += y**2
                rmin[m] = min(rmin[m], x)
                wik = min(wik, x)

    rnew = [100.0 for i in range(nrdfs)]
    elj = [0.0, 0.0, 0.0]
    ecoul = [0.0, 0.0, 0.0]
    ec1 = ec2 = ec3 = 0.0
    ecsx = elsx = 0.0

    if wik > scutsq:
        return 0.0

    scal = 1.0
    if wik >= sl2:
        scal = sul * (su2 - wik)

    if icut != 0:
        j = 0
        for i in range(nsolut):
            if rmin[i] > scutsq:
                for k in range(j, j+nsatno[i]):
                    iyes[k] = 0
            j += nsatno[i]

    k0 = 0
    xedg, yedg, zedg = edge[0], edge[1], edge[2]
    if islab == 1:
        zedg = 0.0

    llim = 3
    if nofep == 1:
        llim = 1

    if movtyp != 2:
        for i in range(nsatm):
            if iztyp[i] == -1 or iyes[i] == 0:
                continue
            
            m = isolut(i)

            for l in range(1, llim + 1):
                xix = xiy = xiz = 0.0
                if l == 1:
                    ii = ityp[i]
                    xx, yy, zz = anew[i][0], anew[i][1], anew[i][2]
                elif l == 2:
                    ii = ityp1[i]
                    xx, yy, zz = anew1[i][0], anew1[i][1], anew1[i][2]
                elif l == 3:
                    ii = ityp2[i]
                    xx, yy, zz = anew2[i][0], anew2[i][1], anew2[i][2]

                aii, bii, qii = a[ii], b[ii], q[ii]

                for j in range(natoms):
                    x = abs(xx - anm[0][j])
                    y = abs(yy - anm[1][j])
                    z = abs(zz - anm[2][j])

                    if j == 0:
                        if x > edg2[0]:
                            xix = xedg
                        if y > edg2[1]:
                            xiy = yedg
                        if z > edg2[2]:
                            xiz = zedg

                    rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                    r1 = pow(rr,0.5)
                    r6 = 1.0 / (rr * rr * rr)
                    el = scal * r6 * (aii * aw[j][n1] * r6 - bii * bw[j][n1])
                    elj[l-1] += el
                    ec = scal * qii * qw[j][n1]

                    if irfon == 0:
                        if iewald == 1:
                            raise NotImplementedError('Ewald not supported')
                        else:
                            ec /= r1
                    else:
                        ec *= (1.0 / r1 + rr * rffac)
                    ecoul[l-1] += ec

                    if l == 1 and m == isolec:
                        ecsx += ec
                        elsx += el

                    for kk in range(k0, nrdfs):
                        if i == nrdfa1[kk]:
                            if j == nrdfa2[kk]:
                                rnew[kk] = r1
                                k0 = kk + 1
                                break
    else:
        for i in range(nsatm):
            if iztyp[i] == -1 or iyes[i] == 0:
                continue

            m = isolut(i)

            for l in range(1, llim + 1):
                xix = xiy = xiz = 0.0
                if l == 1:
                    ii = ityp[i]
                    xx = asol[i][0]
                    yy = asol[i][1]
                    zz = asol[i][2]
                elif l == 2:
                    ii = ityp1[i]
                    xx = asol1[i][0]
                    yy = asol1[i][1]
                    zz = asol1[i][2]
                elif l == 3:
                    ii = ityp2[i]
                    xx = asol2[i][0]
                    yy = asol2[i][1]
                    zz = asol2[i][2]

                aii = a[ii]
                bii = b[ii]
                qii = q[ii]

                for j in range(natoms):
                    x = abs(xx - anew[j][0])
                    y = abs(yy - anew[j][1])
                    z = abs(zz - anew[j][2])

                    if j == 1:
                        if x > edg2[0]:
                            xix = xedg
                        if y > edg2[1]:
                            xiy = yedg
                        if z > edg2[2]:
                            xiz = zedg

                    rr = (x - xix) ** 2 + (y - xiy) ** 2 + (z - xiz) ** 2
                    r1 = pow(rr,0.5)
                    r6 = 1.0 / (rr * rr * rr)

                    el = scal * r6 * (aii * aw[j][n1] * r6 - bii * bw[j][n1])
                    elj[l] += el
                    ec = scal * qii * qw[j][n1]

                    if irfon == 0:
                        if iewald == 1:
                            raise NotImplementedError('Ewald not supported')
                        else:
                            ec /= r1
                    else:
                        ec *= (1.0 / r1 + rr * rffac)

                    ecoul[l] += ec

                    if l == 1:
                        if m == isolec:
                            ecsx += ec
                            elsx += el

                        for kk in range(k0, nrdfs + 1):
                            if i == nrdfa1[kk]:
                                if j == nrdfa2[kk]:
                                    rnew[kk] = r1
                                    k0 = kk + 1
                                    break

    if capfk == 0.0:
        return

    if icapat != 0:
        x = y = z = 0.0
        if movtyp == 1:
            for i in range(3):
                xx = ac[nm][0][i]
                x += (anew[icapat][i] - xx) ** 2
                y += (anew1[icapat][i] - xx) ** 2
                z += (anew2[icapat][i] - xx) ** 2
        else:
            for i in range(3):
                xx = anew[0][i]
                x += (asol[icapat][i] - xx) ** 2
                y += (asol1[icapat][i] - xx) ** 2
                z += (asol2[icapat][i] - xx) ** 2

        ec1 = cappot(x)

        if nofep == 0:
            ec2 = cappot(y)
            ec3 = cappot(z)

    # Final calculations
    sxpot = elj[0] + ecoul[0] + ec1

    if nofep == 0:
        es1 = elj[1] + ecoul[1] + ec2
        es2 = elj[2] + ecoul[2] + ec3
    else:
        es1 = sxpot
        es2 = sxpot
    return


def wxpot(
    nm, movtyp, nstyp, nsvat, ac, nmov, modsv1, modsv2, nsatm, icut, anew, ncent1, ncent2, asol, islab, edg2,
    ncutat, icutat, iztyp, nrdfs, scutsq, nsatno, sl2, sul, su2, nofep, ityp, ityp1, ityp2, anew1, anew2,
    a, b, q, aw, bw, isolec, qw, irfon, iewald, rffac, nrdfa1, nrdfa2, asol1, asol2, edge, nsolut, capfk, icapat,
):
    if movtyp != 2:
        n1 = nstyp[nm]
        natoms = nsvat[n1]
        anm = [[0.0,0.0,0.0] for i in range(natoms)]
        for i in range(3):
            for j in range(natoms):
                anm[i][j] = ac[nm][j][i]
    else:
        n1 = nstyp[nmov]
        natoms = nsvat[n1]

    modsv = modsv1
    if n1 == 2:
        modsv = modsv2

    iyes = [1 for i in range(nsatm)]
    if icut < 2:
        wik = 0.0
        for i in range(3):
            cnew = anew[0][i]
            cnm = anm[i][0]
            if movtyp == 1:
                cnew = 0.5 * (anew[ncent1][i] + anew[ncent2][i])
            elif movtyp == 2:
                cnm = 0.5 * (asol[ncent1][i] + asol[ncent2][i])

            x = abs(cnew - cnm)
            if not (islab == 1 and i == 2):
                if x > edg2[i]:
                    x -= edge[i]
            wik += x**2
    else:
        rmin = [1e20 for i in range(nsolut)]
        wik = 1e20
        for i in range(nsatm):
            if icut >= 3:
                bo = True
                for k in range(ncutat):
                    if i == icutat[k]:
                        bo = False
                        break
                if bo: continue
            else:
                if iztyp[i] == -1:
                    continue

            x = 0.0
            for j in range(3):
                if movtyp != 2:
                    y = abs(anew[i][j] - anm[j][0])
                else:
                    y = abs(anew[0][j] - asol[i][j])

                if not (islab == 1 and j == 2):
                    if y > edg2[j]:
                        y -= edge[j]
                x += y**2

            j = isolut(i)
            rmin[j] = min(rmin[j], x)
            wik = min(wik, x)

    rnew = [100.0 for i in range(nrdfs)]

    ec1 = 0.0
    ec2 = 0.0
    ec3 = 0.0
    ecsx = 0.0
    elsx = 0.0

    elj = [0.0, 0.0, 0.0]
    ecoul = [0.0, 0.0, 0.0]

    if wik > scutsq:
        return

    scal = 1.0
    if wik >= sl2:
        scal = sul * (su2 - wik)

    if icut != 0:
        j = 0
        for i in range(nsolut):
            if rmin[i] > scutsq:
                l = j + nsatno[i]
                for k in range(j, l):
                    iyes[k] = 0
            j += nsatno[i]

    k0 = 0
    xedg = edge[0]
    yedg = edge[1]
    zedg = edge[2]
    if islab == 1:
        zedg = 0.0

    llim = 3
    if nofep == 1:
        llim = 1

    if movtyp != 2:
        for i in range(nsatm):
            if iztyp[i] == -1 or iyes[i] == 0:
                continue

            m = isolut(i)

            for l in range(llim):
                xix = xiy = xiz = 0.0
                if l == 0:
                    ii = ityp[i]
                    xx, yy, zz = anew[i][0], anew[i][1], anew[i][2]
                elif l == 1:
                    ii = ityp1[i]
                    xx, yy, zz = anew1[i][0], anew1[i][1], anew1[i][2]
                else:
                    ii = ityp2[i]
                    xx, yy, zz = anew2[i][0], anew2[i][1], anew2[i][2]

                aii, bii, qii = a[ii], b[ii], q[ii]

                for j in range(natoms):
                    x = abs(xx - anm[0][j])
                    y = abs(yy - anm[1][j])
                    z = abs(zz - anm[2][j])

                    bo = True
                    if j == 0:
                        if x > edg2[0]:
                            xix = xedg
                        if y > edg2[1]:
                            xiy = yedg
                        if z > edg2[2]:
                            xiz = zedg
                        rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                        r1 = pow(rr,0.5)
                        r6 = 1.0 / (rr**3)
                        xlj = scal * r6 * (aii * aw[0, n1] * r6 - bii * bw[0, n1])
                        elj[l] += xlj
                        if l == 0 and m == isolec:
                            elsx += xlj
                        if modsv == 2:
                            #goto 140
                            bo = False
                    else:
                        rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                        r1 = pow(rr,0.5)

                    if bo:
                        coul = qii * qw[j, n1] * scal
                        if irfon == 0:
                            if iewald == 1:
                                raise NotImplementedError('`iewald==1`')
                            else:
                                coul /= r1
                        else:
                            coul = coul * (1.0/r1 + rr*rffac)

                        ecoul[l] += coul
                        if l == 0 and m == isolec:
                            ecsx += coul

                        if j > 1:
                            continue

                    if l == 1:
                        for kk in range(k0, nrdfs):
                            if i == nrdfa1[kk] and j == nrdfa2[kk]:
                                rnew[kk] = r1
                                k0 = kk + 1
    else:
        for i in range(nsatm):
            if iztyp[i] == -1 or iyes[i] == 0:
                continue

            m = isolut(i)

            for l in range(1, llim+1):
                xix = 0.0
                xiy = 0.0
                xiz = 0.0

                if l == 1:
                    ii = ityp[i]
                    xx = asol[i][0]
                    yy = asol[i][1]
                    zz = asol[i][2]
                elif l == 2:
                    ii = ityp1[i]
                    xx = asol1[i][0]
                    yy = asol1[i][1]
                    zz = asol1[i][2]
                elif l == 3:
                    ii = ityp2[i]
                    xx = asol2[i][0]
                    yy = asol2[i][1]
                    zz = asol2[i][2]

                aii = a[ii]
                bii = b[ii]
                qii = q[ii]

                for j in range(natoms):
                    x = abs(xx - anew[j][0])
                    y = abs(yy - anew[j][1])
                    z = abs(zz - anew[j][2])

                    bo = True
                    if j == 1:
                        if x > edg2[0]:
                            xix = xedg
                        if y > edg2[1]:
                            xiy = yedg
                        if z > edg2[2]:
                            xiz = zedg

                        rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                        r1 = rr**0.5
                        r6 = 1.0 / (rr**3)
                        xlj = scal * r6 * (aii * aw[0][n1] * r6 - bii * bw[0][n1])
                        elj[l-1] += xlj

                        if l == 1 and m == isolec:
                            elsx += xlj

                        if modsv == 2:
                            #goto 200
                            bo = False
                    else:
                        rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                        r1 = rr**0.5

                    if bo:
                        coul = qii * qw[j][n1] * scal

                        if irfon == 0:
                            if iewald == 1:
                                raise NotImplementedError('`iewald`')
                            else:
                                coul /= r1
                        else:
                            coul *= (1.0 / r1 + rr * rffac)

                        ecoul[l-1] += coul

                        if l == 1 and m == isolec:
                            ecsx += coul

                        if j >= 2:
                            continue

                    if l == 1:
                        for kk in range(k0, nrdfs):
                            if i == nrdfa1[kk] and j == nrdfa2[kk]:
                                rnew[kk] = r1
                                k0 = kk + 1
                                break

    if capfk == 0.0:
        return

    if icapat != 0:
        x = y = z = 0.0
        if movtyp == 1:
            for i in range(3):
                xx = ac[nm][0][i]
                x += (anew[icapat][i] - xx) ** 2
                y += (anew1[icapat][i] - xx) ** 2
                z += (anew2[icapat][i] - xx) ** 2
        else:
            for i in range(3):
                xx = anew[0][i]
                x += (asol[icapat][i] - xx) ** 2
                y += (asol1[icapat][i] - xx) ** 2
                z += (asol2[icapat][i] - xx) ** 2

        ec1 = cappot(x)

        if nofep == 0:
            ec2 = cappot(y)
            ec3 = cappot(z)

    wxpot = elj[0] + ecoul[0] + ec1

    if nofep == 0:
        es1 = elj[1] + ecoul[1] + ec2
        es2 = elj[2] + ecoul[2] + ec3
    else:
        es1 = wxpot
        es2 = wxpot


def xxpot(
    edge, islab, nbc, nbcat1, nbcat2, rbc, anew, bck0, bck1, bck2, nofep, anew1, anew2, 
    ntdih, phinew, dihv0, rc0, rc1, rc2, isqm, itydi
):
    exxnew = exxne1 = exxne2 = exxnec = exxnel = ebcnew = ebcne1 = ebcne2 = 0.0
    edgg = [edge[0], edge[1], edge[2]]
    
    if islab == 1:
        edgg[2] = 0.0
    
    if nbc != 0:
        for i in range(nbc):
            j = nbcat1[i]
            k = nbcat2[i]
            r = rbc[i]

            if k == 9999:
                rr1 = r * ((anew[j][0] - bck0[i])**2 + (anew[j][1] - bck1[i])**2 + (anew[j][2] - bck2[i])**2)
                ebcnew += rr1

                if nofep == 0:
                    rr2 = r * ((anew1[j][0] - bck0[i])**2 + (anew1[j][1] - bck1[i])**2 + (anew1[j][2] - bck2[i])**2)
                    rr3 = r * ((anew2[j][0] - bck0[i])**2 + (anew2[j][1] - bck1[i])**2 + (anew2[j][2] - bck2[i])**2)
                    ebcne1 += rr2
                    ebcne2 += rr3
            
            elif j == 9999:
                j1 = int(bck0[i])
                j2 = int(bck1[i])
                j3 = int(bck2[i])
                
                if j1 == 0:
                    continue
                
                x = anew[j1][0]
                y = anew[j1][1]
                z = anew[j1][2]
                fac = 1.0
                
                if j2 != 0:
                    x += anew[j2][0]
                    y += anew[j2][1]
                    z += anew[j2][2]
                    fac = 2.0
                
                if j3 != 0:
                    x += anew[j3][0]
                    y += anew[j3][1]
                    z += anew[j3][2]
                    fac += 1.0
                
                x /= fac
                y /= fac
                z /= fac
                
                rr1 = r * ((anew[k][0] - x)**2 + (anew[k][1] - y)**2 + (anew[k][2] - z)**2)
                ebcnew += rr1
                ebcne1 += rr1
                ebcne2 += rr1
            
            elif r < 0.0:
                rr1 = ((anew[j][0] - anew[k][0])**2 + (anew[j][1] - anew[k][1])**2 + (anew[j][2] - anew[k][2])**2)
                rr2 = bck1[i]**2
                
                if rr1 > rr2:
                    rr3 = bck2[i] * (rr1**0.5 - bck1[i])**2
                    ebcnew += rr3
                else:
                    rr2 = bck0[i]**2
                    if rr1 > rr2:
                        continue
                    rr3 = bck2[i] * (rr1**0.5 - bck0[i])**2

                ebcnew += rr3
            else:
                rr1 = ((anew[j][0] - anew[k][0])**2 + (anew[j][1] - anew[k][1])**2 + (anew[j][2] - anew[k][2])**2)
                ebcnew += bck0[i] * (rr1**0.5 - r)**2

                if nofep == 0:
                    rr2 = ((anew1[j][0] - anew1[k][0])**2 + (anew1[j][1] - anew1[k][1])**2 + (anew1[j][2] - anew1[k][2])**2)
                    rr3 = ((anew2[j][0] - anew2[k][0])**2 + (anew2[j][1] - anew2[k][1])**2 + (anew2[j][2] - anew2[k][2])**2)
                    ebcne1 += bck1[i] * (rr2**0.5 - r)**2
                    ebcne2 += bck2[i] * (rr3**0.5 - r)**2
    
    if nofep == 1:
        ebcne1 = ebcnew
        ebcne2 = ebcnew

    if ntdih != 0:
        for i in range(ntdih):
            if itydi[i] == 500:
                p = phinew[i]
                x = dihv0[i][0]
                e = dihv0[i][3] * (p - x) ** 2
                ebcnew += e
                if rc0 != rc1:
                    ebcne1 += e
                if rc0 != rc2:
                    ebcne2 += e

    if isqm == 1:
        return 
    raise NotImplementedError('isqm!=1')


def ranu(irn):
    icnt = 0
    ichg = 1167
    ichg0 = 1167
    imul = 1173
    icon = 458753759
    imod = 1048573
    icn0 = 458753759
    jmul = 1161
    jcon = 458716759
    jmod = 1048573
    jrn = 124690

    icnt += 1
    if icnt == ichg:
        jrn = (jrn * jmul + jcon) % jmod
        rnj = float(jrn) / float(jmod)
        if rnj > 0.5:
            fac = 1.0 + 0.5 * rnj
        else:
            fac = 1.0 - 0.5 * rnj
        fac = float(icn0) * fac
        icon = int(fac)

        jrn = (jrn * jmul + jcon) % jmod
        rnj = float(jrn) / float(jmod)

        if rnj > 0.5:
            fac = 1.0 + 0.5 * rnj
        else:
            fac = 1.0 - 0.5 * rnj
        fac = float(ichg0) * fac
        ichg = int(fac)
        icnt = 0

    irn = (irn * imul + icon) % imod
    return float(irn) / float(imod)


import math
import random  # Used for generating random numbers

def rotate(
    angle,
    movtyp, nstyp, nat, nmov, modcus, ncents, iflxsl, isol, nrota1, nrota2, nfatm, nvdih, nvbnd, nvang,
    nlatm, anew, anew1, anew2
):
    iax = int(3 * ranu() + 1)
    if movtyp != 1:
        iat = 1
        ifat = 2
        ilat = nat
        if nstyp[nmov] == modcus:
            iat = ncents
            ifat = 1
            if iflxsl == 1:
                ilat = 3
    else:
        if isol == 1:
            iat = nrota1
        elif isol == 2:
            iat = nrota2
        elif isol > 2:
            iat = nfatm[isol]

        ifat = nfatm[isol]
        ilat = nlatm[isol]
        n = nvdih + nvbnd + nvang
        if n != 0:
            if (ilat - ifat) > 2:
                ilat = ifat + 2

    csa = math.cos(angle)
    sna = math.sin(angle)

    if iax == 1:
        jax, kax = 1, 2
    elif iax == 2:
        jax, kax = 2, 0
    elif iax == 3:
        jax, kax = 0, 1

    x = anew[iat][jax]
    y = anew[iat][kax]

    for i in range(ifat, ilat):
        if i != iat:
            xo = anew[i][jax] - x
            yo = anew[i][kax] - y
            x1 = xo * csa - yo * sna
            y1 = xo * sna + yo * csa
            anew[i][jax] = x + x1
            anew[i][kax] = y + y1

    if movtyp != 0:
        x = anew1[iat][jax]
        y = anew1[iat][kax]
        for i in range(ifat, ilat):
            if i != iat:
                xo = anew1[i][jax] - x
                yo = anew1[i][kax] - y
                x1 = xo * csa - yo * sna
                y1 = xo * sna + yo * csa
                anew1[i][jax] = x + x1
                anew1[i][kax] = y + y1

        x = anew2[iat][jax]
        y = anew2[iat][kax]
        for i in range(ifat, ilat):
            if i != iat:
                xo = anew2[i][jax] - x
                yo = anew2[i][kax] - y
                x1 = xo * csa - yo * sna
                y1 = xo * sna + yo * csa
                anew2[i][jax] = x + x1
                anew2[i][kax] = y + y1


def movsvn(
    nsvat, nstyp, nmov, modcus, ncents, rdel, rdl2, anew, ac, islab, edg2, edge, adel, adl2, iflxsl,
    modsv1, modsv2, noss, nmol, nbuse, nbor, nrdfs, knew, rinc, rnew, rdlmin, eij,
    esonol, ess, esonco, esonlo, ecsx, essc, essl, elsx, es1, es2, ess1, ess2,
    esol1, esol2, eold, npols, icussl, iewald, nopref, wnew, wi, wisum, wik, wkc, 
):
    natoms = nsvat[nstyp[nmov]]
    nat = natoms
    nc = 1
    if nstyp[nmov] == modcus:
        nc = ncents

    for i in range(3):
        rr = -rdel + ranu() * rdl2
        for j in range(natoms):
            anew[j][i] = ac[nmov][j][i] + rr

        if islab == 1 and i == 2:
            continue

        if anew[nc][i] <= edg2[i]:
            if anew[nc][i] < -edg2[i]:
                for j in range(natoms):
                    anew[j][i] += edge[i]
        else:
            for j in range(natoms):
                anew[j][i] -= edge[i]

    achg = -adel + ranu() * adl2

    rotate(achg)

    if iflxsl == 1:
        if modcus == nstyp[nmov - 1]:
            raise NotImplementedError('iflxsl == 1')

    e2 = 0.0
    emov = [1 for i in range(2000)]

    ntyp = modsv1
    if nstyp[nmov - 1] == 2:
        ntyp = modsv2

    if noss != 1:
        for i in range(nmol):
            if nbuse == 1:
                j = nbor[nmov][i]
                if j == 0:
                    break
            else:
                if i == nmov:
                    continue
                j = i

            if nsvat[1] == 0:
                if modsv1 > 2:
                    e = sspot(j)
                else:
                    e = wwpot(j)
            else:
                jtyp = modsv1
                if nstyp[j] == 2:
                    jtyp = modsv2

                if ntyp <= 2 and jtyp <= 2:
                    e = wwpot(j)
                else:
                    e = sspot(j)

            emov[j] = e
            e2 += e

    movtyp = 2

    if noss == 0:
        if ntyp > 2:
            esol = sxpot(0)
        else:
            esol = wxpot(0)
    else:
        raise NotImplementedError('`sxnoss`')

    movtyp = 0

    if nrdfs != 0:
        for i in range(nrdfs):
            knew[i] = int(rinc * (rnew[i] - rdlmin)) + 1

    e1 = 0.0

    for i in range(nmov-1):
        e1 += eij[nmol*i + nmov - i - (i*(i-1)) // 2]

    j = nmol*nmov - nmov - (nmov * (nmov-1)) // 2
    for i in range(nmov, nmol):
        e1 += eij[j + i + 1]

    emov[nmov] = 0.0
    esone = esonol + esol - ess[nmov]
    esonc = esonco + ecsx - essc[nmov]
    esonl = esonlo + elsx - essl[nmov]
    eson1 = esol1 + es1 - ess1[nmov]
    eson2 = esol2 + es2 - ess2[nmov]
    enew = eold + e2 - e1 + esol - ess[nmov]

    if npols != 0:
        raise NotImplementedError('`polslv`')

    if icussl == 1:
        if nstyp[nmov] == modcus:
            raise NotImplementedError('`esvint`')

    if iewald == 1:
        raise NotImplementedError('`iewald`')

    if nopref == 0:
        wnsum = 0.0
        wnmax = 0.0
        for i in range(nmol):
            wnew[i] = wi[i] * wisum

        wnew[nmov] = 1.0 / (wik + wkc)

        for i in range(nmol):
            wnsum += wnew[i]

        delta = 1.0 / wnsum

        for i in range(nmol):
            wnew[i] *= delta
            wnmax = max(wnmax, wnew[i])


import random  # For generating random numbers

def movmol(
    movvol, nrdfs, nrdl, nsolut, nfatm, nlatm, icapat, anew, asol, asol1, asol2, anew1, anew2,
    ncent1, indsol, rdels, adels, ncent2, islab, edge, edg2, 
):
    if movvol == 1:
        # goto 300
        return

    idists = [ [0.0 for i in range(nrdl)] for j in range(nrdfs)]

    isol = 0
    if nsolut > 1:
        while True:
            isol = int(ranu() * nsolut) + 1
            if nfatm[isol] == nlatm[isol]:
                if icapat == nfatm[isol]:
                    continue
            break

    for j in range(3):
        for i in range(len(asol)):
            anew[i][j] = asol[i][j]
            anew1[i][j] = asol1[i][j]
            anew2[i][j] = asol2[i][j]

    # Determine displacement parameters
    icent = ncent1
    if isol <= 1 and indsol == 0:
        n = 1
        m = nlatm[1]
        if nsolut == 1:
            m = nlatm[0]
        rdl = rdels[0]
        adl = adels[0]
    else:
        rdl = rdels[isol]
        adl = adels[isol]
        if isol <= 1:
            icent = ncent2
        elif isol > 1:
            icent = nfatm[isol]
        n = nfatm[isol]
        m = nlatm[isol]

    for i in range(3):
        rr = -rdl + ranu() * (2 * rdl)
        for j in range(n, m):
            anew[j][i] = asol[j][i] + rr
            anew1[j][i] = asol1[j][i] + rr
            anew2[j][i] = asol2[j][i] + rr

        if islab == 1 and i == 2:
            continue

        if anew[icent][i] <= edg2[i]:
            if anew[icent][i] < -edg2[i]:
                delta = edge[i]
                for j in range(n, m):
                    anew[j][i] += delta
                    anew1[j][i] += delta
                    anew2[j][i] += delta
        else:
            delta = -edge[i]
            for j in range(n, m):
                anew[j][i] += delta
                anew1[j][i] += delta
                anew2[j][i] += delta


