from .constants import SOLVENT_PARS


def init_pars(dpars):
    """In-place editions of input parameter

    Idea:
        line-by-line translation, to be maximumly compatible with BOSS source codes
    """

    dpars['pconv'] = pconv = 0.000014576
    dpars['boltz'] = boltz = 0.00198717

    rc0 = dpars['rc0']
    rc1 = dpars['rc1']
    rc2 = dpars['rc2']
    if abs(rc0-rc1)<0.0001 and abs(rc1-rc2)<0.0001:
        nofep = 1
    else:
        nofep = 0
    if dpars['icalc'] > 1:
        nofep = 1
    dpars['nofep'] = nofep

    tk = dpars['t'] + 273.15
    dpars['beta'] = 1.0 / (boltz*tk)
    dpars['pvcon'] = dpars['p'] * pconv
    dpars['tk'] = tk

    # solvent 1
    svmod1 = dpars['svmod1']
    if svmod1 and svmod1 in SOLVENT_PARS:
        sv1p = SOLVENT_PARS[svmod1]
    else:
        print(f'Fatal: not supported: svmod1: {svmod1}')
        exit()

    dpars['icussl'] = 0
    dpars['svmod1'] = sv1p[0]
    dpars['modsv1'] = modsv1 = sv1p[1]
    dpars['natom1'] = len(sv1p[3])
    dpars['igbsa'] = 1 if modsv1 == 99 else 0
    dpars['modcus'] = 1 if modsv1 == 20 else 0
    if modsv1 == 99:
        dpars['nvchg'] = 999999
        dpars['nschg'] = 1
        dpars['nrdf'] = 0
        dpars['nmol'] = 0
        dpars['ibox'] = 0
        dpars['svmod2'] = ''

    if dpars['scllj'] == 0.0 or dpars['noss'] == 0:
        dpars['scllj'] = 1.0
    if dpars['noss'] == 0 and modsv1 != 99:
        fdielec = 1.0
    if fdielec == 0.0: fdielec = 1.0
    dpars['dielec'] = fdielec

    # solvent 2
    svmod2 = dpars['svmod2']
    if svmod2:
        if svmod2 in SOLVENT_PARS:
            sv2p = SOLVENT_PARS[svmod2]
        else:
            print(f'Fatal: not supported: svmod2: {svmod2}')
            exit()

        dpars['svmod2'] = sv2p[0]
        dpars['modsv2'] = modsv2 = sv2p[1]
        dpars['natom2'] = len(sv2p[3])
        dpars['igbsa'] = 1 if modsv2 == 99 else 0

        if modsv2+modsv1 == 3:
            print('Fatal: do not use TIP3P and TIP4P simultaneously')
            exit()
        elif modsv2 == 20:
            if modsv1 == 20:
                print('Fatal: only one solvent can be flexible')
                exit()
            else:
                dpars['modcus'] = 2
    else:
        dpars['modsv2'] = 0
        dpars['natom2'] = 0

    if dpars['rsolv'] <= 0.0:
        rsolv = 1.4
        if modsv1 < 99:
            rsolv = 3.2
        if modsv1 < 3 or modsv1 == 13:
            rsolv = 1.4
        dpars['rsolv'] = rsolv

    # wait for update
    dpars['nvsdih'] = dpars['nvsang'] = dpars['nvsbnd'] = 0
    if dpars['nvsdih'] + dpars['nvsang'] + dpars['nvsbnd'] != 0:
        dpars['iflxsl'] = 1
    else:
        dpars['iflxsl'] = 0

    dpars['nsvat'] = [dpars['natom1'],dpars['natom2']]
    dpars['natmx'] = max(dpars['nsvat'])

    dpars['nopref'] = 0
    if dpars['noss'] != 0:
        dpars['noss'] = 1
        dpars['nrdf'] = 0
        dpars['nrdfs'] = 0
        dpars['nofep'] = 1
        dpars['nopref'] = 1
        dpars['wkc'] = 20000.0

    qmname = dpars['qmname']
    if qmname in ['am1','am1s']:
        iampm = 1
    elif qmname in ['pm3','pm3s']:
        iampm = 2
    elif qmname in ['pdg','pdgs']:
        iampm = 3
    elif qmname in ['mdg','mdgs']:
        iampm = 4
    elif qmname in ['g09','g09s']:
        iampm = 5
    elif qmname in ['','none']:
        iampm = 0
    else:
        iampm = -1
        print(f'Fatal: unavailable QM choice: QMNAME: {qmname}')
        exit(-1)
    dpars['iampm'] = iampm
    if iampm != 0:
        dpars['isqm'] = 1
    else:
        dpars['isqm'] = 0

    cmname = dpars['cmname']
    if cmname in ['mull']:
        ichmod = 2
    elif cmname in ['cm1','cm1a','cm1p']:
        ichmod = 3
    elif cmname in ['cm3','cm3a','cm3p']:
        ichmod = 4
    elif cmname in ['gpar']:
        ichmod = 5
        iampm = 5
    elif cmname in ['','opls','none']:
        ichmod = 0
    else:
        ichmod = -1
        print(f'Fatal: unavailable QM choice: CMNAME: {cmname}')
        exit(-1)
    dpars['ichmod'] = ichmod
    if dpars['isqm'] == 0: dpars['ichmod'] = 0

    if dpars['qmscale'] <= 0.0:
        qmscale = 1.2
        if dpars['rc0'] != 0.0:
            qmscale = 1.0
        dpars['qmscale'] = qmscale

    dpars['nsatno'] = [0 for i in range(25)]
    dpars['nfatm'] = [0 for i in range(25)]
    dpars['nlatm'] = [0 for i in range(25)]
    dpars['nsatcz'] = [0 for i in range(25)]
    dpars['nfatm'][0] = 1

    dpars['nsolut'] = 1
    if dpars['isolec'] <= 0:
        dpars['isolec'] = 1

    slfmt = dpars['slfmt']
    if slfmt in ['zmat']:
        islfmt = 0
    elif slfmt in ['pdb','pdbp']:
        islfmt = 1
    elif slfmt in ['mind']:
        islfmt = 2
    elif slfmt in ['in']:
        islfmt = 3
    elif slfmt in ['zin']:
        islfmt = 4
    else:
        islfmt = 99
        print('*** UNRECOGNIZED SOLUTE FORMAT IN THE PAR FILE ***')
        exit()
    dpars['islfmt'] = islfmt

    solvent_file_format = dpars['solor']
    if solvent_file_format in ['boxes','cap']:
        isolor = 0
    elif solvent_file_format in ['in']:
        isolor = 1
    elif solvent_file_format in ['','none']:
        isolor = 2
    else:
        isolor = 99
        print('*** UNRECOGNIZED SOLVENT FORMAT IN THE PAR FILE ***')
        exit()
    dpars['isolor'] = isolor

    dpars['ivxyz'] = 1 if dpars['vxyz'] == 'vxyz' else 0
    if dpars['slab'] == 'slab':
        dpars['islab'] = 1
        dpars['ivxyz'] = 1
    else:
        dpars['islab'] = 0

    istart = 99
    dpars['ioptim'] = 0
    if dpars['icalc'] == 1:
        if islfmt == 0:
            if isolor == 1:
                istart = 3
            elif isolor in [0,2]:
                istart = 4
        elif islfmt == 1:
            if isolor == 0:
                istart = 5
            elif isolor == 1:
                istart = 9
        elif islfmt == 2:
            if isolor == 0:
                istart = 4
            elif isolor == 1:
                istart = 7
        elif islfmt == 3:
            if isolor == 0:
                istart = 8
        elif islfmt == 4:
            if isolor == 1:
                istart = 7
    elif dpars['icalc'] == 0:
        istart = 0
        if islfmt == 1:
            istart = 6
    elif dpars['icalc'] >= 2:
        dpars['islfmt'] = 0
        dpars['isolor'] = 2
        istart = 4
        opt = dpars['optim']
        if opt == 'simpl':
            ioptim = 0
        elif opt == 'powel':
            ioptim = 1
        elif opt == 'hooke':
            ioptim = 2
        elif opt == 'flepo':
            ioptim = 3
        elif opt == 'bfgs':
            ioptim = 4
        elif opt in ['conju','cg']:
            ioptim = 5
        elif opt in ['siman','sa']:
            ioptim = 8
        else:
            ioptim = 0
        dpars['ioptim'] = ioptim
        if dpars['ftol'] <= 0.0: dpars['ftol'] = 0.001
        if dpars['nsaevl'] <= 0: dpars['nsaevl'] = 1000000
        if dpars['nsatr'] <= 0: dpars['nsatr'] = 5
        if dpars['satemp'] <= 0: dpars['satemp'] = 5.0
        if dpars['sart'] <= 0 or dpars['sart'] >= 1.0: dpars['sart'] = 0.5
    dpars['istart'] = istart

    newzm = dpars['newzm']
    if newzm == 'full':
        inewzm = 2
    elif newzm == 'fullm':
        inewzm = 3
    elif newzm:
        inewzm = 1
    else:
        inewzm = 0
    dpars['inewzm'] = inewzm

    if dpars['ncents'] <= 0: dpars['ncents'] = 1
    dpars['icutas'] = icutas = [[1,0,0,0,0],[1,0,0,0,0]]
    dpars['ncutas'] = ncutas = [1,1]
    if dpars['icussl'] == 1:
        dpars['icutas'] = icutas = dpars['icutat']
        ncut = len(icutas) - icutas.count(0)
        if ncut == 0:
            icutas[dpars['modcus']][0] = dpars['ncents']
        else:
            ncutas[dpars['modcus']-1] = ncut        # minus 1 to start at 0

    if dpars['ncent1'] <= 0:
        dpars['ncent1'] = 1
    if dpars['ncent2'] <= 0:
        dpars['ncent2'] = dpars['ncent1']

    if dpars['nrota1'] <= 0:
        dpars['nrota1'] = 1
    if dpars['nrota2'] <= 0:
        dpars['nrota2'] = dpars['nrota1']

    if dpars['icapat'] != 0:
        if dpars['icapat'] <= 0:
            dpars['icapat'] = dpars['ncent1']
        if dpars['caprad'] <= 0.0:
            dpars['caprad'] = 15.0
        dpars['capsq'] = dpars['caprad'] ** 2
        dpars['capsq5'] = (dpars['caprad'] + 5.0) ** 2
        dpars['nvchg'] = 999999
        if dpars['icalc'] == 1 and dpars['isolor'] == 0:
            dpars['nmol'] = 9999

    rcut = dpars['rcut']
    if rcut <= 0.0:
        rcut = 8.5
        if dpars['modsv1'] > 2:
            rcut = 10.0
        if dpars['modsv1'] == 13:
            rcut = 8.5
        dpars['rcut'] = rcut

    dielrf = dpars['dielrf']
    if dielrf > 1.0:
        rcut = dpars['rcut']
        dpars['rffac'] = (dielrf-1.0) / ((2.0*dielrf + 1.0) * rcut**3)
        dpars['irfon'] = 1
    else:
        dpars['dielrf'] = 1.0
        dpars['irfon'] = 0
        dpars['rffac'] = 0.0

    scut = dpars['scut']
    if scut <= 0.0:
        scut = 8.5
        if dpars['modsv1'] > 2:
            scut = 10.0
        if dpars['modsv1'] == 13:
            scut = 8.5
        dpars['scut'] = scut

    cutnb = dpars['cutnb']
    if cutnb <= 0.0:
        dpars['cutnb'] = 16.0

    if dpars['wkc'] <= 0.0:
        dpars['wkc'] = 200.0
    if dpars['wkc'] > 1000.0:
        dpars['nopref'] = 1

    if dpars['torcut'] <= 0.0:
        dpars['torcut'] = 1.0

    if dpars['rdel'] <= 0.0:
        dpars['rdel'] = 0.15
        if dpars['modsv1'] > 2: dpars['rdel'] = 0.2

    if dpars['adel'] <= 0.0:
        dpars['adel'] = 15.0
        if dpars['modsv1'] > 2: dpars['adel'] = 20

    if dpars['scl14c'] <= 0.0:
        dpars['scl14c'] = 2.0
    elif dpars['scl14c'] > 99999.0:
        dpars['scl14c'] = 1.0 * 10**30
    dpars['scl41c'] = 1.0 / dpars['scl14c']

    if dpars['scl14l'] <= 0.0:
        dpars['scl14l'] = 2.0
    elif dpars['scl14l'] > 99999.0:
        dpars['scl14l'] = 1.0 * 10**30
    dpars['scl41l'] = 1.0 / dpars['scl14l']

    if dpars['rdlmin'] <= 0.0:
        dpars['rdlmin'] = 1.35

    if dpars['rdlinc'] <= 0.0:
        dpars['rdlinc'] = 0.1
    dpars['rinc'] = 1.0 / dpars['rdlinc']

    if dpars['eprmin'] <= 0.0:
        dpars['eprmin'] = -15.25

    if dpars['eprinc'] <= 0.0:
        dpars['eprinc'] = 0.5

    if dpars['edmin'] <= 0.0:
        dpars['edmin'] = -40.5

    if dpars['edinc'] <= 0.0:
        dpars['edinc'] = 1.0

    if dpars['essmin'] <= 0.0:
        dpars['essmin'] = -40.5

    if dpars['essinc'] <= 0.0:
        dpars['essinc'] = 1.0

    if dpars['ebsmin'] <= 0.0:
        dpars['ebsmin'] = -202.0

    if dpars['ebsinc'] <= 0.0:
        dpars['ebsinc'] = 4.0

    if dpars['nsym'] <= 0:
        dpars['nsym'] = 1

    if dpars['nconf'] <= 0:
        dpars['nconf'] = 1

    if dpars['fscale'] <= 0.0:
        dpars['fscale'] = 1.0

    if dpars['maxvar'] <= 0:
        dpars['maxvar'] = 15

    if dpars['pltfmt'] == 'pdb':
        iplfmt = 1
    elif dpars['pltfmt'] == 'pdbb':
        iplfmt = 2
    elif dpars['pltfmt'] == 'pdb2':
        iplfmt = 3
    elif dpars['pltfmt'] in ['mdlmo', 'mol', 'mdl', 'mdlsd']:
        iplfmt = 4
    else:
        iplfmt = 0
    dpars['iplfmt'] = iplfmt

    dpars['rcutsq'] = rcutsq = dpars['rcut'] ** 2
    dpars['scutsq'] = scutsq = dpars['scut'] ** 2
    dpars['ru2'] = ru2 = rcutsq
    dpars['su2'] = su2 = scutsq
    dpars['rl2'] = rl2 = (rcut-0.5) ** 2
    dpars['sl2'] = sl2 = (scut-0.5) ** 2
    dpars['rul'] = 1.0 / (ru2-rl2)
    dpars['sul'] = 1.0 / (su2-sl2)
    if dpars['nosmth'] == 1 or dpars['irfon'] == 1:
        dpars['rl2'] = 1.0 * 10**20
        dpars['sl2'] = 1.0 * 10**20

    dtorad = 0.017453292519943295
    dpars['adel'] *= dtorad
    dpars['rdl2'] = 2 * dpars['rdel']
    dpars['adl2'] = 2 * dpars['adel']
    dpars['adels'] = [0.14 for i in range(25)]
    dpars['adels'][0] = dpars['adels1'] * dtorad
    dpars['adels'][1] = dpars['adels2'] * dtorad
    dpars['rdels'] = [0.08 for i in range(25)]
    dpars['rdels'][0] = dpars['rdels1']
    dpars['rdels'][1] = dpars['rdels2']

    dpars['istart'] = istart
    assert istart == 4

    #TODO for future updates
    assert dpars['igbsa'] == 0
    assert dpars['iewald'] == 0
    assert dpars['icussl'] == 0
    assert dpars['islab'] == 0
    assert dpars['nofep'] == 0
    assert dpars['icapat'] == 0
    assert dpars['noss'] == 0

    return dpars


