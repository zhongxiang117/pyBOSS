from PyBOSS.constants import SOLVENT_MODE, SOLVENT_PARS


def init_pars(g_x_inner_pars):
    """In-place editions of input parameter"""
    g_c_pconv = 0.000014576
    g_c_boltz = 0.00198717

    rc0 = g_x_inner_pars['rc0']
    rc1 = g_x_inner_pars['rc1']
    rc2 = g_x_inner_pars['rc2']
    if abs(rc0-rc1)<0.0001 and abs(rc1-rc2)<0.0001:
        nofep = 1
    else:
        nofep = 0
    if g_x_inner_pars['icalc'] > 1:
        nofep = 1
    g_x_inner_pars['nofep'] = nofep

    tk = g_x_inner_pars['t'] + 273.15
    beta = 1.0 / (g_c_boltz*tk)
    pvcon = g_x_inner_pars['p'] * g_c_pconv
    g_x_inner_pars['tk'] = tk
    g_x_inner_pars['beta'] = beta
    g_x_inner_pars['pvcon'] = pvcon

    g_x_inner_pars['betlht'] = 1.0 / (g_c_boltz*(g_x_inner_pars['tlht']+273.15))
    g_x_inner_pars['t1'] = g_x_inner_pars['t'] + 2.0
    g_x_inner_pars['t2'] = g_x_inner_pars['t'] - 2.0
    g_x_inner_pars['tk1'] = g_x_inner_pars['t'] + 2.0 + 273.15
    g_x_inner_pars['tk2'] = g_x_inner_pars['t'] - 2.0 + 273.15
    g_x_inner_pars['beta1'] = 1.0 / (g_c_boltz * g_x_inner_pars['tk1'])
    g_x_inner_pars['beta2'] = 1.0 / (g_c_boltz * g_x_inner_pars['tk2'])

    svmod1 = g_x_inner_pars['svmod1']
    if svmod1 and svmod1 in SOLVENT_PARS:
        sv1p = SOLVENT_PARS[svmod1]
    else:
        print(f'Fatal: currently not supported: svmod1: {svmod1}')
        exit()

    solvent1_modenum = sv1p[1]
    solvent1_atomnum = len(sv1p[3])
    g_x_inner_pars['solvent1_modenum'] = solvent1_modenum
    g_x_inner_pars['solvent1_atomnum'] = solvent1_atomnum
    g_x_inner_pars['igbsa'] = 1 if solvent1_modenum == 99 else 0
    g_x_inner_pars['modsv1'] = solvent1_modenum


    if solvent1_modenum > 4:
        x = -5.0 * g_x_inner_pars['nmol']
    else:
        x = -10.1 * g_x_inner_pars['nmol']
    g_x_inner_pars['ct1'] = x * (g_x_inner_pars['beta']-g_x_inner_pars['beta1'])
    g_x_inner_pars['ct2'] = x * (g_x_inner_pars['beta']-g_x_inner_pars['beta2'])


    if g_x_inner_pars['scllj'] == 0.0:
        g_x_inner_pars['scllj'] = 1.0


    svmod2 = g_x_inner_pars['svmod2']
    if svmod2 and svmod2 in SOLVENT_PARS:
        sv2p = SOLVENT_PARS[svmod2]
    else:
        print(f'Fatal: currently not supported: svmod2: {svmod2}')
        exit()

    solvent2_modenum = sv2p[1]
    solvent2_atomnum = len(sv2p[3])
    ncussl = 1 if solvent1_modenum == 20 else 0
    if solvent2_modenum == 0:
        pass
    elif solvent2_modenum+solvent1_modenum == 3:
        print('Do Not Use TIP3P and TIP4P Simultaneously - Aborted')
        exit()
    elif solvent2_modenum == 20:
        if solvent1_modenum == 20:
            print('Only One Solvent Can Be Flexible (Other) - Aborted')
            exit()
        else:
            ncussl = 2
    g_x_inner_pars['icussl'] = 1 if ncussl != 0 else 0
    g_x_inner_pars['i_custom_solvent'] = ncussl
    g_x_inner_pars['solvent2_modenum'] = solvent2_modenum
    g_x_inner_pars['modsv2'] = solvent2_modenum
    g_x_inner_pars['solvent2_atomnum'] = solvent2_atomnum
    g_x_inner_pars['ncutas'] = [1,1]
    if g_x_inner_pars['ncents'] <= 0:
        g_x_inner_pars['ncents'] = 1
    icuta = g_x_inner_pars['icutas']
    g_x_inner_pars['icuta'] = [i for i in icuta]
    icutas = [[0 for i in range(5)] for j in range(2)]  # index starts from zero
    icutas[0][0] = 1
    icutas[1][0] = 1
    ncutas = [1,1]
    if g_x_inner_pars['icussl'] == 1:
        ncut = len(icuta) - icuta.count(0)
        if ncut == 0:
            icutas[ncussl][0] = g_x_inner_pars['ncents']
        else:
            ncutas[ncussl-1] = ncut
    g_x_inner_pars['icutas'] = icutas
    g_x_inner_pars['ncutas'] = ncutas


    if g_x_inner_pars['rsolv'] <= 0.0:
        radius = 1.4
        if solvent1_modenum < 99:
            radius = 3.2
        if solvent1_modenum < 3 or solvent1_modenum == 13:
            radius = 1.4
        g_x_inner_pars['rsolv'] = radius


    if g_x_inner_pars['noss'] != 0:
        g_x_inner_pars['nrdf'] = 0
        g_x_inner_pars['nrdfs'] = 0
        g_x_inner_pars['nofep'] = 1
        g_x_inner_pars['nopref'] = 1
        g_x_inner_pars['wkc'] = 20000.0

    if g_x_inner_pars['wkc'] <= 0.0:
        g_x_inner_pars['wkc'] = 200.0
    if g_x_inner_pars['wkc'] > 1000.0:
        g_x_inner_pars['nopref'] = 1


    qmname = g_x_inner_pars['qmname']
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
        print('Unavailable QM Choice: QMNAME: - Aborted')
        exit(-1)
    g_x_inner_pars['iampm'] = iampm

    if iampm != 0:
        g_x_inner_pars['isqm'] = 1
    else:
        g_x_inner_pars['isqm'] = 0

    cmname = g_x_inner_pars['cmname']
    if cmname in ['mull']:
        ichargemode = 2
    elif cmname in ['cm1','cm1a','cm1p']:
        ichargemode = 3
    elif cmname in ['cm3','cm3a','cm3p']:
        ichargemode = 4
    elif cmname in ['gpar']:
        ichargemode = 5
        iampm = 5
    elif cmname in ['','opls','none']:
        ichargemode = 0
    else:
        ichargemode = -1
        print('Unavailable QM Choice: CMNAME: - Aborted')
        exit(-1)
    g_x_inner_pars['ichargemode'] = ichargemode


    qmscale = g_x_inner_pars['qmscale']
    if qmscale <= 0.:
        qmscale = 1.2
        if g_x_inner_pars['rc0'] != 0.0:
            qmscale = 1.0
    g_x_inner_pars['qmscale'] = qmscale


    isolec = g_x_inner_pars['isolec']
    if isolec <= 0:
        isolec = 1
    g_x_inner_pars['isolec'] = isolec


    slfmt = g_x_inner_pars['slfmt']
    if slfmt in ['zmat']:
        solute_filetype = 0
    elif slfmt in ['pdb','pdbp']:
        #TODO
        solute_filetype = 1
    elif slfmt in ['mind']:
        #TODO
        solute_filetype = 2
    elif slfmt in ['in']:
        solute_filetype = 3
    elif slfmt in ['zin']:
        solute_filetype = 4
    else:
        solute_filetype = 99
        print('*** UNRECOGNIZED SOLUTE FORMAT IN THE PAR FILE ***')
        exit()
    g_x_inner_pars['solute_filetype'] = solute_filetype


    solvent_file_format = g_x_inner_pars['solor']
    if solvent_file_format in ['boxes','cap']:
        solvent_filetype = 0
    elif solvent_file_format in ['in']:
        solvent_filetype = 1
    elif solvent_file_format in ['','none']:
        solvent_filetype = 2
    else:
        solvent_filetype = 99
        print('*** UNRECOGNIZED SOLVENT FORMAT IN THE PAR FILE ***')
        exit()
    g_x_inner_pars['solvent_filetype'] = solvent_filetype


    istart = 99
    ioptimize = 0
    if g_x_inner_pars['icalc'] == 1:
        if solute_filetype == 0:
            if solvent_filetype == 1:
                istart = 3
            elif solvent_filetype in [0,2]:
                istart = 4
        elif solute_filetype == 1:
            if solvent_filetype == 0:
                istart = 5
            elif solvent_filetype == 1:
                istart = 9
        elif solute_filetype == 2:
            if solvent_filetype == 0:
                istart = 4
            elif solvent_filetype == 1:
                istart = 7
        elif solute_filetype == 3:
            if solvent_filetype == 0:
                istart = 8
        elif solute_filetype == 4:
            if solvent_filetype == 1:
                istart = 7
    elif g_x_inner_pars['icalc'] == 0:
        istart = 0
        if solute_filetype == 1:
            istart = 6
    elif g_x_inner_pars['icalc'] >= 2:
        solute_filetype = 0
        solvent_filetype = 2
        istart = 4
        opt = g_x_inner_pars['optim']
        if opt == 'simpl':
            ioptimize = 0
        elif opt == 'powel':
            ioptimize = 1
        elif opt == 'hooke':
            ioptimize = 2
        elif opt == 'flepo':
            ioptimize = 3
        elif opt == 'bfgs':
            ioptimize = 4
        elif opt in ['conju','cg']:
            ioptimize = 5
        elif opt in ['siman','sa']:
            ioptimize = 8
        else:
            ioptimize = 0
        g_x_inner_pars['optim'] = ioptimize
        #TODO simulation annealing
    g_x_inner_pars['istart'] = istart

    solvent_solvent_cutoff = g_x_inner_pars['rcut']
    if solvent_solvent_cutoff <= 0.0:
        solvent_solvent_cutoff = 8.5
        if solvent1_modenum > 2:
            solvent_solvent_cutoff = 10.0
        if solvent1_modenum == 13:      # specific
            solvent_solvent_cutoff = 8.5
        g_x_inner_pars['rcut'] = solvent_solvent_cutoff

    solute_solvent_cutoff = g_x_inner_pars['scut']
    if solute_solvent_cutoff <= 0.0:
        solute_solvent_cutoff = 8.5
        if solvent1_modenum > 2:
            solute_solvent_cutoff = 10.0
        if solvent1_modenum == 13:      # specific
            solute_solvent_cutoff = 8.5
        g_x_inner_pars['scut'] = solute_solvent_cutoff

    fdielec = g_x_inner_pars['dielec']
    if g_x_inner_pars['noss'] == 0:
        g_x_inner_pars['scllj'] = 1.0
        if solvent1_modenum != 99:
            fdielec = 1.0
    if fdielec <= 0.0:
        fdielec = 1.0
    g_x_inner_pars['dielec'] = fdielec

    if fdielec > 1.0:
        g_x_inner_pars['dielrf_factor'] = (fdielec-1.0) / ((2.*fdielec+1.)*solvent_solvent_cutoff**3)
        g_x_inner_pars['irfon'] = 1
    else:
        g_x_inner_pars['dielrf'] = 1.0
        g_x_inner_pars['dielrf_factor'] = 0.0
        g_x_inner_pars['irfon'] = 0

    nonbonded_cutoff = g_x_inner_pars['cutnb']
    if nonbonded_cutoff <= 0.0:
        g_x_inner_pars['cutnb'] = 16

    if g_x_inner_pars['torcut'] <= 0.0:
        g_x_inner_pars['torcut'] = 1.0

    g_x_inner_pars['ivxyz'] = 1 if g_x_inner_pars['vxyz'] == 'vxyz' else 0
    if g_x_inner_pars['slab'] == 'slab':
        g_x_inner_pars['islab'] = 1
        g_x_inner_pars['ivxyz'] = 1
    else:
        g_x_inner_pars['islab'] = 0
    
    if g_x_inner_pars['nvchg'] <= 0:
        g_x_inner_pars['nvchg'] = 125 * (g_x_inner_pars['nmol']//20)
    if g_x_inner_pars['nvchg'] <= 0:
        g_x_inner_pars['nvchg'] = 999999
    
    if g_x_inner_pars['nschg'] <= 0:
        g_x_inner_pars['nschg'] = 10 * (g_x_inner_pars['nmol']//40)
    if g_x_inner_pars['nschg'] <= 0:
        g_x_inner_pars['nschg'] = g_x_inner_pars['nmol'] + 1

    if g_x_inner_pars['nconsv']  == 0:
        g_x_inner_pars['nconsv']  = 999999

    if g_x_inner_pars['vdel'] <= 0.0:
        g_x_inner_pars['vdel'] = float(g_x_inner_pars['nmol']//20*10)

    if g_x_inner_pars['isolec'] <= 0:
        g_x_inner_pars['isolec'] = 1
    
    if g_x_inner_pars['ncent2'] <= 0:
        g_x_inner_pars['ncent2'] = g_x_inner_pars['ncent1']

    icutat = g_x_inner_pars['icutat']
    if g_x_inner_pars['icut'] == 3:
        g_x_inner_pars['ncutat'] = 2
        g_x_inner_pars['icutat'][0] = g_x_inner_pars['ncent1']
        g_x_inner_pars['icutat'][1] = g_x_inner_pars['ncent2']
    elif g_x_inner_pars['icut'] == 4:
        g_x_inner_pars['ncutat'] = len(icutat) - icutat.count(0)
    else:
        g_x_inner_pars['ncutat'] = 0
    
    # make sure atom pairs of RDFs are set correctly
    n = g_x_inner_pars['nrdf']
    if n == 0:
        l1 = list(set(g_x_inner_pars['nrdfs1']))
        assert len(l1) == 1 and l1[0] == 0
        l2 = list(set(g_x_inner_pars['nrdfs2']))
        assert len(l2) == 1 and l2[0] == 0
    else:
        l1 = g_x_inner_pars['nrdfs1']
        m = max([l1.count(i) for i in set(l1) if i != 0])   # ignore `zero`s 
        assert m <= n and 0 not in l1[:n]
        l2 = g_x_inner_pars['nrdfs2']
        m = max([l2.count(i) for i in set(l2) if i != 0])
        assert m <= n and 0 not in l2[:n]
    g_x_inner_pars['rdfs_solvent_solvent'] = [[l1[i],l2[i]] for i in range(n)]

    n = g_x_inner_pars['nrdfs']
    if n == 0:
        l1 = list(set(g_x_inner_pars['nrdfa1']))
        assert len(l1) == 1 and l1[0] == 0
        l2 = list(set(g_x_inner_pars['nrdfa2']))
        assert len(l2) == 1 and l2[0] == 0
    else:
        l1 = g_x_inner_pars['nrdfa1']
        m = max([l1.count(i) for i in set(l1) if i != 0])
        assert m <= n and 0 not in l1[:n]
        l2 = g_x_inner_pars['nrdfa2']
        m = max([l2.count(i) for i in set(l2) if i != 0])
        assert m <= n and 0 not in l2[:n]
    g_x_inner_pars['rdfs_solute_solvent'] = [f'{l1[i]}-{l2[i]}' for i in range(n)]

    if g_x_inner_pars['ebsmin'] <= 0.0:
        g_x_inner_pars['ebsmin'] = -202.0


    #TODO for future updates
    assert g_x_inner_pars['igbsa'] == 0
    assert g_x_inner_pars['iewald'] == 0
    assert g_x_inner_pars['icussl'] == 0
    assert g_x_inner_pars['islab'] == 0
    assert g_x_inner_pars['nofep'] == 0
    assert g_x_inner_pars['icapat'] == 0

    return g_x_inner_pars


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

    g_x_solvents_data = {
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
    g_x_solvents_data.update(spars)
    g_x_solvents_data['natmx'] = max(len(solvent1_pars[2]),len(solvent2_pars[2]))

    sqrtesq = pow(g_c_esq,0.5)
    for s in ['1','2']:
        p = g_x_solvents_data[s]
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
    n = len(g_x_solvents_data['1']['QW'])
    for i in range(n):
        q = g_x_solvents_data['1']['QW'][i]
        a = g_x_solvents_data['1']['AW'][i]
        b = g_x_solvents_data['1']['BW'][i]
        for j in range(n):
            qq = g_x_solvents_data['1']['QW'][j] * q
            aa = g_x_solvents_data['1']['AW'][j] * a
            bb = g_x_solvents_data['1']['BW'][j] * b
            key = f'{i}-{j}'
            inters[key] = [qq,aa,bb]
    g_x_solvents_data['1-1'] = inters

    return g_x_solvents_data






