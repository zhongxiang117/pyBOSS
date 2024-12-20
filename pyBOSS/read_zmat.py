import os


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


