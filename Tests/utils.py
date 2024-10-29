

def allclose(a,b,tol=None):
    """emulate numpy.allclose: recursively compare a and b on given tolerance"""
    tol = tol if tol else 0.000001
    if a is None and b is None:
        return True
    elif a is True and b is True:
        return True
    elif a is False and b is False:
        return True
    elif isinstance(a, (int,float)) and isinstance(b, (int,float)):
        if abs(a-b) > tol:
            return False
        return True
    elif isinstance(a, str) and isinstance(b, str):
        if a == b:
            return True
        return False
    elif isinstance(a, (list,tuple)) and isinstance(b, (list,tuple)):
        if len(a) != len(b):
            return False
        if not a:
            return True
        for i,vi in enumerate(a):
            vj = b[i]
            if not allclose(vi,vj,tol):
                return False
        return True
    elif isinstance(a, dict) and isinstance(b, dict):
        akeys = sorted(a.keys())
        bkeys = sorted(b.keys())
        if len(akeys) != len(bkeys):
            return False
        if not allclose(akeys,bkeys,tol):
            return False
        for k in akeys:
            if not allclose(a[k],b[k],tol):
                return False
        return True
    else:
        return False



