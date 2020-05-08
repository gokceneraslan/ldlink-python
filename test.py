def test():
    api = LDlink(token='30e92210fa1b')
    
    res = api.ldproxy(var='rs3', pop='CEU')
    assert res is not None
    
    freq, haplo = api.ldhap(snps=['rs3', 'rs4'], pop='ALL')
    assert freq is not None
    assert haplo is not None
    
    resp = api.ldmatrix(snps=['rs3', 'rs4', 'rs148890987'], pop='CEU')
    assert resp is not None
    
    resp = api.ldmatrix(snps=['rs3', 'rs4', 'rs148890987'], pop='CEU', method='post')
    assert resp is not None
    
    resp = api.ldpair(var1='rs3', var2='rs4', pop='CEU')
    assert resp is not None
    
    resp = api.ldpop(var1='rs3', var2='rs4', pop='CEU')
    assert resp is not None
    
    resp = api.ldtrait(snps=['rs9834970', 'rs4'])
    assert resp is not None
    
    resp = api.snpclip(snps=['rs3', 'rs4'])
    assert resp is not None
    
    resp = api.snpchip(snps=['rs3', 'rs4'], platforms=['A_10X', 'A_250N', 'A_250S'])
    assert resp is not None