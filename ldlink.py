import requests
import pandas as pd
import io, json


class LDlink:
    def __init__(self, token):
        self.token = token
        self.url = 'https://ldlink.nci.nih.gov/LDlinkRest/'
        
    def _isjson(self, resp):
        return 'json' in resp.headers.get('content-type')
    
    def _check_r2_d(self, r2):
        return r2 in ('r2', 'd')
        
    def ldproxy(self, var, pop='CEU', r2_d='r2'):
        assert self._check_r2_d(r2_d)
        
        url = self.url + 'ldproxy'
        resp = requests.get(url, params=dict(
            token=self.token,
            var=var,
            pop=pop,
            r2_d=r2_d,
        ))

        if resp.ok and not self._isjson(resp):
            return pd.read_table(io.StringIO(resp.text))
        
    def ldhap(self, snps, pop='CEU'):
        assert isinstance(snps, (list, tuple))
        assert len(snps) <= 30
        
        snps = '\n'.join(snps)
        url = self.url + 'ldhap'
        resp = requests.get(url, params=dict(
            token=self.token,
            snps=snps,
            pop=pop,
        ))
        split = '#####################################################################################'
        
        if resp.ok and not self._isjson(resp) and split in resp.text:
            text = resp.text
            allele_freq = text.split(split)[0]
            haplo = text.split(split)[1]
        
            return pd.read_table(io.StringIO(allele_freq)), pd.read_table(io.StringIO(haplo))
        
    def ldmatrix(self, snps, pop='CEU', r2_d='r2', method='auto'):
        assert self._check_r2_d(r2_d)
        assert isinstance(snps, (list, tuple))
        assert len(snps) <= 1000

        url = self.url + 'ldmatrix'
        snps = '\n'.join(snps)
        
        if method == 'auto':
            if len(snps) <= 300:
                method = 'get'
            else:
                method = 'post'
            
        if method == 'get':
            assert len(snps) <= 300
            resp = requests.get(url, params=dict(
                token=self.token,
                snps=snps,
                pop=pop,
                r2_d=r2_d,
            ))
        elif method == 'post':
            url = url + '?token=' + self.token
            resp = requests.post(url, json=dict(
                snps=snps,
                pop=pop,
                r2_d=r2_d,
            ))
        else:
            return

        if resp.ok and not self._isjson(resp):
            return pd.read_table(io.StringIO(resp.text))
        
    def ldpair(self, var1, var2, pop='CEU'):
        url = self.url + 'ldpair'
        resp = requests.get(url, params=dict(
            token=self.token,
            var1=var1,
            var2=var2,
            pop=pop,
        ))
        if resp.ok and not self._isjson(resp):
            return resp.text
        
    def ldpop(self, var1, var2, pop='CEU', r2_d='r2'):
        assert self._check_r2_d(r2_d)

        url = self.url + 'ldpop'
        resp = requests.get(url, params=dict(
            token=self.token,
            var1=var1,
            var2=var2,
            pop=pop,
            r2_d=r2_d,            
        ))
        if resp.ok and not self._isjson(resp):
            return pd.read_table(io.StringIO(resp.text))
        
    def ldtrait(self, snps, pop='CEU', r2_d='r2', r2_d_threshold=0.1, window=500000):
        assert isinstance(snps, (list, tuple))
        assert self._check_r2_d(r2_d)
        assert len(snps) <= 300

        snps = '\n'.join(snps)
        url = self.url + 'ldtrait?token=' + self.token
        jsn = dict(
            snps=snps,
            pop=pop,
            r2_d=r2_d,
            r2_d_threshold=str(r2_d_threshold),
            window=str(window),
        )
        resp = requests.post(url, json=jsn)

        if resp.ok and not self._isjson(resp):
            return pd.read_table(io.StringIO(resp.text))
        
        
    def snpclip(self, snps, pop='CEU', r2_threshold=0.1, maf_threshold=0.01):
        assert isinstance(snps, (list, tuple))
        
        snps = '\n'.join(snps)
        url = self.url + 'snpclip?token=' + self.token
        jsn = dict(
            snps=snps,
            pop=pop,
            r2_threshold=str(r2_threshold),
            maf_threshold=str(maf_threshold),            
        )
        resp = requests.post(url, json=jsn)

        if resp.ok and not self._isjson(resp):
            return pd.read_table(io.StringIO(resp.text))
        
    def snpchip(self, snps, platforms):
        assert isinstance(snps, (list, tuple))
        assert isinstance(platforms, (list, tuple))

        snps = '\n'.join(snps)
        platforms = '+'.join(platforms)
        url = self.url + 'snpchip?token=' + self.token
        jsn = dict(
            snps=snps,
            platforms=platforms,
        )

        resp = requests.post(url, json=jsn)
        if resp.ok and not self._isjson(resp):
            return pd.read_table(io.StringIO(resp.text))