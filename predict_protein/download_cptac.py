# -*- coding: utf-8 -*-

""" download cptac data"""

import cptac
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler

def download_cptac(n_tumors: int = 2,
                   ):
    # List current CPTAC datasets
    cptac.list_datasets()

    cptac.download(dataset="endometrial")
    cptac.download(dataset="ovarian")
    cptac.download(dataset="colon")
    cptac.download(dataset="brca")
    cptac.download(dataset="luad")
    cptac.download(dataset="ccrcc")
    cptac.download(dataset="gbm")
    cptac.download(dataset="lscc")
    cptac.download(dataset="hnscc")

    # %%
    en = cptac.Endometrial()
    ov = cptac.Ovarian()
    co = cptac.Colon()
    br = cptac.Brca()
    lu = cptac.Luad()
    cc = cptac.Ccrcc()
    gb = cptac.Gbm()
    ls = cptac.Lscc()
    hn = cptac.Hnscc()

    # For endometrial, try getting the RNA and protein data
    en_rna = en.get_transcriptomics()
    en_pro = en.get_proteomics()
    a = en.join_omics_to_omics('transcriptomics', 'proteomics')

    ov_rna = ov.get_transcriptomics()
    ov_pro = ov.get_proteomics()
    b = ov.join_omics_to_omics('transcriptomics', 'proteomics')
    b.columns = b.columns.droplevel(1)

    co_rna = co.get_transcriptomics()
    co_pro = co.get_proteomics()
    c = co.join_omics_to_omics('transcriptomics', 'proteomics')

    br_rna = br.get_transcriptomics()
    br_pro = br.get_proteomics()
    d = br.join_omics_to_omics('transcriptomics', 'proteomics')
    d.columns = d.columns.droplevel(1)

    lu_rna = lu.get_transcriptomics()
    lu_pro = lu.get_proteomics()
    e = lu.join_omics_to_omics('transcriptomics', 'proteomics')
    e.columns = e.columns.droplevel(1)

    cc_rna = cc.get_transcriptomics()
    cc_pro = cc.get_proteomics()
    f = cc.join_omics_to_omics('transcriptomics', 'proteomics')
    f.columns = f.columns.droplevel(1)

    gb_rna = gb.get_transcriptomics()
    gb_pro = gb.get_proteomics()
    g = gb.join_omics_to_omics('transcriptomics', 'proteomics')
    g.columns = g.columns.droplevel(1)

    ls_rna = ls.get_transcriptomics()
    ls_pro = ls.get_proteomics()
    h = ls.join_omics_to_omics('transcriptomics', 'proteomics')
    h.columns = h.columns.droplevel(1)

    hn_rna = hn.get_transcriptomics()
    hn_pro = hn.get_proteomics()
    i = hn.join_omics_to_omics('transcriptomics', 'proteomics')

    # Transform
    a_std = a.copy()
    a_tx_cols = [col for col in a_std.columns if col.endswith('transcriptomics')]
    a_std[a_tx_cols] = RobustScaler().fit_transform(a_std[a_tx_cols])
    a_std.index = 'EN' + a_std.index

    # Log transform
    b_std = b.copy()
    b_std = b_std.loc[:, ~b_std.columns.duplicated(keep='first')]
    b_tx_cols = [col for col in b_std.columns if col.endswith('transcriptomics')]
    b_std[b_tx_cols] = np.log2(b_std[b_tx_cols])
    b_std[np.isinf(b_std)] = np.nan
    b_std[b_tx_cols] = RobustScaler().fit_transform(b_std[b_tx_cols])
    b_std.index = 'OV' + b_std.index

    c_std = c.copy()
    c_tx_cols = [col for col in c_std.columns if col.endswith('transcriptomics')]
    c_std[c_tx_cols] = RobustScaler().fit_transform(c_std[c_tx_cols])
    c_std.index = 'CO' + c_std.index

    d_std = d.copy()
    d_std = d_std.loc[:, ~d_std.columns.duplicated(keep='first')]
    d_tx_cols = [col for col in d_std.columns if col.endswith('transcriptomics')]
    d_std[d_tx_cols] = RobustScaler().fit_transform(d_std[d_tx_cols])
    d_std.index = 'BR' + d_std.index

    e_std = e.copy()
    e_std = e_std.loc[:, ~e_std.columns.duplicated(keep='first')]
    e_tx_cols = [col for col in e_std.columns if col.endswith('transcriptomics')]
    e_std[e_tx_cols] = RobustScaler().fit_transform(e_std[e_tx_cols])
    e_std.index = 'LU' + e_std.index

    # Log transform
    f_std = f.copy()
    f_std = f_std.loc[:, ~f_std.columns.duplicated(keep='first')]
    f_tx_cols = [col for col in f_std.columns if col.endswith('transcriptomics')]
    f_std[f_tx_cols] = np.log2(f_std[f_tx_cols])
    f_std[np.isinf(f_std)] = np.nan
    f_std[f_tx_cols] = RobustScaler().fit_transform(f_std[f_tx_cols])
    f_std.index = 'CC' + f_std.index

    # Log transform
    g_std = g.copy()
    g_std = g_std.loc[:, ~g_std.columns.duplicated(keep='first')]
    g_tx_cols = [col for col in g_std.columns if col.endswith('transcriptomics')]
    g_std[g_tx_cols] = np.log2(g_std[g_tx_cols])
    g_std[np.isinf(g_std)] = np.nan
    g_std[g_tx_cols] = RobustScaler().fit_transform(g_std[g_tx_cols])
    g_std.index = 'GB' + g_std.index

    h_std = h.copy()
    h_std = h_std.loc[:, ~h_std.columns.duplicated(keep='first')]
    h_tx_cols = [col for col in h_std.columns if col.endswith('transcriptomics')]
    h_std[h_tx_cols] = RobustScaler().fit_transform(h_std[h_tx_cols])
    h_std.index = 'LS' + h_std.index

    i_std = i.copy()
    i_tx_cols = [col for col in i_std.columns if col.endswith('transcriptomics')]
    i_std[i_tx_cols] = RobustScaler().fit_transform(i_std[i_tx_cols])
    i_std.index = 'HN' + i_std.index

    return [h_std, f_std, e_std, g_std, a_std, d_std, b_std, c_std, ][:n_tumors]