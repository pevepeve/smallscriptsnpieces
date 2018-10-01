
########################
#        imports       #
########################

import os
import re
import shlex
import sys
import math
import numpy as np
import pandas as pd
import pickle

########################
#        functions     #
########################

def parse_xvg(fname, sel_columns='all'):
    """Parses XVG file legends and data"""
    
    _ignored = set(('legend', 'view'))
    _re_series = re.compile('s[0-9]+$')
    _re_xyaxis = re.compile('[xy]axis$')

    metadata = {}
    num_data = []
    
    metadata['labels'] = {}
    metadata['labels']['series'] = []

    ff_path = os.path.abspath(fname)
    if not os.path.isfile(ff_path):
        raise IOError('File not readable: {0}'.format(ff_path))
    
    with open(ff_path, 'r') as fhandle:
        for line in fhandle:
            line = line.strip()
            if line.startswith('@'):
                tokens = shlex.split(line[1:])
                if tokens[0] in _ignored:
                    continue
                elif tokens[0] == 'TYPE':
                    if tokens[1] != 'xy':
                        raise ValueError('Chart type unsupported: \'{0}\'. Must be \'xy\''.format(tokens[1]))
                elif _re_series.match(tokens[0]):
                    metadata['labels']['series'].append(tokens[-1])
                elif _re_xyaxis.match(tokens[0]):
                    metadata['labels'][tokens[0]] = tokens[-1]
                elif len(tokens) == 2:
                    metadata[tokens[0]] = tokens[1]
                else:
                    print('Unsupported entry: {0} - ignoring'.format(tokens[0]), file=sys.stderr)
            elif line[0].isdigit():
                num_data.append(map(float, line.split()))
    
    if not metadata['labels']['series']:
        for series in range(len(num_data) - 1):
            metadata['labels']['series'].append('')
    if sel_columns != 'all':
        sel_columns = map(int, sel_columns)
        x_axis = num_data[0]
        num_data = [x_axis] + [num_data[col] for col in sel_columns]
        metadata['labels']['series'] = [metadata['labels']['series'][col - 1] for col in sel_columns]
    
    return metadata, num_data




#############################################
#                 MAIN CODE                 #
#############################################


receptors = os.listdir('.')
allinone = []
allinone_4pickle = {}
for folder in receptors:
    if folder[0:2] == 're':
        code = folder[-4:]
        display(code)
########################
#     binding site     #
########################
        metadata_bs05, data_bs05 = parse_xvg(str('./'+folder+'/'+'rmsd_binding_site_05.xvg'), 'all')
        metadata_bs07, data_bs07 = parse_xvg(str('./'+folder+'/'+'rmsd_binding_site_07.xvg'), 'all')
        metadata_bs10, data_bs10 = parse_xvg(str('./'+folder+'/'+'rmsd_binding_site_10.xvg'), 'all')
        df_bs05 = pd.DataFrame(data_bs05, columns= ['time','rmsd'])
        df_bs07 = pd.DataFrame(data_bs07, columns= ['time','rmsd'])
        df_bs10 = pd.DataFrame(data_bs10, columns= ['time','rmsd'])
        df_bs_all = pd.merge(pd.merge(df_bs05, df_bs07, how='inner', on = 'time'), df_bs10, how='inner', on = 'time')
        df_bs_all.columns = ['time', 'rmsd_bs_05', 'rmsd_bs_07', 'rmsd_bs_10']
        
        
########################
#   receptor & ligand  #
########################
        metadata_rec, data_rec = parse_xvg(str('./'+folder+'/'+'rmsd_receptor.xvg'), 'all')
        metadata_lig, data_lig = parse_xvg(str('./'+folder+'/'+'rmsd_ligand.xvg'), 'all')
        df_rec = pd.DataFrame(data_rec, columns= ['time','rmsd'])
        df_lig = pd.DataFrame(data_lig, columns= ['time','rmsd'])
        df_reclig = pd.merge(df_rec, df_lig, how='inner', on = 'time')
        df_reclig.columns = ['time', 'rmsd_rec', 'rmsd_lig']
        df_reclig.time.apply(lambda x: x //1)
        df_reclig['time'] = df_reclig['time'] * 1000
        df_reclig.time = df_reclig.time.astype('int')

########################
#         SASA         #
########################
        metadata_sasa, data_sasa = parse_xvg(str('./'+folder+'/'+'sasa.xvg'), 'all')
        labels_sasa = ['time']
        for m in metadata_sasa['labels']['series']:
            labels_sasa.append(m)
        df_sasa = pd.DataFrame(data_sasa, columns=labels_sasa)

########################
#         Energies     #
########################
        metadata_ene, data_ene = parse_xvg(str('./'+folder+'/'+'energy_MM.xvg'), 'all')
        labels_ene = ['time']
        for m in metadata_ene['labels']['series']:
            labels_ene.append(m)
        df_ene = pd.DataFrame(data_ene, columns=labels_ene
########################
#         Final merge  #
########################
        df_all1 = pd.merge(pd.merge(pd.merge(df_bs_all, df_reclig, how='inner', on = 'time'),
                                   df_sasa, how='inner', on = 'time'),
                           df_ene, how='inner', on = 'time')


        df1 = df_all1.loc[:, (df_all1 != 0).any(axis=0)]
        df2 = df1.set_index('time')
        dfdict = df2.to_dict(orient='index')
        ndr = df1.values
        display(ndr.shape)
        allinone.append(ndr)
        allinone_4pickle[code] = dfdict
nall = np.array(allinone)
display(nall.shape)
outfile = 'numpy1.npy'
np.save(outfile, nall)
display(allinone_4pickle)
with open('data1.pkl', 'wb') as output:
    pickle.dump(allinone_4pickle, output, 2)
