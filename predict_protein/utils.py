# -*- coding: utf-8 -*-

""" utility functions for notebooks and various scripts """

import pandas as pd
import pickle

def get_dataframe(path):
    with open(path, 'rb') as f:
        results = pickle.load(f)
        df_ = pd.concat([result['metrics'] for result in results])
    return df_