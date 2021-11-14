import cptac
import numpy as np

#%%
import sys
import os
import argparse
import pickle

import pandas as pd
from predict_protein import select_features, train_model, download_cptac, __version__


#%%
def predict_protein(args):
    #%%
    cptac_list = download_cptac.download_cptac(n_tumors=args.n_tumors)

    #%%
    # Example combining 2 tumors then learn against self using an elastic net
    # TODO: Can we speed this up?

    tumor_df = pd.concat(cptac_list)
    tm = train_model.LearnCPTAC(tumor_df)
    #self_elastic_result = comb_2tumors.learn_all_proteins(tx_to_include="string",
    #                                                  train_method="elastic")

    #%%
    tm.included_features=args.features
    tm.train_method=args.method
    pred = tm.learn_all_proteins()
    pickle.dump(pred, open(args.out, 'wb'))

    return sys.exit(os.EX_OK)

def main():
    """
    parsing command line arguments
    :return:
    """

    # Main command
    parser = argparse.ArgumentParser(description='Predicts proteins from transcripts')

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))

    parser.add_argument('-n', '--n_tumors', type=int,
                        choices=range(0, 10),
                        help='number of types of tumors',
                        )

    parser.add_argument('-o', '--out', type=str,
                        help='path to output file',
                        default='out.p')

    parser.add_argument('-f', '--features',
                        type=str,
                        choices=['single', 'all', 'string', 'corum'],
                        help='which transcripts to use for prediction',
                        default='single',
                        )

    parser.add_argument('-m', '--method',
                        type=str,
                        choices=['elastic', 'voting', 'forest', 'linreg'],
                        help='model to use for fitting',
                        default='linreg',)

    parser.set_defaults(func=predict_protein)

    # Print help message if no arguments are given
    import sys
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # gc.enable()
    # gc.set_debug(gc.DEBUG_LEAK)

    # Parse all the arguments
    args = parser.parse_args()

    # Run the function in the argument
    args.func(args)
