# -*- coding: utf-8 -*-

""" downloads cptac data and predict with scikit-learn """

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import re
import tqdm
import numpy as np
import os
import datetime

from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor, VotingRegressor, GradientBoostingRegressor
from sklearn.linear_model import LinearRegression, ElasticNetCV

from . import select_features, params
from .select_features import GetProtein, GetComplex


class LearnCPTAC(object):
    """
    df
    """

    def __init__(self, cptac_df):
        """

        :param cptac_df:
        """
        self.df = cptac_df

        self.all_proteomics = [re.sub('_proteomics', "", protein) for protein in self.df.columns if
                               protein.endswith('_proteomics')]
        self.all_transcriptomics = [re.sub('_transcriptomics', "", transcript) for transcript in self.df.columns if
                                    transcript.endswith('_transcriptomics')]

        # these are the proteins with both transcriptomics and proteomics data in the cptac-df
        self.shared_proteins = [protein for protein in self.all_proteomics if protein in self.all_transcriptomics]

        # Remove zero variance threshold
        self.var_threshold = 0.2

        self._features = 'single'
        self._method = 'linreg'

        self.stringdb = GetProtein()  # Disabled for now, need to figure out import path to get_proteins
        self.corumdb = GetComplex()

        pass

    @property
    def train_method(self):
        return self._method

    @train_method.setter
    def train_method(self, value):
        """
        Set train method
        :param value:   linreg, elastic, forest, boosting, or voting
        :return:
        """

        self._method = value

    @property
    def included_features(self):
        return self._features

    @included_features.setter
    def included_features(self, value):
        """
        Set features
        :param value:     single, all, string, stringhi, or corum
        :return:
        """
        self._features = value

    def learn_all_proteins(
            self,
            n_threads=1,
            # tx_to_include=self._features,
            # train_method=self._method,
    ):
        """
        Wrapper for learning one protein

        :param n_threads:     int: number of threads
        :return:
        """

       # self.included_features = tx_to_include
       # self.train_method = train_method



        # TODO: 2021-11-17 trying to parallelize this. Not sure if it is ok to simply let the function read the
        # class or whether we need to wrap it a partial function

        # learning_results = []
        # for i, protein in enumerate(tqdm.tqdm(self.shared_proteins,
        #                                       desc=f'Running {self._method} model with {self._features} transcripts')):
        #     tqdm.tqdm.write(f'Doing {i}: {protein}')
        #
        #     res = self.learn_one_protein(protein, n_threads=n_threads)
        #
        #     if res:
        #         learning_results.append(res)
        #         corr = res['metrics'].corr_test.values
        #         r2 = res['metrics'].r2_test.values
        #         tqdm.tqdm.write(f'corr: {corr}, r2: {r2}')

        loop_ = self.shared_proteins

        from concurrent import futures
        with futures.ThreadPoolExecutor(max_workers=n_threads) as ex:
            res = list(tqdm.tqdm(ex.map(self.learn_one_protein, loop_),
                       total=len(loop_),
                       desc=f'Running {self._method} model with {self._features} transcripts'))

        learning_results = [r for r in res if r is not None]
        return learning_results

    def get_train_test(self,
                       protein_to_do,
                       ):
        """
        get the train-test data frames for a particular protein
        :param protein_to_do:       str: protein
        :return:  x_train, x_test, y_train, y_test
        """

        y_df = self.df[[protein_to_do + '_proteomics']]
        y_df = y_df.dropna(subset=[protein_to_do + '_proteomics'])

        # skip proteins with fewer than 20 samples
        # 2021-11-12 this should be filtered at the protein step (y_df) rather than the
        # transcript-protein joined set (xy_df). We will do a simple impute for the transcripts instead
        if len(y_df) < params.min_proteins:
            tqdm.tqdm.write('Not enough proteomics observations. Skipping protein.')
            return None

        # Decide how many proteins to use:
        if self._features == "single":
            proteins_to_include = [protein_to_do]

        elif self._features == "all":
            proteins_to_include = self.shared_proteins

        elif self._features == "string":
            string_interactors = self.stringdb.find_interactor(protein_to_do,
                                                               combined_score=200,
                                                               max_proteins=1000,
                                                               )
            proteins_to_include = [p for p in string_interactors if p in self.shared_proteins]

        elif self._features == "stringhi":
            string_interactors = self.stringdb.find_interactor(protein_to_do,
                                                               combined_score=800,
                                                               max_proteins=1000,
                                                               )
            proteins_to_include = [p for p in string_interactors if p in self.shared_proteins]

        elif self._features == "corum":
            corum_interactors = self.corumdb.find_cosubunits(protein_to_do)
            proteins_to_include = [p for p in corum_interactors if p in self.shared_proteins]

        else:
            raise Exception('Invalid features')

        # TODO: 2021-11-13 I found the problem of why the string method is returning fewer model,
        #  it is because we are doing an inner join of all the transcript features with the protein
        #  and taking only complete rows. Because some features could have very few observations,
        #  this causes the entire data frame to not have enough
        #  observations for training.
        #  Currently I have implemented a temporary solution to impute and remove near zero var columns

        # create transcript dataframe with the appropriate features
        x_df = self.df[[tx + '_transcriptomics' for tx in proteins_to_include]]

        # simple impute for missing values
        # TODO: use knn for better performance.
        x_impute = SimpleImputer(missing_values=np.nan, strategy='median').fit_transform(x_df)
        x_impute_df = pd.DataFrame(x_impute, columns=x_df.columns, index=x_df.index)

        # remove low variance columns
        try:
            vt = VarianceThreshold().fit(x_impute_df)

        # catch instances where no column has any variance
        except ValueError:
            tqdm.tqdm.write('No feature with non-zero variance. Skipping protein.')
            return None

        x_impute_var_df = x_impute_df.iloc[:, np.where(vt.variances_ > self.var_threshold)[0]]

        # skip proteins where there are insufficient features
        if x_impute_var_df.shape[1] == 0:
            tqdm.tqdm.write('Not enough features. Skipping protein.')
            return None

        # Join x and y data frames
        xy_df = x_impute_var_df.join(y_df, how='inner').copy().dropna()

        # perform train-test split
        x = xy_df.iloc[:, :-1]  # .values
        y = xy_df.iloc[:, -1]  # .values

        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=2)

        # perform pca on the train and test
        # from sklearn.decomposition import PCA
        # pca = PCA(n_components=15)
        # x_train_pca = pca.fit_transform(x_train)
        # x_test_pca = pca.transform(x_test)

        return x_train, x_test, y_train, y_test

    def learn_one_protein(
            self,
            protein_to_do,
            n_threads: int = 1,
    ):
        """
        train model for one protein of interest

        :param protein_to_do:
        :param n_threads:       int: threads
        :return: dict of models and metrics
        """

        try:
            x_train, x_test, y_train, y_test = self.get_train_test(protein_to_do=protein_to_do)

        # if the get_train_test method returns none, return empty result
        except TypeError:
            return None

        if self._method == 'linreg':
            vreg = LinearRegression(n_jobs=-1)

        elif self._method == 'voting':
            # Train model
            elastic = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9, 0.95],
                                   cv=5,
                                   fit_intercept=False,
                                   n_jobs=4,
                                   tol=1e-3,
                                   max_iter=2000,
                                   random_state=2,
                                   )

            # Shallow Forest was 100 estimators 5 depth. Deep is 1000 estimators 10 depth.
            forest = RandomForestRegressor(n_estimators=500,
                                           criterion='squared_error',
                                           max_depth=4,
                                           random_state=2,
                                           # oob_score=True,
                                           n_jobs=-1)

            vreg = VotingRegressor(estimators=[('en', elastic), ('rf', forest), ])

        elif self._method == 'forest':
            vreg = RandomForestRegressor(n_estimators=500,
                                         criterion='squared_error',
                                         max_depth=4,
                                         random_state=2,
                                         # oob_score=True,
                                         n_jobs=-1)

        elif self._method == 'elastic':
            vreg = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9, 0.95],
                                cv=5,
                                fit_intercept=False,
                                n_jobs=-1,
                                tol=1e-3,
                                max_iter=2000,
                                random_state=2,
                                )

        elif self._method == 'boosting':
            vreg = GradientBoostingRegressor(n_estimators=1000,
                                             max_depth=4,
                                             subsample=0.8,
                                             min_samples_split=5,
                                             learning_rate=0.04,
                                             random_state=2,
                                             loss="huber", )

        else:
            raise Exception('Invalid method')


        # Quadratic transform
        # from sklearn.preprocessing import PolynomialFeatures
        # quad = PolynomialFeatures(degree=3)
        # X_train = quad.fit_transform(X_train)
        # X_test = quad.fit_transform(X_test)

        '''
        #plot coefficients
        elastic.fit(X_train, y_train)
        pd.DataFrame([_ for _ in zip(X_train.columns, elastic.coef_)]).to_csv(protein_to_do + 'coef.txt')

        '''

        vreg.fit(x_train, y_train)
        y_train_pred = vreg.predict(x_train)
        y_test_pred = vreg.predict(x_test)

        # End

        # print('train R2 {0}, test R2 {1}, test corr {2}'.format(r2_train, r2_test, corr))

        # Plot scatter plot
        '''
        plt.scatter(y=y_train, x=y_train_pred, label="train")
        plt.scatter(y=y_test, x=y_test_pred, label="test")
        plt.legend(loc=2)
        plt.xlabel("Predicted from RNA-seq")
        plt.ylabel("Actual Protein Conc.")
        plt.title('protein {0} with train R2: {1} and test R2: {2}'.format(protein_to_do,
                                                                           r2_train,
                                                                           r2_test))
        plt.savefig(fname=os.path.join(directory, protein_to_do + '.png'))
        plt.clf()
        #
        '''

        # Write train and test table
        x_train_out = x_train.copy()
        x_train_out['y_train'] = y_train
        x_train_out['y_train_pred'] = y_train_pred
        # X_train_out.to_csv(os.path.join(directory, protein_to_do + '_train.txt'))

        x_test_out = x_test.copy()
        x_test_out['y_test'] = y_test
        x_test_out['y_test_pred'] = y_test_pred
        # X_test_out.to_csv(os.path.join(directory, protein_to_do + '_test.txt'))

        # Mean square error and R2
        if np.std(y_test_pred) == 0:
            corr_test = 0
        else:
            corr_test = np.corrcoef(y_test, y_test_pred)[0][1]

        if np.std(y_train_pred) == 0:
            corr_train = 0
        else:
            corr_train = np.corrcoef(y_train, y_train_pred)[0][1]

        r2_train = r2_score(y_train, y_train_pred)
        r2_test = r2_score(y_test, y_test_pred)
        nrmse = np.sqrt(mean_squared_error(y_test, y_test_pred) / (np.max(y_test) - np.min(y_test)))

        # baseline nrmse by simply guessing the median value of the y_train
        baseline = np.sqrt(mean_squared_error(y_test, np.repeat(np.median(y_train),
                                                                len(y_test))) / (np.max(y_test) - np.min(y_test)))

        metrics_df = pd.DataFrame(data={'corr_train': [corr_train],
                                        'corr_test': [corr_test],
                                        'r2_train': [r2_train],
                                        'r2_test': [r2_test],
                                        'num_obs': x_train.shape[0],
                                        'num_features': x_train.shape[1],
                                        'nrmse': [nrmse],
                                        'baseline_nrmse': [baseline]
                                        },

                                  index=[protein_to_do])

        return {'model': vreg, 'metrics': metrics_df}



