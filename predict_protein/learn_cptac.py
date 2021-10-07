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
from multiprocessing import Pool, cpu_count
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor, VotingRegressor
from sklearn.linear_model import LinearRegression, ElasticNet, LassoCV, Lasso, ElasticNetCV


class LearnCPTAC(object):
    """
    df
    """

    def __init__(self, cptac_df):

        from get_proteins import GetProtein, GetComplex

        self.df = cptac_df

        self.all_proteomics = [re.sub('_proteomics', "", protein) for protein in self.df.columns if
                               protein.endswith('_proteomics')]

        self.all_transcriptomics = [re.sub('_transcriptomics', "", transcript) for transcript in self.df.columns if
                                    transcript.endswith('_transcriptomics')]

        self.shared_proteins = [protein for protein in self.all_proteomics if protein in self.all_transcriptomics]

        self.tx_to_include = "self"
        self.train_method = "simple"

        self.stringdb = GetProtein()
        self.corumdb = GetComplex()

        pass

    def learn_all_proteins(
            self,
            tx_to_include="self",
            train_method="simple",
    ):
        """
        Wrapper for learning one protein

        :param tx_to_include: transcript to include (self, all, string, or corum)
        :param train_method: simple or voting
        :return:
        """

        self.tx_to_include = tx_to_include
        self.train_method = train_method

        learning_results = []

        for i, protein in enumerate(tqdm.tqdm(self.shared_proteins)):
            learning_result = self.learn_one_protein(protein)

            if learning_result:
                learning_results.append(learning_result)

                if i % 100 == 0:
                    corr_values = [metric[2].corr_test.values[0] for metric in learning_results]
                    r2_values = [metric[2].r2_test.values[0] for metric in learning_results]
                    nrmse = [metric[2].nrmse.values[0] for metric in learning_results]

                    tqdm.tqdm.write('{0}: {1}, r: {2}, R2: {3}, med.r: {4}, med.R2: {5}, med.NRMSE: {6}'.format(
                        i,
                        protein,
                        round(list(learning_result[2].corr_test.values)[0], 3),
                        round(list(learning_result[2].r2_test.values)[0], 3),
                        round(np.median([r for r in corr_values if not np.isnan(r)]), 3),
                        round(np.median([r2 for r2 in r2_values if not np.isnan(r2)]), 3),
                        round(np.median([nr for nr in nrmse if not np.isnan(nr)]), 3),
                    ))

        return learning_results

    def learn_one_protein(
            self,
            protein_to_do,
            returnModel=False,
    ):
        """

        :param protein_to_do:
        :param returnModel:  Return the model rather than the results, for examining coefficients
        :return:
        """
        # return [X_train_out, X_test_out, comb_df]

        y_df = self.df[[protein_to_do + '_proteomics']]
        y_df = y_df.dropna(subset=[protein_to_do + '_proteomics'])

        # Decide how many proteins to use:
        if self.tx_to_include == "self":
            proteins_to_include = [protein_to_do]

        elif self.tx_to_include == "all":
            proteins_to_include = self.shared_proteins

        elif self.tx_to_include == "string":
            string_interactors = self.stringdb.find_interactor(protein_to_do, combined_score=200, max_proteins=1000)
            proteins_to_include = [p for p in string_interactors if p in self.shared_proteins]

        elif self.tx_to_include == "corum":
            corum_interactors = self.corumdb.find_cosubunits(protein_to_do)
            proteins_to_include = [p for p in corum_interactors if p in self.shared_proteins]

        # Join X and Y
        xy_df = self.df[[tx + '_transcriptomics' for tx in proteins_to_include]].join(y_df, how='inner').copy().dropna()

        # Skip proteins with fewer than 20 samples
        if len(xy_df) < 20:
            return []

        # Mean impute the missing values
        # en_XY_imp = SimpleImputer(missing_values=np.nan, strategy='mean').fit_transform(en_XY)
        # en_XY = pd.DataFrame(en_XY_imp, columns=en_XY.columns, index=en_XY.index)
        # End

        # Do train-test split
        x = xy_df.iloc[:, :-1]  # .values
        y = xy_df.iloc[:, -1]  # .values
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=2)
        # End

        # PCA the train and test
        # from sklearn.decomposition import PCA
        # pca = PCA(n_components=15)
        # X_train_pca = pca.fit_transform(X_train)
        # X_test_pca = pca.transform(X_test)
        # End

        if self.train_method == 'simple':
            vreg = LinearRegression(n_jobs=-1)

        elif self.train_method == 'voting':
            # Train model

            elastic = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1],
                                   cv=5,
                                   fit_intercept=False,
                                   n_jobs=-1,
                                   tol=1e-3,
                                   )

            # Shallow Forest was 100 estimators 5 depth. Deep is 1000 estimators 10 depth.
            forest = RandomForestRegressor(n_estimators=1000,
                                           criterion='mse',
                                           max_depth=10,
                                           random_state=1,
                                           n_jobs=-1)

            vreg = VotingRegressor(estimators=[('en', elastic), ('rf', forest), ])

        elif self.train_method == 'forest':
            vreg = RandomForestRegressor(n_estimators=100,
                                         criterion='mse',
                                         max_depth=5,
                                         random_state=1,
                                         n_jobs=-1)

        elif self.train_method == 'elastic':
            vreg = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1],
                                cv=5,
                                fit_intercept=False,
                                n_jobs=-1,
                                )

        # slr = Lasso()

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

        if returnModel:
            return (vreg)
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
            corr = 0
        else:
            corr = round(np.corrcoef(y_test, y_test_pred)[0][1], 4)
        r2_test = round(r2_score(y_test, y_test_pred), 4)
        nrmse = round(np.sqrt(mean_squared_error(y_test, y_test_pred) / (np.max(y_test) - np.min(y_test))), 4)
        metrics_df = pd.DataFrame(data={'corr_test': [corr],
                                        'r2_test': [r2_test],
                                        'num_proteins': [len(xy_df)],
                                        'nrmse': [nrmse],
                                        },

                                  index=[protein_to_do])

        return [round(x_train_out, 5), round(x_test_out, 5), round(metrics_df, 5)]



