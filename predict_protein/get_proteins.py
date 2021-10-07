# -*- coding: utf-8 -*-

""" get interacting protein from a downloaded String DB file.
The String DB file was modified by an R script  turn the Ensembl peptide ID into gene name."""

import pandas as pd


class GetProtein(object):
    """

    """

    def __init__(self):
        self.df = pd.read_csv('data/string/tidy_stringdb_homosapiens_250.txt',
                              sep='\t')
        self.df = self.df.sort_values('combined_score')

        pass

    def find_interactor(self,
                        bait,
                        combined_score=500,
                        max_proteins=100
                        ):
        """

        :param bait:
        :param combined_score: Minimum STRINGdb score to fetch
        :param max_proteins:  Maximum number of entries with max combined_score to return
        :return:
        """

        #interactors = [row[2] for row in self.df.itertuples() if row[1] == bait and row[3] >= combined_score]

        interactors = self.df[(self.df.p1 == bait) & (self.df.combined_score > combined_score)]

        interactors = interactors.sort_values('combined_score', ascending=False).p2

        return [bait] + list(interactors)[:max_proteins]


class GetComplex(object):
    """
    Get interacting complex from a downloaded CORUM file. The CORUM DB file was modified by an R script to turn
    Uniprot into GeneName
    """

    def __init__(self):
        self.corum = pd.read_csv('data/corum/tidy_corum_homosapiens.txt', sep='\t')
        pass

    def find_cosubunits(self,
                        bait):
        """
        Find other proteins that belong to the same complex as the bait
        :param bait:
        :return:
        """

        # Get all the complexes associated with the query protein
        complexes = list(self.corum[(self.corum.SubunitGN == bait)]['ComplexID'].values)

        # Get all the subunits associated
        cosubunits = list(self.corum[(self.corum.ComplexID.isin(complexes))].SubunitGN.values)

        return list(set(cosubunits + [bait]))