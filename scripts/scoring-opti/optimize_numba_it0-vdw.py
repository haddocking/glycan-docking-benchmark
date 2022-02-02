#!/usr/bin/python

# Returns an OPTIMIZATION VECTOR: Eelec, Evdw, Edesolv, Eair, BSA

from __future__ import division

from __future__ import with_statement

import os
import numpy as np
from scipy.integrate import simps
import itertools
import glob
import sys
from numba import jit
from sys import stderr


def print_err(*args, **kwargs):
	sys.stderr.write(' '.join(map(str,args)) + '\n')
	

def frange(start, stop, step):

	'''Returns a list of floats given a start, stop and step'''

	weights = []
	x = start
	
	while x <= stop:
		weights.append(round(x,3))
		x += step

	#combinations = list(itertools.combinations_with_replacement(weights, 5))
	#return np.array(combinations)
	return weights


def parse_scores(scores_file):

	'''Parses the scores file (X.scores). It DOES NOT parse the identifier (csv)'''

	with open(scores_file) as inp:
		matrix = np.loadtxt(inp, delimiter=',', usecols=(1, 2, 3, 4, 5, 6, 7))
    
	return matrix
	

@jit(nopython=True)
def calculate_recall(irmsd, topN):

	'''Calculates the recall for a given top (topN) for a list of i-RMSD'''
	
	criteria = 3.0
	N_acceptable = 0
	T_acceptable = 0
	
	for r in list(irmsd):
		if r <= criteria:
			T_acceptable += 1
		else:
			continue
	
	for r in list(irmsd[0:topN]):
		if r <= criteria:
			N_acceptable += 1
		else:
			continue
	
	if N_acceptable == topN or N_acceptable == T_acceptable:
		recall = 1.0
	elif T_acceptable != topN:
		recall = N_acceptable / min(topN,T_acceptable)
	else:
		recall = 0.0
	
	return recall


@jit(nopython=True)
def calculate_hits(irmsd, topN):

	'''Calculates the number of hits (good models according to a given criteria). It returns 1 if 
	#hits > minimum. 0 if #hits < minimum. Code for clustering purposes'''

	N_acceptable = 0
	criteria = 3.0
	minimum = 1

	for r in list(irmsd[:topN]):
		if r <= criteria:
			N_acceptable += 1
		
		if N_acceptable >= minimum:
			x = 1
		else:
			x = 0
	
	return x
			

@jit(nopython=True)
def optimize_weights(energy_terms, matrix, combinations):

	'''MAIN FUNCTION. You may optimize for recall, AUC or min number of good models. NOTE: only 4 terms taken
		into consideration'''

	#nr_models = [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 120, 140, 160, 180, 200, 225, 250, 275, 300, 325, 350, 375, 400]
	#auc = np.zeros((len(combinations),))
	#recall_200 = np.zeros((len(combinations),))
	#recall_400 = np.zeros((len(combinations),))
	
	cases_1 = np.zeros((len(combinations),))
	cases_5 = np.zeros((len(combinations),))
	
	for n, combi in enumerate(combinations):
		#recall = []
		#energy_terms = matrix[:, [2,3,4,5,6]]
		irmsd = matrix[:,0]
		hscore = (energy_terms[:,0] * combi[0] + energy_terms[:,1] * combi[1] + energy_terms[:,2] * combi[2] + \
			energy_terms[:,3] * combi[3] + energy_terms[:,4] * combi[4])
		new = np.column_stack((irmsd, hscore))
		sorted = new[np.argsort(new[:, 1])]	
		sorted_irmsd = sorted[:, 0]
		#for nr in nr_models:
		#	recall.append(calculate_recall(sorted_irmsd, nr))
		
		cases_1[n] = calculate_hits(sorted_irmsd, 1)
		cases_5[n] = calculate_hits(sorted_irmsd, 5)
		
		
		'''This piece calculates the AUC for a given set of points (recalls)'''
		#i = 1
		#total = recall[0] + recall[-1]
		
		#for r in recall[1:-1]:
		#	if i%2 == 0:
		#		total += 2*r
		#	else:
		#		total += 4*r
		#	i += 1
		#auc[n] = total * (1/3.0)
		
	#return auc, recall_200, recall_400
	return cases_1, cases_5
			

	
if __name__ == '__main__':

	#combinations = frange(-1, 1, 0.1)
	Eelec = frange(0, 2, 0.1) #20	
	Edesolv = frange(0, 2, 0.1) #20
	Evdw = [1] * 20 #20  
	Eair = frange(-2, 2, 0.2) #20
	BSA = frange(-1, 0, 0.05) #20
	#Edesolv = [2.0] * 20
	#Eair = [0.0] * 20
	#BSA = [-0.2] * 20
	Etotal = [Eelec, Evdw, Edesolv, Eair, BSA]
	combinations = list(itertools.product(*Etotal))
	#combinations.append(tuple([1.9, 0.0, 1.9, -0.05]))
	#print len(combinations)
	#nr_models = [1, 10, 25, 50, 100, 200, 400]
	#matrix = parse_scores('test.txt')
	
	file_list = [each for each in os.listdir(sys.argv[1]) if each.endswith('.scores')]
	#bench_auc = []
	#bench_recalls200 = []
	#bench_recalls400 = []
	bench_1 = []
	bench_5 = []
	
	for file in file_list:
		matrix = parse_scores(sys.argv[1] + file)
		print(file)
		energy_terms = matrix[:, [2,3,4,5,6]]
		#auc, recall_200, recall_400 = optimize_weights(energy_terms, matrix, combinations)
		#bench_auc.append(auc)
		#bench_recalls200.append(recall_200)
		#bench_recalls400.append(recall_400)
		cases_1, cases_5 = optimize_weights(energy_terms, matrix, combinations)
		bench_1.append(cases_1)
		bench_5.append(cases_5)
		sys.stderr.write(os.path.basename(file)[:7]+"\n")
	
	#matrix_auc = np.array(bench_auc)
	#matrix_recall200 = np.array(bench_recalls200)
	#matrix_recall400 = np.array(bench_recalls400)
	
	matrix_1 = np.array(bench_1)
	matrix_5 = np.array(bench_5)
	
	output = open("all_combi.txt", 'a+')
	for n, combi in enumerate(combinations):
		#mean_auc = np.mean(matrix_auc[:,n])
		#mean_recall200 = np.mean(matrix_recall200[:,n])
		#mean_recall400 = np.mean(matrix_recall400[:,n])
		N_1 = np.sum(matrix_1[:,n])
		N_5 = np.sum(matrix_5[:,n])
		
		#print "STEP:", n, "AUC:", '%f' % mean_auc, "RECALL_200:", '%f' % mean_recall200, "RECALL_400:", '%f' % mean_recall400, "OPTIMIZATION_VECTOR:", combi
		output.write("STEP: {} N_1: {} N_5: {} OPTIMIZATION_VECTOR: {}".format(n, N_1, N_5, combi) + '\n')
	output.close()
		
	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
