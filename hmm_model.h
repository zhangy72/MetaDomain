// Copyright (C) 2011 Yuan Zhang, Yanni Sun.
// You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
// header file of HmmModel class.

#ifndef HMM_MODEL_H_
#define HMM_MODEL_H_
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>

using namespace std;
class HmmModel { // All the numbers are scores. No probabilities.
public:
	void set_length(const int length) {length_ = length;}
	void set_accession(const char* accession) {accession_ = string(accession);}
	void set_accession(const string& accession) {accession_ = accession;}
	void set_name(const char* name) {name_ = string(name);}
	void set_name(const string name) {name_ = name;}
	void set_nj(const double nj) {nj_ = nj;}
	const double nj() const {return nj_;}
	const int length() const {return length_;}
	const string accession() const {return accession_;}
	const string name() const {return name_;}
	// Static const variables
	static const double kInf_; // The possibly lowest score
	static const int kAlphaSize_; // The number of alphabet, i.e. the number of possible valid amino acid characters
	static const int p7P_N, p7P_E, p7P_C, p7P_B, p7P_J;
	static const int p7P_LOOP, p7P_MOVE;
	static const int M_STATE, I_STATE, D_STATE, B_STATE, E_STATE, J_STATE, N_STATE, C_STATE;
	static const int kTransitionNum_;
	static const int kSpecialStateNum_; // There are 4 special states
	static const int M2M, M2I, M2D, I2M, I2I, D2M, D2D, B2M, M2E;
	static const int kBeginTransitionNum_; // Transitions at the beginning stage of the HMM model
	static const int B2M1, B2I0,B2D, I02M1, I02I0;

	vector<vector<double> > match_emission_; // Emission scores for all match states. Match states are from 0 to length_. M0 is for background. 
											//Amino acids are from 0 to 20
										   // The last emssion score is for 'X'
	vector<vector<double> > insert_emission_; // Emission scores of 20 amino acids for all insertion states.
	vector<vector<double> > transition_; // Transition index begins with 1. 0 in the vector is vacant. Transition scores in the order of M->M,M->I,M->D,I->M,I->I,D->M,D->D,B->M,M->E
	vector<double> begin_transition_;    // B->M1,B->I0,B->D,I0->M1,I0->I0
	vector<vector<double> > special_transition_;			// transition scores for 5 special states, p7P_B, p7P_C, p7P_E, p7P_N and p7P_J.
														// Each state includes p7P_LOOP and p7P_MOVE
	vector<double> max_emission_; // The average emission score for each match state

private:
	int length_;       
	string accession_; // Accession number of the Pfam domain
	string name_; // Name of the Pfam domain
	double nj_;	// Expected number of uses of J states
};
void ConstructHmmModel(ifstream& in_file, HmmModel& hmm); // Set the parameters which have not been set in constructor
int p7_hmm_CalculateOccupancy(const HmmModel& hmm, float *mocc, float *iocc);
double DoubleVectorAverage(vector<double>::const_iterator begin, vector<double>::const_iterator end);
void TokenizeString(const char* line, vector<string>& tokens, const string& delimiters);// Tokenize the input line and output a vector of all separate strings
#endif