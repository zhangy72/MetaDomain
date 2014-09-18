// Copyright (C) 2011 Yuan Zhang, Yanni Sun.
// You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
// implementation of HmmModel class.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "hmm_model.h"

using namespace std;
const double HmmModel::kInf_ = -32768;
const int HmmModel::kAlphaSize_ = 20; // The number of alphabet, i.e. the number of possible valid amino acid characters
const int HmmModel::kSpecialStateNum_ = 5;// B, C, E, N and J are 4 special states
const int HmmModel::p7P_B = 0;
const int HmmModel::p7P_C = 1;
const int HmmModel::p7P_E = 2;
const int HmmModel::p7P_N = 3; 
const int HmmModel::p7P_J = 4;
const int HmmModel::p7P_LOOP = 0;
const int HmmModel::p7P_MOVE = 1;
const int HmmModel::M_STATE = 0;
const int HmmModel::I_STATE = 1;
const int HmmModel::D_STATE = 2;
const int HmmModel::B_STATE = 3;
const int HmmModel::E_STATE = 4;
const int HmmModel::J_STATE = 5;
const int HmmModel::N_STATE = 6;
const int HmmModel::C_STATE = 7;

const int HmmModel::kTransitionNum_ = 9;
const int HmmModel::M2M = 0;
const int HmmModel::M2I = 1;
const int HmmModel::M2D = 2;
const int HmmModel::I2M = 3;
const int HmmModel::I2I = 4;
const int HmmModel::D2M = 5;
const int HmmModel::D2D = 6;
const int HmmModel::B2M = 7;
const int HmmModel::M2E = 8;

const int HmmModel::kBeginTransitionNum_ = 5; 
const int HmmModel::B2M1 = 0;
const int HmmModel::B2I0 = 1;
const int HmmModel::B2D = 2;
const int HmmModel::I02M1 = 3;
const int HmmModel::I02I0 = 4;

void ConstructHmmModel(ifstream& hmm_file, HmmModel& hmm) {
	if (!hmm_file.is_open()) {
		cout << "HMM file is not open!" << endl;
		hmm_file.close();
		exit(1);
	}
	const int kMaxLineLength = 4096;
	char line[kMaxLineLength];
	while (hmm_file.getline(line, kMaxLineLength)) {
		if (line[0] == '/' && line[1] == '/') {
			break;
		} else {
			vector<string> tokens;
			TokenizeString(line, tokens, " ");
			if (tokens[0] == "LENG") {
				hmm.set_length(atoi(tokens[1].c_str()));
				hmm.insert_emission_ = vector<vector<double> >(hmm.length() + 1, vector<double>(HmmModel::kAlphaSize_ + 1)); // Note that we have I0 state and 'X'
				hmm.match_emission_ = vector<vector<double> >(hmm.length() + 1, vector<double>(HmmModel::kAlphaSize_ + 1)); // M0 is background emissions
				hmm.transition_ = vector<vector<double> >(hmm.length() + 1, vector<double>(HmmModel::kTransitionNum_)); // Transitions start with 1.Transition_[0] is vacant			
				hmm.max_emission_ = vector<double>(hmm.length() + 1); // The average emssion score for each state. 0 is for background			
			} else if (tokens[0] == "NAME") {
				hmm.set_name(tokens[1]);
			} else if (tokens[0] == "ACC") {
				hmm.set_accession(tokens[1]);
			} else if ( tokens[0] == "COMPO") { // This block reads the three lines after the COMPO line inclusively.
				hmm.max_emission_[0] = 0; // This is not useful since other max_emission means the difference between emission and background
				hmm.match_emission_[0][HmmModel::kAlphaSize_] = log(0.05);
				for (int i = 0; i < HmmModel::kAlphaSize_; i++) {
					if (tokens[i + 1] == "*") {
						hmm.match_emission_[0][i] = HmmModel::kInf_;
					} else {
						// Background probability
						hmm.match_emission_[0][i] = 0.0 - atof(tokens[i + 1].c_str()); // Convert to log(probability)
					}
				}		
				hmm_file.getline(line, kMaxLineLength); // I0 insertion scores
				TokenizeString(line, tokens, " ");
				hmm.insert_emission_[0][HmmModel::kAlphaSize_] = log(0.05);
				for (int i = 0; i < HmmModel::kAlphaSize_; i++) {
					if (tokens[i] == "*") {
						hmm.insert_emission_[0][i] = HmmModel::kInf_;
					} else {
						hmm.insert_emission_[0][i] = 0.0 - atof(tokens[i].c_str());
					}
				}
				hmm_file.getline(line, kMaxLineLength); // Begin transition line
				hmm.begin_transition_ = vector<double>(HmmModel::kBeginTransitionNum_);
				for (int i = 0; i < HmmModel::kBeginTransitionNum_; i++) { // Only consider the first 5 fields. The last two are for deletion state,
																			// which does not exist in the beginning stage
					hmm.begin_transition_[i] = 0.0 - atof(tokens[i].c_str());
				}
				for (int i = 1; i < hmm.length(); i ++) {
					hmm_file.getline(line, kMaxLineLength); // Match state emissions
					TokenizeString(line, tokens," ");
					hmm.max_emission_[i] = HmmModel::kInf_;
					hmm.match_emission_[i][HmmModel::kAlphaSize_] = log(0.05);	
					for (int j = 0; j < HmmModel::kAlphaSize_; j++) {
						hmm.match_emission_[i][j] = 0.0 - atof(tokens[j + 1].c_str()); // The first column is the state index
						// Give the maximum score difference between emssion and background to max_emission
						if (hmm.max_emission_[i] < hmm.match_emission_[i][j] - hmm.match_emission_[0][j]) {
							hmm.max_emission_[i] = hmm.match_emission_[i][j]  - hmm.match_emission_[0][j];
						}
					}	
					hmm_file.getline(line, kMaxLineLength); // Insertion state emissions
					TokenizeString(line, tokens," ");
					hmm.insert_emission_[i][HmmModel::kAlphaSize_] = log(0.05);
					for (int j = 0; j < HmmModel::kAlphaSize_; j++) {
						hmm.insert_emission_[i][j] = 0.0 - atof(tokens[j].c_str());
					}
					hmm_file.getline(line, kMaxLineLength); // Transition scores
					TokenizeString(line, tokens," ");
					for (int j = 0; j < 7; j++) {
						hmm.transition_[i][j] = 0.0 - atof(tokens[j].c_str());
					}
					hmm.transition_[i][HmmModel::M2E] = 0.0; // No penalty to from match state to END state
					hmm.transition_[i][HmmModel::B2M] = 0.0; // No penalty to BEGIN state to match state
					/*
					// Print out match emssion scores
					for (int j = 0; j <= HmmModel::kAlphaSize_; j++) {
						cout << hmm.match_emission_[i][j] << " "; // Convert to log(probability)
					}
					cout << "max:" << hmm.max_emission_[i] << endl;	*/		
				}
				// Consider the last block, which has the transition from the last state to END
				hmm_file.getline(line, kMaxLineLength); // Match state emissions
				TokenizeString(line, tokens," ");
				// Give 'X' average score
				hmm.match_emission_[hmm.length()][HmmModel::kAlphaSize_] = log(0.05);
				for (int j = 0; j < HmmModel::kAlphaSize_; j++) {
					hmm.match_emission_[hmm.length()][j] = 0.0 - atof(tokens[j + 1].c_str());
					if (hmm.max_emission_[hmm.length()] < hmm.match_emission_[hmm.length()][j] - hmm.match_emission_[0][j]) {
						hmm.max_emission_[hmm.length()] = hmm.match_emission_[hmm.length()][j] - hmm.match_emission_[0][j];
					}
				}	
				hmm_file.getline(line, kMaxLineLength); // Insertion state emissions
				TokenizeString(line, tokens," ");
				// Insertion score for 'X' should be the average value
				hmm.insert_emission_[hmm.length()][HmmModel::kAlphaSize_] = log(0.05);
				for (int j = 0; j < HmmModel::kAlphaSize_; j++) {
					hmm.insert_emission_[hmm.length()][j] = 0.0 - atof(tokens[j].c_str());
				}
				hmm_file.getline(line, kMaxLineLength); // Transition scores
				TokenizeString(line, tokens," ");
				for (int j = 0; j < 7; j++) {
					if (tokens[j] == "*") {
						hmm.transition_[hmm.length()][j] = HmmModel::kInf_;
					} else {
						hmm.transition_[hmm.length()][j] = 0.0 - atof(tokens[j].c_str());
					}
				}
			}
			else {
				continue;
			}
		} 
	} 
	// Allocate space for special transitions
	hmm.special_transition_ = vector<vector<double> >(HmmModel::kSpecialStateNum_, vector<double>(2));
	hmm.special_transition_[HmmModel::p7P_E][HmmModel::p7P_MOVE] = 0.0;   
	hmm.special_transition_[HmmModel::p7P_E][HmmModel::p7P_LOOP] = HmmModel::kInf_;  
   	hmm.set_nj(0.0);
	// Transition between BEGIN state and match states
	/*float Z = 0.0f;
	float *occ = new float[hmm.length() + 1];
    if ((p7_hmm_CalculateOccupancy(hmm, occ, NULL)) != 0) {
		cout << "B2M score calculation error!" << endl;
		hmm_file.close();
		exit(1);
	}
    for (int k = 1; k <= hmm.length(); k++) {
	  Z += occ[k] * float(hmm.length() - k + 1);
    }
	vector<double> B2M_penalty(hmm.length()); // Print out the BEGIN to match state penalty
	//ofstream out_file("B2M_penalty.txt");
    for (int k = 1; k <= hmm.length(); k++) {
		if (occ[k] <= 0) {
			B2M_penalty[k - 1] = HmmModel::kInf_;
		} else {
			B2M_penalty[k - 1] = log(occ[k] / Z);
		}
    }
	//out_file << DoubleVectorAverage(B2M_penalty.begin(), B2M_penalty.end()) << endl;
	//out_file.close();
	delete[] occ;*/
}
double DoubleVectorAverage(vector<double>::const_iterator begin, vector<double>::const_iterator end) {
	vector<double>::const_iterator it;
	double average = 0;
	int count = 0;
	for (it = begin; it != end; it++) {
			average += *it;
			count++;
	}
	return average / count;
}	
void TokenizeString(const char* line, vector<string>& tokens, const string& delimiters) {
	tokens.clear(); // Initialize tokens because tokens.push_back() will be used.
	string str(line);
    string::size_type last_pos = str.find_first_not_of(delimiters, 0);
    string::size_type pos     = str.find_first_of(delimiters, last_pos);
    while (string::npos != pos || string::npos != last_pos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(last_pos, pos - last_pos));
        // Skip delimiters.  Note the "not_of"
        last_pos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, last_pos);
    }
}
int p7_hmm_CalculateOccupancy(const HmmModel& hmm, float *mocc, float *iocc)
{
  int k;
  const int length = hmm.length();
  mocc[0] = 0.;			                    /* no M_0 state */
  mocc[1] = exp(hmm.begin_transition_[HmmModel::M2I]) + exp(hmm.begin_transition_[HmmModel::M2M]);   /* initialize w/ 1 - B->D_1 */
  for (k = 2; k <= length; k++){
	  mocc[k] = mocc[k-1] * (exp(hmm.transition_[k-1][HmmModel::M2M]) + exp(hmm.transition_[k-1][HmmModel::M2I])) +
		  (1.0-mocc[k-1]) * exp(hmm.transition_[k-1][HmmModel::D2M]);  }

  if (iocc != NULL) {
	  iocc[0] = exp(hmm.begin_transition_[HmmModel::M2I] - hmm.begin_transition_[HmmModel::I2I]);
    for (k = 1; k <= length; k++)
		iocc[k] = mocc[k] * exp(hmm.transition_[k][HmmModel::M2I] - hmm.transition_[k][HmmModel::I2I]);
  }
  return 0;
}
