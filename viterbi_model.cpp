// Copyright (C) 2011 Yuan Zhang, Yanni Sun.
// You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
// implemenation of Alignment class and Viterbi algorithm.
#include <vector>
#include <algorithm>
#include <stack>
#include <cstring>
#include "viterbi_model.h"
using namespace std;
Alignment::Alignment(bool reverse_complement, const int& length_threshold, const double& rate)
:reverse_complement_(reverse_complement), 
 length_threshold_(length_threshold),
 rate_ (rate),
 hmm_acc_("none")
{
}
int Codon2AminoAcidIndex(const char* codon);
// In this function all indices should begin with 0 to be consistent
void Viterbi(HmmModel& hmm, DNASeq& dna, Alignment& alignment) {
	const int M = hmm.length(); // M is the number of states
	const int L = dna.length();
	char* seq = new char[L + 2];
	seq[0] = '$';
	if (!alignment.reverse_complement()) { // Original sequence
		alignment.set_seq_name(dna.seq_name());
    strncpy(seq + 1,dna.seq().c_str(), dna.length());
	} else {
		alignment.set_seq_name(dna.seq_name());
		strncpy(seq + 1, dna.ReverseComplement().c_str(), dna.length());
	}
	p7_ReconfigLength(hmm, L / 3); // Configure some special parameters for hmm model.
	//The following matrix for DP begins with 1. The first position is unused.
	vector<vector<double> > mat_M(M + 1, vector<double>((L + 1), HmmModel::kInf_)); // Score matrix for match staes
	vector<vector<double> > mat_I(M + 1, vector<double>((L + 1), HmmModel::kInf_)); // Score matrix for insertion states
	vector<vector<double> > mat_D(M + 1, vector<double>((L + 1), HmmModel::kInf_)); // Score matrix for deletion states
	vector<vector<double> > mat_X(4, vector<double>(L + 1, HmmModel::kInf_)); // B, C, E and N states
	vector<vector<vector<int> > > flag(3, vector<vector<int> >(M + 1, vector<int>(L + 1, -1))); // Flags for M, I, D and G states
	vector<vector<int> >E_state_flag(L + 1, vector<int>(2, -1)); //  E_state_flag[i][0] stores the value of j which gives largest END score; E_state_flag[i][1] indicates whether E state leads to
																										// largest score 
	vector<double> background = hmm.match_emission_[0]; // Background emission scores

	// Do some initialization
	mat_X[HmmModel::p7P_N][0] = mat_X[HmmModel::p7P_N][1] = mat_X[HmmModel::p7P_N][2] = 0; // S -> N, p = 1
	mat_X[HmmModel::p7P_B][0] = mat_X[HmmModel::p7P_B][1] = mat_X[HmmModel::p7P_B][2] 
		= hmm.special_transition_[HmmModel::p7P_N][HmmModel::p7P_MOVE]; // S->N->B, no N-tail 
	char codon[4];
	codon[3] = '\0';
	//amino_acid_index is used to store the translated amino acid index
	int amino_acid_index; 
	// Begin the main DP
	for (int i = 3; i < L + 1; i++) {
		codon[0] = seq[i - 2];
		codon[1] = seq[i - 1];
		codon[2] = seq[i];
		amino_acid_index = Codon2AminoAcidIndex(codon);
		double max, temp_score;
		for (int j = 1; j < M + 1; j++) {
			// Match state
			// The first block contains cases without sequencing errors in the current codon
			max = mat_X[HmmModel::p7P_B][i - 1] + hmm.transition_[j - 1][HmmModel::B2M] + hmm.match_emission_[j][amino_acid_index]
				- background[amino_acid_index];  // Case 0: BEGIN state
			flag[HmmModel::M_STATE][j][i] = 0;
			temp_score = mat_M[j - 1][i - 3] + hmm.transition_[j - 1][HmmModel::M2M] + hmm.match_emission_[j][amino_acid_index]
				- background[amino_acid_index]; // Case 1: previous match state without error
			if (temp_score > max) {
				max = temp_score;
				flag[HmmModel::M_STATE][j][i] = 1; 
			}
			temp_score = mat_I[j - 1][i - 3] + hmm.transition_[j - 1][HmmModel::I2M] + hmm.match_emission_[j][amino_acid_index]
				- background[amino_acid_index];
			if (temp_score > max) {
				max = temp_score;
				flag[HmmModel::M_STATE][j][i] = 2; // Case 2: previous insertion state without error
			}
			temp_score = mat_D[j - 1][i - 3] + hmm.transition_[j - 1][HmmModel::D2M] + hmm.match_emission_[j][amino_acid_index]
				- background[amino_acid_index];
			if (temp_score > max) {
				max = temp_score;
				flag[HmmModel::M_STATE][j][i] = 3; // Case 3: previous deletion state without error
			}
			mat_M[j][i] = max;
			//Insertion state
			max = mat_M[j][i - 3] + hmm.transition_[j][HmmModel::M2I] + hmm.insert_emission_[j][amino_acid_index] 
				- background[amino_acid_index];
			flag[HmmModel::I_STATE][j][i] = 0; // Case 0: from previous match state
			temp_score = mat_I[j][i] + hmm.transition_[j][HmmModel::I2I] + hmm.insert_emission_[j][amino_acid_index] 
				- background[amino_acid_index];
			if (temp_score > max) {
				max = temp_score;
				flag[HmmModel::I_STATE][j][i] = 1; // Case 1: a insertion state loop
			}
			mat_I[j][i] = max;
			// Deletion state
			max = mat_M[j - 1][i] + hmm.transition_[j - 1][HmmModel::M2D];
			flag[HmmModel::D_STATE][j][i] = 0; // Case 0: from previous match state
			temp_score = mat_D[j - 1][i] + hmm.transition_[j - 1][HmmModel::D2D];
			if (temp_score > max) {
				max = temp_score;
				flag[HmmModel::D_STATE][j][i] = 1; // Case 1: from previous deletion state
			}
			mat_D[j][i] = max;
			// E state. Here do not consider the last state to be a deletion state
			max = mat_X[HmmModel::p7P_E][i];
			temp_score = mat_M[j][i];
			if (temp_score > max) {
				E_state_flag[i][0] = j;
				max = temp_score;
			}
			mat_X[HmmModel::p7P_E][i] = max;
		}
		// Other special states
		double escore = mat_X[HmmModel::p7P_E][i];
		max = mat_X[HmmModel::p7P_C][i - 1] + hmm.special_transition_[HmmModel::p7P_C][HmmModel::p7P_LOOP];
		E_state_flag[i][1] = 0; // 0 indicates that the sequence has not ended
		temp_score = mat_X[HmmModel::p7P_E][i] + hmm.special_transition_[HmmModel::p7P_E][HmmModel::p7P_MOVE];
		if (temp_score > max) {
			max = temp_score;
			E_state_flag[i][1] = 1; // This is the place that an END state leads to better score
		}
		mat_X[HmmModel::p7P_C][i] = max;
		mat_X[HmmModel::p7P_N][i] = mat_X[HmmModel::p7P_N][i-1] + hmm.special_transition_[HmmModel::p7P_N][HmmModel::p7P_LOOP];
		mat_X[HmmModel::p7P_B][i] = mat_X[HmmModel::p7P_N][i] + hmm.special_transition_[HmmModel::p7P_N][HmmModel::p7P_MOVE];
	}
	double score = mat_X[HmmModel::p7P_C][L] + hmm.special_transition_[HmmModel::p7P_C][HmmModel::p7P_MOVE];
	score = score / log(2.0f);
	string out_seq;
	int begin_state;
	int end_state;
	int begin_base;
	int end_base;
	TraceBack(flag, E_state_flag, L, M, seq, out_seq, begin_state, end_state, begin_base, end_base);

	// Begin to set alignment values
	alignment.set_score(score);
	alignment.set_DNA_length(dna.length());
	alignment.set_begin_state(begin_state);
	alignment.set_end_state(end_state);
	alignment.set_begin_base(begin_base);
	alignment.set_end_base(end_base);
	alignment.set_align_length((out_seq.length() - 1) / 3); // The aligned length is in amino acid.
	alignment.set_hmm_acc(hmm.accession());
	alignment.set_hmm_length(hmm.length());
	alignment.set_hmm_name(hmm.name());
	alignment.set_upper_bound(CalculateUpperBound(hmm, begin_state, end_state, dna.length() / 3));
	alignment.set_score_threshold();
  delete[] seq;
}

void TraceBack(const vector<vector<vector<int> > >& flag, const vector<vector<int> >& E_state_flag,  int L, int M, const char* seq, string& out_seq,
	int& begin_state, int& end_state, int& begin_base, int& end_base) {
	end_base = L;
	for (int i = L; i > 0; i--) {
		if(E_state_flag[i][1] == 1)
		{
			end_base = i;
			break;
		}
	}
	end_state = E_state_flag[end_base][0];
	out_seq= string(seq, seq + end_base + 1); // The output sequence. Begin base not known yet. Index begins with 1
	int choice = 0; // choice[0]: match state; choice[1]: insertion state; choice[2]: deletion state; choice[3]: G state; choice[4]:BEGIN
	stack<string> path; // The traceback path
	int index1 = end_base;
	int index2 = end_state;
	const int kLabelLength = 10; // The length of each label of the DNA base
	char label[kLabelLength];
	char buff[kLabelLength];
	// Handle end base
	strcpy(label,"E");
	sprintf(buff,"%d",end_base);
	strcat(label,buff);
	path.push(label);
	while (index1 >= 0 && index2 >= 0) {
			if(choice == 0) {
				switch(flag[0][index2][index1]) {
					case 0: // B->M
						strcpy(label,"M");
						sprintf(buff,"%d",index2);
						strcat(label,buff);
						path.push(label);
						path.push(label);
						path.push(label);
						begin_base = index1 - 2;
						begin_state = index2;
						sprintf(buff,"%d",begin_base);
						strcpy(label, "B");
						strcat(label, buff); 
						path.push(label);
						if (begin_base != 1) {
							out_seq.erase(1, begin_base - 1);
						}
						choice = 4;
						break;
					case 1: // M->M
						strcpy(label,"M");
						sprintf(buff,"%d",index2);
						strcat(label,buff);
						path.push(label);
						path.push(label);
						path.push(label);
						index1 -= 3;
						index2--;
						choice = 0;
						continue;
					case 2: //I->M
						strcpy(label,"M");
		     			sprintf(buff,"%d",index2);
						strcat(label,buff);
						path.push(label);
						path.push(label);
						path.push(label);
						index1 -= 3;
						index2--;
						choice = 1;
						continue;
					case 3: //D->M
						strcpy(label,"M");
		     			sprintf(buff,"%d",index2);
						strcat(label,buff);
						path.push(label);
						path.push(label);
						path.push(label);
						index1 -= 3;
						index2--;
						choice = 2;	
			}
		}
		if (choice == 1)
		{	
			switch (flag[1][index2][index1])
			{
				case 0:
					strcpy(label, "I");
					sprintf(buff, "%d", index2);
					strcat(label, buff);	
					path.push(label);
					path.push(label);	
					path.push(label);
					index1 -= 3;
					choice = 0;
					continue;
				case 1:
					strcpy(label, "I");
					sprintf(buff, "%d", index2);
					strcat(label, buff);	
					path.push(label);
					path.push(label);
					path.push(label);
					index1 -= 3;
					choice = 1;
					continue;
			}
		}
		if (choice == 2)
		{
			switch (flag[2][index2][index1]) {
				case 0:
					strcpy(label, "D");
					sprintf(buff, "%d", index2);
					strcat(label, buff);	
					path.push(label);
					index2--;
					choice = 0;
					continue;
				case 1:
					strcpy(label, "D");
					sprintf(buff, "%d", index2);
					strcat(label, buff);	
					path.push(label);
					index2--;
					choice = 2;
					continue;
			}
		}
		if (choice == 4)
		{
			break;
		}		
	}
	// Print out the path
       /*	cout << "The path is:" << endl;
	while (!path.empty()) {
		cout <<  path.top() << " ";
		path.pop();
	}
	cout << endl;*/
}
// Calculate the average score for a given interval in the HMM model
double CalculateUpperBound(const HmmModel& hmm, const int begin_state, const int end_state, const int read_length) {
	double score = 0;
	int real_begin, real_end; // These two states adjust the position of the sub-model
	if (begin_state + read_length - 1 <= hmm.length()) {
		real_begin = begin_state;
		real_end = begin_state + read_length - 1;
	} else if (end_state - read_length > 0) {
		real_begin = end_state + 1 - read_length;
		real_end = end_state;
	} else {
		real_begin = begin_state;
                real_end = end_state;
	}
	for (int i = real_begin; i < real_end; i++) {
		score = score + hmm.max_emission_[i] + hmm.transition_[i][0] ;
	}
	score += hmm.max_emission_[real_end]; // average emission minus average background
	//cout << "begin_state: " << begin_state << " end_state:" << end_state << endl;
	return score;
}
// This function is called in main function due to L
int p7_ReconfigLength(HmmModel& hmm, int seq_length)
{
  double ploop, pmove;  
  /* Configure N,J,C transitions so they bear L/(2+nj) of the total
   * unannotated sequence length L. 
   */
  pmove = (2.0f + hmm.nj()) / ((double) seq_length + 2.0f + hmm.nj()); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;
  hmm.special_transition_[HmmModel::p7P_N][HmmModel::p7P_LOOP] =  hmm.special_transition_[HmmModel::p7P_C][HmmModel::p7P_LOOP]
	= hmm.special_transition_[HmmModel::p7P_J][HmmModel::p7P_LOOP] = log(ploop);
  hmm.special_transition_[HmmModel::p7P_N][HmmModel::p7P_MOVE] =  hmm.special_transition_[HmmModel::p7P_C][HmmModel::p7P_MOVE] 
	= hmm.special_transition_[HmmModel::p7P_J][HmmModel::p7P_MOVE] = log(pmove);
  return 0;
}

 double HiddenCodon2AminoAcidIndex(char* codon, int missing_position, vector<double> background, vector<double> match_emission, char& missing_base) {
	 
	 double max = HmmModel::kInf_; // max_missing_base_score
	 double temp; 
	 int guess_amino_acid_index; // Guess a base
	 int max_amino_acid_index; // The base which maximizes the emission score minus background score
	 // Begin to guess the missing base
	 codon[missing_position] = 'A';
	 guess_amino_acid_index = Codon2AminoAcidIndex(codon);
	 temp = match_emission[guess_amino_acid_index] - background[guess_amino_acid_index];
	 if (temp > max) {
		 max = temp;
		 max_amino_acid_index = guess_amino_acid_index;
		 missing_base = 'A';
	 }
	 codon[missing_position] = 'C';
	 guess_amino_acid_index = Codon2AminoAcidIndex(codon);
	 temp = match_emission[guess_amino_acid_index] - background[guess_amino_acid_index];
	 if (temp > max) {
		 max = temp;
		 max_amino_acid_index = guess_amino_acid_index;
		 missing_base = 'C';
	 }
	 codon[missing_position] = 'G';
	 guess_amino_acid_index = Codon2AminoAcidIndex(codon);
	 temp = match_emission[guess_amino_acid_index] - background[guess_amino_acid_index];
	 if (temp > max) {
		 max = temp;
		 max_amino_acid_index = guess_amino_acid_index;
		 missing_base = 'G';
	 }
	 codon[missing_position] = 'T';
	 guess_amino_acid_index = Codon2AminoAcidIndex(codon);
	 temp = match_emission[guess_amino_acid_index] - background[guess_amino_acid_index];
	 if (temp > max) {
		 max = temp;
		 max_amino_acid_index = guess_amino_acid_index;
		 missing_base = 'T';
	 }
	 return max;
 }
