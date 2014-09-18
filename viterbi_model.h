// Copyright (C) 2011 Yuan Zhang, Yanni Sun.
// You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
// header file of Alignment class and Viterbi algorithm

#ifndef VITERBI_MODEL_H_
#define VITERBI_MODEL_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cmath>
#include "hmm_model.h"
#include "DNA_seq.h"

using namespace std;
class Alignment {
public: 
	Alignment(bool reverse_complement, const int& length_threshold, const double& rate); // These two parameters are set in main
	string seq_name() const {return seq_name_;} // The input sequence name
	bool reverse_complement() const {return reverse_complement_;}
	string hmm_name() const {return hmm_name_;}
	string hmm_acc() const {return hmm_acc_;}
	double score() const { return score_;}
	double score_threshold() const {return score_threshold_;}
	int hmm_length() const {return hmm_length_;}
	int begin_state() const {return begin_state_;}
	int end_state() const {return end_state_;}
	int begin_base() const {return begin_base_;}
	int end_base() const {return end_base_;}
	int DNA_length() const {return DNA_length_;}
	int align_length() const {return align_length_;}
	void set_seq_name(const string& seq_name) {seq_name_ = seq_name;}
	void set_score(const double& score) {score_ = score;}
	void set_upper_bound(const double& upper_bound) {upper_bound_ = upper_bound;}
	void set_DNA_length(const int& length) {DNA_length_ = length;}
	void set_align_length(const int& length) {align_length_ = length;}
	void set_hmm_length(const int& length) {hmm_length_ = length;}
	void set_begin_state(const int& begin_state) {begin_state_ = begin_state;}
	void set_end_state(const int& end_state) {end_state_ = end_state;}
	void set_begin_base(const int& begin_base) {begin_base_ = begin_base;}
	void set_end_base(const int& end_base) {end_base_ = end_base;}
	void set_hmm_name(const string& name) {hmm_name_ = name;}
	void set_hmm_acc(const string& acc) {hmm_acc_ = acc;}
	void set_score_threshold() {score_threshold_ = upper_bound_ * rate_;}
  void set_reverse_complement(bool reverse_complement) {reverse_complement_ = reverse_complement;}
private:
	string seq_name_;
	bool reverse_complement_; // False indicate that the sequence used is the original one
	double score_; // Alignment score
	double score_threshold_; // The alignment score threshold desinged for the interval of alignment state
	double rate_; // The percentage of the upper bound used as the threshold
	double upper_bound_; // The upper bound of alignment score for the aligned interval
	int length_threshold_; // The threshold of the alignment length
	int DNA_length_; // The length of the query DNA sequence
	int align_length_; // Aligned sequence length. The scale is amino acid
	int hmm_length_; // The length of the pfam domain. The scale is amino acid
	int begin_state_; // begin state of the alignment. The scale is amino acid
	int end_state_; // end state of the alignment
	int begin_base_; // The index of the beginning base of the original DNA sequence. 1 means starts from the first base
	int end_base_;
	string hmm_name_; // hmm domain name 
	string hmm_acc_; // hmm accession number if available
};
void Viterbi(HmmModel& hmm, DNASeq& dna, Alignment& alignment);
int p7_ReconfigLength(HmmModel& hmm, int seq_length);
double HiddenCodon2AminoAcidIndex(char* codon, int missing_position, vector<double> background, vector<double> match_emission, char& missing_base);
void TraceBack(const vector<vector<vector<int> > >& flag, const vector<vector<int> >& E_state_flag,  int L, int M, const char* seq, string& out_seq,
	int& begin_state, int& end_state, int& begin_base, int& end_base);
double CalculateUpperBound(const HmmModel& hmm, const int begin_state, const int end_state, const int read_length);
#endif