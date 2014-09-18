// Copyright (C) 2011 Yuan Zhang, Yanni Sun.
// You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
// header file of DNASeq class.

#ifndef DNASEQ_H_
#define DNASEQ_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

class DNASeq {
public:
	DNASeq(const string& seq_name, const string& seq);
	string seq_name() const {return seq_name_;}
	string seq() const {return seq_;}
	int length() const {return length_;}
	void set_corrected_seq(string out_seq) {corrected_seq_ = out_seq;}
	string corrected_seq() const {return corrected_seq_;} 
	vector<string> protein_seqs(); // Return 6-frame translations of seq_
	bool error_model() const {return error_model_;}
	bool set_error_model(bool error_model) {error_model_ = error_model;}
  string ReverseComplement() const;
	friend int Codon2AminoAcidIndex(const char* codon);
	friend char Codon2AminoAcid(const char* codon);
private:
	string seq_name_; 
	string seq_; 
	int length_;
	string corrected_seq_; // The output sequence after error correction
	bool error_model_;
	static const int kCodonNum_;
	static const char* kCodonTable_; // Amino acids for 64 possible codons plus 'X' for codons wtih unexpected characters other than A, C, G and T
	static const int kCodonIndexTable_[65]; // Indices for each amino acid according to their position in hmm files
	static const int kMaxLineLength_;
};
vector<vector<double> > CalculateErrorScore(string seq, const bool mode); // The index begins with 1 to be consistent with Viterbi
int Codon2AminoAcidIndex(const char* codon);
char Codon2AminoAcid(const char* codon);

#endif
  