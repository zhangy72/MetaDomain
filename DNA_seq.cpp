// Copyright (C) 2011 Yuan Zhang, Yanni Sun.
// You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
// implementation of DNASeq class.

#include "DNA_seq.h"
using namespace std;
const int DNASeq::kMaxLineLength_ = 2048;
const int DNASeq::kCodonNum_ = 64;
const char* DNASeq::kCodonTable_ = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF";
const int DNASeq::kCodonIndexTable_[65] = {8,11,8,11,16,16,16,16,14,15,14,15,7,7,10,7,13,6,13,6,12,12,12,12,14,14,14,14,9,9,9,9,3,2,3,2,0,0,0,0,5,5,5,5,17,17,17,17,20,19,20,19,15,15,15,15,20,1,18,1,9,4,9,4,20};
DNASeq::DNASeq(const string& seq_name, const string& seq)
:seq_(seq),
 seq_name_(seq_name)
{
	length_ = seq_.length();
  for (int i = 0; i < length_; ++i) seq_[i] = toupper(seq_[i]);
}

string DNASeq::ReverseComplement() const {
  string reverse_complement = seq_;
	for (int i = 0; i < length_; i++) {
		switch (seq_[i]) {
			case 'A':
				reverse_complement[length_ - 1 - i]= 'T';
				break;
			case 'C':
				reverse_complement[length_ - 1 - i] = 'G';
				break;
			case 'G':
				reverse_complement[length_ - 1 - i] = 'C';
				break;
			case 'T':
				reverse_complement[length_ - 1 - i] = 'A';
				break;
			default:
				reverse_complement[length_ - 1 - i] = seq_[i];
				break;
		}
	}
  return reverse_complement;
}

//6-frame translation
vector<string> DNASeq::protein_seqs() {
	vector<string> protein_seqs(6);
	const char* seq = seq_.c_str();
	const char* reverse_complement = ReverseComplement().c_str();
	for (int i = 0; i < 6; i ++) {
		if (i < 3) {
			protein_seqs[i].reserve((length_ - i) / 3);
			for (int j = i; j + 3 <= length_; j += 3) {
				protein_seqs[i].push_back(Codon2AminoAcid(seq + j));
			}
		} else {
			int k = i - 3;
			protein_seqs[i].reserve((length_ - k) / 3);
			for (int j = k; j + 3 <= length_; j += 3) {
				protein_seqs[i].push_back(Codon2AminoAcid(reverse_complement + j));
			}
		}
	}
	return protein_seqs;
}
vector<vector<double> > CalculateErrorScore(string seq, const bool mode) { 
	// In this function all probabilities are converted to scores. Probabilities are hidden from users
	double error_score[7];
	const int length = seq.length();
	vector<vector<double> > seq_error_score(2, vector<double>(length + 1)); // 0 position is vacant. 
	if ( mode == false) { // mode false means we use the reference error model
	  const double kNonHomo = log(0.0007f);
	  const double kHomo = log(0.0044f);
	  error_score[0] = kNonHomo;
	  for (int i = 1; i < 7; i++) {
	    error_score[i] = kHomo;
	  }
	} else {
	  error_score[0] = log(0.000532f);
	  error_score[1] = log(0.000698f);
	  error_score[2] = log(0.00102f);
	  error_score[3] = log(0.000688f);
	  error_score[4] = log(0.0372f);
	  error_score[5] = log(0.00167f);
	  error_score[6] = log(0.143f);
	}
	int* occurence = new int[length]; // This begins with 0 to be consistent with seq_
	for(int i = 0; i < length;) {
		if (i == length - 1) {
			occurence[i] = 1;
			break;
		}
		int count = 1; // Count for each repeated region
		int j;
		for (j = i; j < length - 1; j++) {
			if(seq[j] == seq[j+1]){
				count++;
				if (j == length - 2) {
					for (int k = i; k <= j + 1; k++) {
						occurence[k] = count;
					}
					break;
				}
			} else {
				for (int k = i; k < j + 1; k++) {
					occurence[k] = count;
				}
				break;
			}
		}
		i = j + 1;
	}
	// Note that occurence begins with 0 and seq_error_score begins with 1
	for (int i = 0; i < length; i++) {
		if (i == 0){
			seq_error_score[0][i + 1] = error_score[min(occurence[i + 1] - 1, 6)];
			seq_error_score[1][i + 1] = error_score[min(occurence[i + 1] - 1, 6)];
		}
		else if (i == length - 1) {
			seq_error_score[0][i + 1] = error_score[min(occurence[i-1] - 1, 6)];
			seq_error_score[1][i + 1] = error_score[min(occurence[i] - 1, 6)];
		}
		else{
			seq_error_score[0][i + 1] = error_score[min(max(max(occurence[i - 1],occurence[i]), occurence[i + 1]) - 1, 6)];
			seq_error_score[1][i + 1] = error_score[min(max(occurence[i],occurence[i + 1]) - 1, 6)];
		}
	}
	delete[] occurence;
	return seq_error_score;
}
int Codon2AminoAcidIndex(const char* codon) { // Convert a codon into its corresponding amino acid character
		int index = 0;
		for (int i =0; i < 3; i++) {
			switch (codon[2 - i]) {
				case 'A':
					//index += pow(4.,i) * 0;
					break;
				case 'C':
					index += (1<<(i<<1));
					break;
				case 'G':
					//index += pow(4.,i) * 2;
					index |= (1<<(i<<1))<<1;
					break;
				case 'T': 
					//index += pow(4.,i) * 3;
					index |= (1<<((i<<1)|1))+(1<<(i<<1));
					break;
				default:
					return 20; // Invalid characters lead to X amino acid
			}
		}
		return DNASeq::kCodonIndexTable_[index];
}
char Codon2AminoAcid(const char* codon) { // Convert a codon into its corresponding amino acid character
		double index = 0;
		for (int i =0; i < 3; i++) {
			switch (codon[2 - i]) {
				case 'A':
					index += pow(4.,i) * 0;
					break;
				case 'C':
					index += pow(4.,i) * 1;
					break;
				case 'G':
					index += pow(4.,i) * 2;
					break;
				case 'T': 
					index += pow(4.,i) * 3;
					break;
				default:
					return 'X'; // Invalid characters lead to X amino acid
			}
		}
		return DNASeq::kCodonTable_[static_cast<int>(index)];
	}
