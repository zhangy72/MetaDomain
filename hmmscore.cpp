// Copyright (C) 2011 Yuan Zhang, Yanni Sun.
// You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
// HMMSCORE is used to align short reads with a Pfam HMM file and generate significant hits. 

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "DNA_seq.h"
#include "hmm_model.h"
#include "viterbi_model.h"
using namespace std;
int main(int argc, char* argv[])
{
  if (argc != 5) {
    cout << "Usage: hmmscore <hmm file> <fasta_file> "
         << "<threshold rate> <0=no reverse complement; 1=consider both>" << endl;
    return 1; 	
  }
  const char* hmm_file_name = argv[1]; 
  const char* seq_file_name = argv[2];
  const double kRate =atof(argv[3]);
  bool consider_rc = static_cast<bool>(atoi(argv[4]));
  const int kLengthThreshold = 0;
  ifstream seq_file(seq_file_name);
  if (!seq_file) {
    cerr << "Cannot open fasta file!" << endl;
    seq_file.close();
    return 1;
  }
  ifstream hmm_file(hmm_file_name);
  if (!hmm_file) {
    cerr << "Cannot open HMM file!" << endl;
    hmm_file.close();
    return 1;
  }
  HmmModel hmm;
  ConstructHmmModel(hmm_file, hmm);
  hmm_file.close();
  //Begin the Viterbi process. Each time store one sequence and process it. 
  string line;
  string seq;
  string seq_name;
  while (getline(seq_file, line)) {
    if (line.find_first_not_of(" \t") == string::npos) {
      continue;
    }
    if (line[0] == '>') {
      if (!seq_name.empty() && !seq.empty()) {
        DNASeq dna(seq_name, seq);
        Alignment alignment(false, kLengthThreshold, kRate);
        Viterbi(hmm, dna, alignment);
        if (alignment.score() >= alignment.score_threshold()) {
          cout << ">" << alignment.seq_name() << " hmm_acc=" << alignment.hmm_acc() 
               << " score=" << alignment.score() << " score_threshold=" 
               << alignment.score_threshold() << " begin_state=" << alignment.begin_state() 
	       << " end_state=" << alignment.end_state() << " alignment_length=" 
               << alignment.align_length() << " read_aa_length_=" 
               << alignment.DNA_length() / 3 << endl;	
        } else {
         if (consider_rc) { 
           alignment.set_reverse_complement(true);  // use reverse complement of the read
           Viterbi(hmm, dna, alignment);			
           if (alignment.score() >= alignment.score_threshold()) {
             cout << ">" << alignment.seq_name() << " hmm_acc=" << alignment.hmm_acc() 
                  << " score=" << alignment.score() << " score_threshold=" 
                  << alignment.score_threshold() << " begin_state=" 
                  << alignment.begin_state() << " end_state=" << alignment.end_state()
                  << " alignment_length=" << alignment.align_length()
                  << " read_aa_length_=" << alignment.DNA_length() / 3 << endl;
             }
           }
        }
        seq.clear();
      }
      seq_name = line.substr(1); // seq_name will not include the prefix ">"
    } else {
      seq.append(line);
    }
  }
  seq_file.close();

  // handle the last sequence.
  DNASeq dna(seq_name, seq);
  Alignment alignment(false, kLengthThreshold, kRate);
  Viterbi(hmm, dna, alignment);
  if (alignment.score() >= alignment.score_threshold()) {
    cout << ">" << alignment.seq_name() << " hmm_acc=" << alignment.hmm_acc() 
         << " score=" << alignment.score() << " score_threshold=" 
         << alignment.score_threshold() << " begin_state=" << alignment.begin_state() 
         << " end_state=" << alignment.end_state() << " alignment_length=" 
         << alignment.align_length() << " read_aa_length_=" 
         << alignment.DNA_length() / 3 << endl;	
  } else {
    if (consider_rc) { 
      alignment.set_reverse_complement(true);  // use reverse complement of the read
      Viterbi(hmm, dna, alignment);			
      if (alignment.score() >= alignment.score_threshold()) {
        cout << ">" << alignment.seq_name() << " hmm_acc=" << alignment.hmm_acc() 
             << " score=" << alignment.score() << " score_threshold=" 
             << alignment.score_threshold() << " begin_state=" 
             << alignment.begin_state() << " end_state=" << alignment.end_state()
             << " alignment_length=" << alignment.align_length()
             << " read_aa_length_=" << alignment.DNA_length() / 3 << endl;
      }
    }
  }
  return 0;
}

