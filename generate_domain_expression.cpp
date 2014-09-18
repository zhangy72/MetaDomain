// Copyright (C) 2011 Yuan Zhang, Yanni Sun.
// You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
// classifies input protein domain based on alignment result by hmmscore.

#include <cstdio>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <string>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
using namespace std;
enum GeneExpression {kUntranscribed, kUnknown, kLowTranscribed, kMediumTranscribed, kHighTranscribed};
class AlignedTag { // Note that in this class all the positions are in the scale of amino acid residues
public:
  int length_;
	int begin_state_;
	int end_state_;
	int alignment_length_;
	int weight_;
	double score_;
	double score_threshold_;
	bool is_duplicate_;
	string name_;
	string domain_;
	AlignedTag(const string& name);
};
class AlignedDomain {
public:
	int length_;
	int tag_num_;
	int read_num_;
	double average_depth_;
	double coverage_;
	bool transcribed_;
	string name_;
	vector<int> depths_; // depth for all positions. 0 position is vacant
	vector<AlignedTag> aligned_tags_;
	AlignedDomain(const string& name);
	void set_tag_num();
	void set_read_num();
	void set_average_depth();
	void set_coverage();
	double coverage_threshold(const int& read_num_threshold, const int& read_length) const;
};
void TokenizeString(const char* line, vector<string>& tokens, const string& delimiters);
void GetDomainLengths(const char* file_name, map<string, int>& domain_lengths);
void GetTags(const char* file_name, vector<AlignedTag>& tags);
void ConstructDomains(const map<string, int>& domain_lengths, map<string, AlignedDomain>& domains);
void GetDomains(const vector<AlignedTag>& tags, map<string, AlignedDomain>& domains);
void LabelDomains(const int& read_num_threshold, const double& average_depth_threshold, const double& coverage_depth,
	              map<string, AlignedDomain>& domains);
void RemoveDuplicateTags(map<string, AlignedDomain>& domains);
// set basic parameters except for the domain name and length
void SetDomainParameters(map<string, AlignedDomain>& domains);

// the format of the file:
// domain_name length aligned_tags aligned_reads coverage average_depth
void OutputDomains(const map<string, AlignedDomain>& domains, 
                   const char* out_file_name); 

int main(int argc, char* argv[]) {
	if (argc != 6) {
    printf("Parameters: <domain length file> <hmmscore file>");
    printf(" <read number threshold> <average depth threshold>"); 
    printf(" <coverage threshold> <output file name>\n");
		exit(1);
	}
	const char* length_file_name = argv[1]; // note that the length here is in residue not bases
	const char* hmmscore_file_name = argv[2];  // this file contains hmmscore results for all domains.
	const int read_num_threshold = atoi(argv[3]); 
	const double coverage_threshold = atof(argv[4]);
	const char* out_file_name = argv[5];  // the output file name.
	double average_depth_threshold = 0;
	map<string, int> domain_lengths;
	vector<AlignedTag> tags;
	map<string, AlignedDomain> domains;
	GetDomainLengths(length_file_name, domain_lengths);
	GetTags(hmmscore_file_name, tags);
	ConstructDomains(domain_lengths, domains);
	GetDomains(tags, domains);
	RemoveDuplicateTags(domains);
	SetDomainParameters(domains);
	LabelDomains(read_num_threshold, average_depth_threshold, coverage_threshold, domains);
	OutputDomains(domains, out_file_name);
	return 0;
}
AlignedTag::AlignedTag(const string& name) {
	name_ = name;
	begin_state_ = 0;
	end_state_ = 0;
	alignment_length_ = 0;
	weight_ = 1;
	score_ = 0.0;
	score_threshold_ = 0.0;
	is_duplicate_ = false;
}
AlignedDomain::AlignedDomain(const string& name) {
	name_ = name;
	length_ = 0;
	tag_num_ = 0;
	read_num_ = 0;
	average_depth_ = 0.0;
	coverage_ = 0.0;
	transcribed_ = false;
}
double AlignedDomain::coverage_threshold(const int& read_num_threshold, const int& read_length) const {
	return (static_cast<double>(read_num_threshold * read_length + read_length) / 2 / length_);
}
void AlignedDomain::set_tag_num() {
	tag_num_ = aligned_tags_.size();
}
void AlignedDomain::set_read_num() {
	for (int i = 0; i < aligned_tags_.size(); i++) {
		read_num_ += aligned_tags_[i].weight_;
	}
}
void AlignedDomain::set_average_depth() {
	depths_ = vector<int>(length_ + 1, 0); // the 0 position is not used
	for (int i = 0; i < aligned_tags_.size(); i++) {
		int first_position = max(1, aligned_tags_[i].begin_state_);
		int last_position = min(length_, aligned_tags_[i].end_state_);
		int tag_weight = aligned_tags_[i].weight_;
		for (int j = first_position; j <= last_position; j++) {
			depths_[j] += tag_weight;
		}
	}
	for (int i = 1; i < depths_.size(); i++) {
		average_depth_ += depths_[i];
	}
	average_depth_ = average_depth_ / length_;
}
void AlignedDomain::set_coverage() {
	for (int i = 1; i < depths_.size(); i++) {
		if (depths_[i] > 0) {
			coverage_++;
		}
	}
	coverage_ = coverage_ / length_;
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
void GetDomainLengths(const char* file_name, map<string, int>& domain_lengths) {
	ifstream input_file(file_name);
	if (!input_file.is_open()) {
		cerr << "Cannot open domain label file!" << endl;
		exit(1);
	}
	const int kMaxLineLength = 4096;
	char line[kMaxLineLength];
	while (input_file.getline(line, kMaxLineLength)) {
		if (strlen(line) == 0) {
			break;
		}
		vector<string> tokens;
		TokenizeString(line, tokens, " ");
		domain_lengths.insert(pair<string, int>(tokens[0].substr(0, 7), atoi(tokens[1].c_str())));
	}
  input_file.close();
}
void GetTags(const char* file_name, vector<AlignedTag>& tags) {
	ifstream input_file(file_name);
	if (!input_file.is_open()) {
		cerr << "Cannot open domain label file!" << endl;
		exit(1);
	}
	const int kMaxLineLength = 4096;
	char line[kMaxLineLength];
	while (input_file.getline(line, kMaxLineLength)) {
		if (strlen(line) == 0) {
			break;
		}
		vector<string> tokens;
		TokenizeString(line, tokens, " ");
		AlignedTag tag(tokens[0].substr(1));
		tag.domain_ = tokens[1].substr(tokens[1].find("=") + 1, 7); 
		tag.score_ = atof(tokens[2].substr(tokens[2].find("=") + 1).c_str());
		tag.score_threshold_ = atof(tokens[3].substr(tokens[3].find("=") + 1).c_str());
		tag.begin_state_ = atoi(tokens[4].substr(tokens[4].find("=") + 1).c_str());
		tag.end_state_ = atoi(tokens[5].substr(tokens[5].find("=") + 1).c_str());
		tag.alignment_length_ = tag.end_state_ - tag.begin_state_ + 1;
		tag.length_ = atoi(tokens[6].substr(tokens[6].find("=") + 1).c_str());
		tags.push_back(tag);
	}
	input_file.close();
}
// this function only initializes domain's name and length
void ConstructDomains(const map<string, int>& domain_lengths, map<string, AlignedDomain>& domains) {
	map<string, int>::const_iterator it;
	for (it = domain_lengths.begin(); it != domain_lengths.end(); it++) {
		AlignedDomain domain(it->first);
		domain.length_ = it->second;
		domains.insert(pair<string, AlignedDomain>(it->first, domain));
	}
}
// this function purely adds tags into the domain without setting its parameters
void GetDomains(const vector<AlignedTag>& tags, map<string, AlignedDomain>& domains) {
	map<string, AlignedDomain>::iterator it;
	for (int i = 0; i < tags.size(); i++) {
		it = domains.find(tags[i].domain_);
		assert(it != domains.end());
		it->second.aligned_tags_.push_back(tags[i]);
	}
}
void LabelDomains(const int& read_num_threshold, const double& average_depth_threshold, const double& coverage_threshold,
	              map<string, AlignedDomain>& domains) {
	map<string, AlignedDomain>::iterator it;
	double new_coverage_threshold = 0.0; // this coverage threshold is domain independent. It takse the average coverage
	const int kReadLength = 13;
	for (it = domains.begin(); it != domains.end(); it++) {
		new_coverage_threshold = it->second.coverage_threshold(read_num_threshold, kReadLength);
		if (it->second.read_num_ >= read_num_threshold && it->second.average_depth_ 
			>= average_depth_threshold && it->second.coverage_ >= coverage_threshold) { // Now use domain independent threshold
				it->second.transcribed_ = true;
		}
	}
}
// For tags with the same alignment beginning state only leave the one with largest weight
void RemoveDuplicateTags(map<string, AlignedDomain>& domains) {
	map<string, AlignedDomain>::iterator it;
	for (it = domains.begin(); it != domains.end(); it++) {
		int tag_num =  it->second.aligned_tags_.size();
		if (tag_num == 0) {
			continue;
		}
		map<int, pair<int, int> > aligned_positions;
		pair<map<int, pair<int, int> >::iterator, bool> ret;
		for (int i = 0; i <tag_num; i++) {
			int begin_state = it->second.aligned_tags_[i].begin_state_;
			int end_state = it->second.aligned_tags_[i].end_state_;
			int tag_index = i;
			ret = aligned_positions.insert(pair<int, pair<int, int> >(begin_state, pair<int, int>(tag_index, end_state)));
			if (ret.second == false) {
				if (end_state > ret.first->second.second) {
					it->second.aligned_tags_[ret.first->second.first].is_duplicate_ = true; // the original tag is duplicate
					ret.first->second.second = end_state; // update the ending state
					ret.first->second.first = tag_index; // replace the current tag index to i
				} else {
					it->second.aligned_tags_[tag_index].is_duplicate_ = true;
				}
			}
		}
		for (int i = tag_num - 1; i >=0; i--) {
			if (it->second.aligned_tags_[i].is_duplicate_ == true) {
				it->second.aligned_tags_.erase(it->second.aligned_tags_.begin() + i);
			}
		}
	}
}
void SetDomainParameters(map<string, AlignedDomain>& domains) {
	map<string, AlignedDomain>::iterator it;
	for (it = domains.begin(); it != domains.end(); it++) {
		it->second.set_tag_num();
		it->second.set_read_num();
		it->second.set_average_depth();
		it->second.set_coverage();
	}
}
void OutputDomains(const map<string, AlignedDomain>& domains, 
                   const char* out_file_name) {
  ofstream out_file(out_file_name);
  assert(out_file.is_open());
  map<string, AlignedDomain>::const_iterator domain_it;
  string status;  // whether the domain is transcribed or not
  for (domain_it = domains.begin(); domain_it != domains.end(); ++domain_it) {
    if (domain_it->second.transcribed_) {
      status = "transcribed";
    } else {
      status = "nontranscribed";
    }
    out_file << domain_it->first << " " << domain_it->second.length_ << " " <<  
				domain_it->second.read_num_ << " " << domain_it->second.coverage_ << " " << status << endl;
  }
  out_file.close();
}
