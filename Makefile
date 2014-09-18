all:
	g++ -O3  -w -o hmmscore hmmscore.cpp viterbi_model.cpp DNA_seq.cpp hmm_model.cpp
	g++ -O3  -w -o generate_domain_expression generate_domain_expression.cpp
	bash -c "bash chmod_scripts.sh"
