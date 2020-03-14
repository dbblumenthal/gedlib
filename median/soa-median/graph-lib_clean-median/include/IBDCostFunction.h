#ifndef __IBD_COSTFUNCTION_H__
#define __IBD_COSTFUNCTION_H__

#include "GraphEditDistance.h"
#include "IBDGraph.h"

void
tokenizeB(const std::string & sentence, char sep, std::vector<std::string> & words) {
	bool outside_quotes{true};
	std::size_t word_length{0};
	std::size_t pos_word_start{0};
	for (std::size_t pos{0}; pos < sentence.size(); pos++) {
		if (sentence.at(pos) == '\'') {
			if (not outside_quotes and pos < sentence.size() - 1) {
				if (sentence.at(pos + 1) != sep) {
					throw std::runtime_error("Sentence contains closing single quote which is followed by a char different from " + std::to_string(sep) + ".");
				}
			}
			word_length++;
			outside_quotes = not outside_quotes;
		}
		else if (outside_quotes and sentence.at(pos) == sep) {
			if (word_length > 0) {
				words.push_back(sentence.substr(pos_word_start, word_length));
			}
			pos_word_start = pos + 1;
			word_length = 0;
		}
		else {
			word_length++;
		}
	}
	if (not outside_quotes) {
		throw std::runtime_error("Sentence contains unbalanced single quotes.");
	}
	if (word_length > 0) {
		words.push_back(sentence.substr(pos_word_start, word_length));
	}
}

class IBDDistanceCost:public EditDistanceCost<int,double>
{

private:

  double _alpha;
  std::vector<std::vector<double>> node_rel_costs_;
  std::vector<std::size_t> otus_;
  
public:

  double NodeSubstitutionCost(GNode<int,double> * n1,GNode<int,double> * n2,
			      Graph<int,double> * g1,Graph<int,double> * g2)
  {
    return _alpha * node_rel_costs_.at(n1->attr).at(n2->attr);
  }
  double NodeDeletionCost(GNode<int,double> * n1,Graph<int,double> * g1)
  {
    return _alpha;
  }
  double NodeInsertionCost(GNode<int,double> * n2,Graph<int,double> * g2)
  {
    return _alpha;
  }
  double EdgeSubstitutionCost(GEdge<double> * e1,GEdge<double> * e2,
			      Graph<int,double> * g1,Graph<int,double> * g2)
  {
    return (1 - _alpha) * std::fabs(e1->attr - e2->attr);
  }
  double EdgeDeletionCost(GEdge<double> * e1,Graph<int,double> * g1)
  {
    return 1.0-_alpha;
  }
  double EdgeInsertionCost(GEdge<double> * e2,Graph<int,double> * g2)
  {
    return 1.0-_alpha;
  }
  double SubstitutionCost(int label_1, int label_2) // for nodes
  {
    return _alpha * node_rel_costs_.at(label_1).at(label_2);
  }
  double EdgeSubstitutionCost(double label_1, double label_2)
  {
    return (1 - _alpha) * std::fabs(label_1 - label_2);
  }

  IBDDistanceCost * clone() const {return new IBDDistanceCost(*this);}

  IBDDistanceCost(const std::string & otu_distances, double alpha = 0.5) :
  _alpha(alpha), node_rel_costs_(), otus_()
  {
	std::ifstream csv_file(otu_distances);
	std::string row;
	std::getline(csv_file, row);
	std::vector<std::string> row_as_vector;
	tokenizeB(row, ',', row_as_vector);
	std::size_t max_otu{0};
	std::vector<std::string> otu_as_vector;
	for (const std::string & otu_with_prefix : row_as_vector) {
		otu_as_vector.clear();
		tokenizeB(otu_with_prefix, '_', otu_as_vector);
		otus_.emplace_back(std::stoul(otu_as_vector.at(1)));
		max_otu = std::max(max_otu, otus_.back());
	}
	node_rel_costs_ = std::vector<std::vector<double>>(max_otu + 1, std::vector<double>(max_otu + 1, 0));

	std::size_t otu_1;
	std::size_t otu_2;
	while(std::getline(csv_file, row)) {
		row_as_vector.clear();
		tokenizeB(row, ',', row_as_vector);
		otu_as_vector.clear();
		tokenizeB(row_as_vector.at(0), '_', otu_as_vector);
		otu_1 = std::stoul(otu_as_vector.at(1));
		for (std::size_t pos{1}; pos < row_as_vector.size(); pos++) {
			otu_2 = otus_.at(pos - 1);
			node_rel_costs_.at(otu_1).at(otu_2) = std::stod(row_as_vector.at(pos));
		}
	}
	double max_node_rel_cost(0);
	for (std::size_t otu_1(0); otu_1 <= max_otu; otu_1++) {
		for (std::size_t otu_2(0); otu_2 <= max_otu; otu_2++) {
			max_node_rel_cost = std::max(max_node_rel_cost, node_rel_costs_.at(otu_1).at(otu_2));
		}
	}
	for (std::size_t otu_1(0); otu_1 <= max_otu; otu_1++) {
		for (std::size_t otu_2(0); otu_2 <= max_otu; otu_2++) {
			node_rel_costs_.at(otu_1).at(otu_2) /= max_node_rel_cost;
		}
	}
}


};


#endif
