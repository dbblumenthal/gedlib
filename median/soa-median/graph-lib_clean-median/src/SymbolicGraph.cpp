/*
 * @file SymbolicGraph.cpp
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>>
 * @version     0.0.1 - Sat Feb  4 2017
 *
 */
#include <string>
#include <map>

#include "SymbolicGraph.h"
using namespace std;

/*
 * Map to convert atom symbol to int
 */
static std::map<std::string, int> AtomTable;

void fillAtomTable(std::map<std::string, int> & AtomTable){
	AtomTable["H"] = 1;
	AtomTable["He"] = 2;
	AtomTable["Li"] = 3;
	AtomTable["Be"] = 4;
	AtomTable["B"] = 5;
	AtomTable["C"] = 6;
	AtomTable["N"] = 7;
	AtomTable["O"] = 8;
	AtomTable["F"] = 9;
	AtomTable["Ne"] = 10;
	AtomTable["Na"] = 11;
	AtomTable["Mg"] = 12;
	AtomTable["Al"] = 13;
	AtomTable["Si"] = 14;
	AtomTable["P"] = 15;
	AtomTable["S"] = 16;
	AtomTable["Cl"] = 17;
	AtomTable["Ar"] = 18;
	AtomTable["K"] = 19;
	AtomTable["Ca"] = 20;
	AtomTable["Sc"] = 21;
	AtomTable["Ti"] = 22;
	AtomTable["V"] = 23;
	AtomTable["Cr"] = 24;
	AtomTable["Mn"] = 25;
	AtomTable["Fe"] = 26;
	AtomTable["Co"] = 27;
	AtomTable["Ni"] = 28;
	AtomTable["Cu"] = 29;
	AtomTable["Zn"] = 30;
	AtomTable["Ga"] = 31;
	AtomTable["Ge"] = 32;
	AtomTable["As"] = 33;
	AtomTable["Se"] = 34;
	AtomTable["Br"] = 35;
	AtomTable["Kr"] = 36;
	AtomTable["Rb"] = 37;
	AtomTable["Sr"] = 38;
	AtomTable["Y"] = 39;
	AtomTable["Zr"] = 40;
	AtomTable["Nb"] = 41;
	AtomTable["Mo"] = 42;
	AtomTable["Tc"] = 43;
	AtomTable["Ru"] = 44;
	AtomTable["Rh"] = 45;
	AtomTable["Pd"] = 46;
	AtomTable["Ag"] = 47;
	AtomTable["Cd"] = 48;
	AtomTable["In"] = 49;
	AtomTable["Sn"] = 50;
	AtomTable["Sb"] = 51;
	AtomTable["Te"] = 52;
	AtomTable["I"] = 53;
	AtomTable["Xe"] = 54;
	AtomTable["Cs"] = 55;
	AtomTable["Ba"] = 56;
	AtomTable["La"] = 57;
	AtomTable["Ce"] = 58;
	AtomTable["Pr"] = 59;
	AtomTable["Nd"] = 60;
	AtomTable["Pm"] = 61;
	AtomTable["Zr"] = 62;
	AtomTable["Eu"] = 63;
	AtomTable["Gd"] = 64;
	AtomTable["Tb"] = 65;
	AtomTable["Dy"] = 66;
	AtomTable["Ho"] = 67;
	AtomTable["Er"] = 68;
	AtomTable["Tm"] = 69;
	AtomTable["Yb"] = 70;
	AtomTable["Lu"] = 71;
	AtomTable["Hf"] = 72;
	AtomTable["Ta"] = 73;
	AtomTable["W"] = 74;
	AtomTable["Re"] = 75;
	AtomTable["Os"] = 76;
	AtomTable["Ir"] = 77;
	AtomTable["Pt"] = 78;
	AtomTable["Au"] = 79;
	AtomTable["Hg"] = 80;
	AtomTable["Tl"] = 81;
	AtomTable["Pb"] = 82;
	AtomTable["Bi"] = 83;
	AtomTable["Po"] = 84;
	AtomTable["At"] = 85;
	AtomTable["Rn"] = 86;
	AtomTable["Fr"] = 87;
	AtomTable["Ra"] = 88;
	AtomTable["Ac"] = 89;
	AtomTable["Th"] = 90;
	AtomTable["Pa"] = 91;
	AtomTable["U"] = 92;
	AtomTable["Np"] = 93;
	AtomTable["Pu"] = 94;
	AtomTable["Am"] = 95;
	AtomTable["Cm"] = 96;
	AtomTable["Bk"] = 97;
	AtomTable["Cf"] = 98;
	AtomTable["Es"] = 99;
	AtomTable["Fm"] = 100;
	AtomTable["Md"] = 101;
	AtomTable["No"] = 102;
	AtomTable["Lr"] = 103;
	AtomTable["Rf"] = 104;
	AtomTable["Db"] = 105;
	AtomTable["Sg"] = 106;
	AtomTable["Bh"] = 107;
	AtomTable["Hs"] = 108;
	AtomTable["Mt"] = 109;
	AtomTable["Ds"] = 110;
	AtomTable["Rg"] = 111;
	AtomTable["Cn"] = 112;
	AtomTable["Uut"] = 113;
	AtomTable["Uuq"] = 114;
	AtomTable["Uup"] = 115;
	AtomTable["Uuh"] = 116;
	AtomTable["Uus"] = 117;
	AtomTable["Uuo"] = 118;
	AtomTable["D"] = 119; // Deuterium (isotope de H)
}


int SymbolicGraph::readChemicalEdgeLabel(TiXmlElement *elem){
	int bond_type = -1;
	TiXmlElement* child = elem->FirstChildElement();
	while ( child ){
		std::string childName =  child->Attribute("name");
		TiXmlElement * child2 = child->FirstChildElement();
		if ( child2 ){
			if(childName.compare("valence")==0){
				bond_type = std::stoi(child2->GetText());
			}
		}
		child = child->NextSiblingElement(); // iteration
	}
	return bond_type;
}

int SymbolicGraph::readChemicalNodeLabel(TiXmlElement *elem){
	int atom = -1;
	TiXmlElement* child = elem->FirstChildElement();
	while ( child ){
		std::string childName =  child->Attribute("name");
		TiXmlElement * child2 = child->FirstChildElement();
		if ( child2 ){
			if(childName.compare("chem")==0){
				try {
					atom = std::stoi(child2->GetText());
				}
				catch (const std::invalid_argument& ia) {
					atom = AtomTable[child2->GetText()];
				}
			}
		}
		child = child->NextSiblingElement(); // iteration au prochain attribut
	}
	return atom;
}



int SymbolicGraph::readGraphmlEdgeLabel(TiXmlElement *elem){
	return 1;
}

int SymbolicGraph::readGraphmlNodeLabel(TiXmlElement *elem){
	int atom = -1;
	TiXmlElement* child = elem->FirstChildElement();
	if ( child ){
		std::string label =  child->GetText();
		atom = atoi(label.c_str());
	}
	return atom;
}




SymbolicGraph::SymbolicGraph(const char * filename):Graph<int,int>(false){
	fillAtomTable(AtomTable);
	const char * ext = strrchr(filename,'.');
	if (strcmp(ext,".ct") == 0){

		std::ifstream file(filename,std::ios::in);
		if (file.is_open())
		{
			char * s = new char[255];
			file.getline(s, 255); // The first line is useless
			int  mynbNodes, mynbEdges;
			file >> mynbNodes;
			file >> mynbEdges;

			float coord;
			std::string index;
			for(int i=0; i<mynbNodes; i++)
			{
				// ignore x y and z coordinates
				file >> coord; file >> coord; file >> coord;
				file >> index;

				Add(new GNode<int,int>(i,AtomTable[index]));
				// int X=std::rand()%10+1;// à commenter pour avoir les vrais labels sur les noeuds (et decommenter celle au dessus) => v9
				//  if (X<=8) X=1; else X=std::rand()%5+1; // à commenter pour v9, à décommenter pout v10
				// Add(new GNode<int,int>(i,X));// à commenter pour avoir les vrais labels sur les noeuds
			}

			// Creation of the edges
			int start, end, label;
			for (int i=0; i<mynbEdges; i++)
			{
				file >> start; file >> end; file >> label;
				//int X=std::rand()%5+1;// à commenter pour avoir les vrais labels sur les arêtes => v9
				//  if (X<=4) X=1; else X=std::rand()%2+1; // à commenter pour v9, à décommenter pout v10
				// label = X; // à commenter pour avoir les vrais labels sur les aretes
				Link(start-1, end-1, label);
				file >> start; //ignore last value
				//file.ignore(255, '\n'); // go to the next line
			}
			delete [] s;
		}
	}else if (strcmp(ext,".gxl") == 0){
		GraphLoadGXL(filename,readChemicalNodeLabel,readChemicalEdgeLabel);
	}else if (strcmp(ext, ".graphml") == 0){
		GraphLoadGXL(filename,readGraphmlNodeLabel,readGraphmlEdgeLabel);
	}else{
		std::cerr << "Unsupported file format. filename = "<<filename << std::endl;
	}
}

SymbolicGraph::SymbolicGraph(int * am, int nb_nodes, bool directed):Graph<int,int>(directed){
	for(int n = 0; n<nb_nodes; n++){ //We traverse diagonal for node labels
		Add(new GNode<int,int>(n,am[sub2ind(n,n, nb_nodes)]));
	}

	for(int n = 0; n<nb_nodes; n++){
		int start= (directed)?0:n+1;
		for(int m = start; m<nb_nodes; m++) {
			int edge = am[n+nb_nodes*m];
			if(edge > 0){
				Link(n,m,edge);
			}
		}
	}
}

int * SymbolicGraph::getLabeledAdjacencyMatrix(){
	int * am=new int[Size()*Size()];
	memset(am,0,sizeof(int)*Size()*Size());
	map<int,int> conv_nodes;
	int index = 1;
	for(int i=0;i<Size();i++)
		if(conv_nodes.find((*this)[i]->attr) == conv_nodes.end())
			conv_nodes[(*this)[i]->attr] = index++;

	for(int i=0;i<Size();i++){
		am[sub2ind(i,i,Size())] = conv_nodes[(*this)[i]->attr];
		GEdge<int> * e = (*this)[i]->getIncidentEdges();
		while(e){
			int j = e->IncidentNode();
			am[sub2ind(i,j,Size())] = e->attr;
			e=e->Next();
		}
	}
	return am;
}


bool writeCTfile(Graph<int,int>& graph, std::ofstream& output){
	try{
		output << "Generated graph" << std::endl;
		output << graph.Size() << " " << graph.getNbEdges() << std::endl;
		for (int i=0; i<graph.Size(); i++){
			int label = graph[i]->attr;
			std::map<std::string, int>::iterator it = AtomTable.begin();
			while (it != AtomTable.end() && it->second != label) it++;
			output << "   0.0    0.0    0.0   ";
			if (it != AtomTable.end())
				output << it->second << std::endl;
			else
				output << "?" << std::endl;
		}

		for (int i=0; i<graph.Size(); i++){
			for (int j=0; j<i; j++){
				if (graph.isLinked(i,j))
					output << i+1 << " " << j+1 << "  " << graph.getEdge(i,j)->attr << "  0" << std::endl;
			}
		}
	}
	catch(std::exception & e){
		std::cerr << "[E]  Error while writing file : " << e.what() << std::endl;
		return false;
	}
	return true;
}

bool writeCTfile_std(Graph<int,int>& graph){
	try{
		std::cout << graph.Size() << " " << graph.getNbEdges() << std::endl;
		for (int i=0; i<graph.Size(); i++){
			int label = graph[i]->attr;
			//std::cout << "label sommet " << i << " = " << label << std::endl;
			std::map<std::string, int>::iterator it = AtomTable.begin();
			while (it != AtomTable.end() && it->second != label) it++;
			std::cout << "   0.0    0.0    0.0   ";
			if (it != AtomTable.end())
				std::cout << it->first << std::endl;
			else
				std::cout << "?" << std::endl;
		}

		for (int i=0; i<graph.Size(); i++){
			for (int j=0; j<i; j++){
				if (graph.isLinked(i,j))
					std::cout << i+1 << " " << j+1 << "  " << graph.getEdge(i,j)->attr << "  " << graph.getEdge(i,j)->attr  << std::endl;
			}
		}
	}
	catch(std::exception & e){
		std::cerr << "[E]  Error while writing file : " << e.what() << std::endl;
		return false;
	}
	return true;
}


