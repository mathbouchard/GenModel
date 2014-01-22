/*
 * HyperGraph.h
 *
 *  Created on: 2011-03-17
 *      Author: mbouchard
 */

#ifndef HYPERGRAPH_H_
#define HYPERGRAPH_H_

#include <vector>
#include <limits>
#include <string>
#include <stdlib.h>
#define uint unsigned int

using namespace std;

class HyperGen;
class GenFunction;

class HyperGraph
{
public:
	HyperGraph() {hn=0; he = 0;}
	~HyperGraph() {}
	uint hn;
	uint he;

	int SetQty(uint n, uint e);
	int AddNodes(uint n);
	int SetNodeUb(uint i, double nub);
	int SetNodeLb(uint i, double nlb);
	int SetNodeVal(uint i, double nval);
	int SetNodeName(uint i, string name);
	int SetEdgeUb(uint i, double nub);
	int SetEdgeLb(uint i, double nlb);
	int SetEdgeVal(uint i, double nval);
	int SetEdgeName(uint i, string name);
	int AddEdges(uint e);
	int RemoveEdge(uint e);
	int RemoveNode(uint e);
	//int Contract();
	int Clean(vector<bool>& exin, vector<bool>& exout);
	int PairNodeAndEdge(uint ni, uint ei, double t, bool force = false);

	vector<HyperGen> hnodes;
	vector<HyperGen> hedges;
};

class HyperGen
{
public:
	HyperGen() {val = 0.0; dall=0; din=0; dout = 0; ub = numeric_limits<double>::infinity();}
	HyperGen(double nval) {val = nval; dall=0; din=0; dout = 0; ub = numeric_limits<double>::infinity();}
	~HyperGen() {}
	uint dall;
	uint din;
	uint dout;
	// max val x; Ax=b; lb <= x <= ub; min b y + ub c+ + lb c- A'y + c+ - c- = cost
	double val; // cost for hedges, demand for hnodes
	double ub;  // ub for hedges, positive complement cost for hnodes
	double lb;  // lb for hedges, negative complement cost for hnodes
	string name;
//	GenFunction* f;
	vector<uint> all;
	vector<uint> in;
	vector<uint> out;
	vector<double> tau;
};

/*
class GenFunction
{
public:
	GenFunction() {};
	~GenFunction() {};
	virtual double f(double x) = 0;
//	virtual double invf(double x) = 0;
	virtual double auc(double x) = 0;
};
*/

#endif /* HYPERGRAPH_H_ */
