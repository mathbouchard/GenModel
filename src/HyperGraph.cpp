/*
 * HyperGraph.cpp
 *
 *  Created on: 2011-03-17
 *      Author: mbouchard
 */

#include <stdio.h>
#include <stdlib.h>
#include "HyperGraph.h"

void checkbound(uint i, uint ub);

int HyperGraph::SetQty(uint n, uint e)
{
	hnodes.resize(n);
	hnodes.assign(n, HyperGen(0.0));

	hedges.resize(n);
	hedges.assign(n, HyperGen(0.0));

	hn = n;
	he = e;

	return 0;
}


int HyperGraph::PairNodeAndEdge(uint ni, uint ei, double t, bool force)
{
	if(!force && t == 0.0)
		return 0;

	checkbound(ei, he);
	checkbound(ni, hn);
	hedges[ei].all.push_back(ni);
	hedges[ei].tau.push_back(t);
	hedges[ei].dall++;
	hnodes[ni].all.push_back(ei);
	hnodes[ni].tau.push_back(t);
	hnodes[ni].dall++;

	if(t < 0)
	{
		hedges[ei].in.push_back(hedges[ei].dall-1);
		hedges[ei].din++;
		hnodes[ni].out.push_back(hnodes[ni].dall-1);
		hnodes[ni].dout++;
	}
	else
	{
		hedges[ei].out.push_back(hedges[ei].dall-1);
		hedges[ei].dout++;
		hnodes[ni].in.push_back(hnodes[ni].dall-1);
		hnodes[ni].din++;
	}

	return 0;
}

int AddGens(uint nq, vector<HyperGen>* vhg, uint* cq)
{
	if(nq == 1)
		vhg->push_back(HyperGen(0.0));
	else
		vhg->resize(*cq+nq, HyperGen(0.0));
	(*cq)+=nq;
	return 0;
}

int SetGenVal(uint i, double nval, vector<HyperGen>* vhg, uint s)
{
	checkbound(i, s);
	(*vhg)[i].val = nval;
	return 0;
}

int SetGenUb(uint i, double nub, vector<HyperGen>* vhg, uint s)
{
	checkbound(i, s);
	(*vhg)[i].ub = nub;
	return 0;
}

int SetGenLb(uint i, double nlb, vector<HyperGen>* vhg, uint s)
{
	checkbound(i, s);
	(*vhg)[i].lb = nlb;
	return 0;
}

int SetGenName(uint i, string name, vector<HyperGen>* vhg, uint s)
{
	checkbound(i, s);
	(*vhg)[i].name = name;
	return 0;
}

int RemoveGen(uint n1, vector<HyperGen>* hv1, vector<HyperGen>* hv2, uint* hq)
{
	for(uint i = 0; i < (*hv1)[n1].dall; i++)
	{
		uint n2 = (*hv1)[n1].all[i];
		vector<uint>::iterator ita = (*hv2)[n2].all.begin();
		vector<double>::iterator itt = (*hv2)[n2].tau.begin();
		vector<uint>::iterator iti = (*hv2)[n2].in.begin();
		vector<uint>::iterator ito = (*hv2)[n2].out.begin();

		uint remindex = 0;
		uint indexcnt = 0;
		while(ita != (*hv2)[n2].all.end())
		{
			if(n1 == *ita)
			{
				ita = (*hv2)[n2].all.erase(ita);
				itt = (*hv2)[n2].tau.erase(itt);
				remindex=indexcnt;
			}
			else
			{
				ita++;
				itt++;
				indexcnt++;
			}
		}

		while(iti != (*hv2)[n2].in.end())
		{
			if(remindex == *iti)
				iti = (*hv2)[n2].in.erase(iti);
			else
			{
				if(*iti > remindex)
					(*iti)--;
				iti++;
			}
		}

		while(ito != (*hv2)[n2].out.end())
		{
			if(remindex == *ito)
				ito = (*hv2)[n2].out.erase(ito);
			else
			{
				if(*ito > remindex)
					(*ito)--;
				ito++;
			}
		}
	}

	vector<HyperGen>::iterator it = hv1->begin();
	for(uint i = 0; i < n1; i++)
		it++;

	hv1->erase(it);
	(*hq)--;

	return 0;
}

int HyperGraph::AddNodes(uint n)
{
	return AddGens(n, &hnodes, &hn);
}

int HyperGraph::SetNodeVal(uint i, double nval)
{
	return SetGenVal(i, nval, &hnodes, hn);
}

int HyperGraph::SetNodeUb(uint i, double nub)
{
	return SetGenUb(i, nub, &hnodes, hn);
}

int HyperGraph::SetNodeLb(uint i, double nlb)
{
	return SetGenLb(i, nlb, &hnodes, hn);
}


int HyperGraph::SetNodeName(uint i, string name)
{
	return SetGenName(i, name, &hnodes, hn);
}

int HyperGraph::AddEdges(uint e)
{
	return AddGens(e, &hedges, &he);
}

int HyperGraph::SetEdgeVal(uint i, double nval)
{
	return SetGenVal(i, nval, &hedges, he);
}

int HyperGraph::SetEdgeUb(uint i, double nub)
{
	return SetGenUb(i, nub, &hedges, he);
}

int HyperGraph::SetEdgeLb(uint i, double nlb)
{
	return SetGenLb(i, nlb, &hedges, he);
}

int HyperGraph::SetEdgeName(uint i, string name)
{
	return SetGenName(i, name, &hedges, he);
}

int HyperGraph::RemoveEdge(uint e)
{
	return RemoveGen(e, &hedges, &hnodes, &he);
}
int HyperGraph::RemoveNode(uint n)
{
	return RemoveGen(n, &hnodes, &hedges, &hn);
}

int HyperGraph::Clean(vector<bool>& exin, vector<bool>& exout)
{
	vector<uint> nstatus;
	vector<uint> ndin;
	vector<uint> ndout;
	for(long i = 0; i < hn; i++)
	{
		ndin.push_back(hnodes[i].din);
		ndout.push_back(hnodes[i].dout);
	}
	nstatus.resize(hn, 0);
	vector<uint> estatus;
	vector<uint> edin;
	vector<uint> edout;
	for(long i = 0; i < he; i++)
	{
		edin.push_back(hedges[i].din);
		edout.push_back(hedges[i].dout);
	}
	estatus.resize(he, 0);

	vector<uint> nfinal;
	for(long i = 0; i < hn; i++)
		nfinal[i] = i;


	uint erem = 1;
	uint nrem = 1;
	uint terem = 1;
	uint tnrem = 1;

//	list<double> rn;
//	for(long i = 0; i < hn; i++)
//		rn.push_back(i);
//	list<double> re;
//	for(long i = 0; i < he; i++)
//		re.push_back(i);

	while(erem > 0 || nrem > 0)
	{
		erem = 0;
		nrem = 0;
		// Analysing nodes
		for(uint i = 0; i < hn; i++)
		{
			if(nstatus[i] > 0)
				continue;
			if(ndin[i] == 0 && !exout[i])
			{
				nstatus[i] = 1;
				for(uint j = 0; j < hnodes[i].dout; j++)
					edin[hnodes[i].all[hnodes[i].out[j]]]--;
				nrem++;
				tnrem++;
			}
			else if(ndout[i] == 0 && !exin[i])
			{
				nstatus[i] = 2;
				for(uint j = 0; j < hnodes[i].din; j++)
					edin[hnodes[i].all[hnodes[i].in[j]]]--;
				nrem++;
				tnrem++;
			}
		}

		// Analysing edges
		for(uint i = 0; i < he; i++)
		{
			if(estatus[i] > 0)
				continue;
			if(edin[i] == 0 || edout[i] == 0)
			{
				estatus[i] = 1;
				for(uint j = 0; j < hedges[i].dout; j++)
					ndin[hedges[i].all[hedges[i].out[j]]]--;
				for(uint j = 0; j < hedges[i].din; j++)
					ndout[hedges[i].all[hedges[i].in[j]]]--;

				erem++;
				terem++;
			}
		}
	}

	printf("%d nnodes removed, %d edges removed\n", tnrem, terem);

	return 0;
}

void checkbound(uint i, uint ub)
{
	if(i > ub)
	{
		printf("indice %d out of range (>= %d) \n", i, ub);
		exit(0);
	}
	return;
}
