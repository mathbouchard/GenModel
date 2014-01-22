#include "GenModelHG.h"
//#include "ProblemReader.h"
#include <limits>

#define NUMMP 1
#if defined WIN64 || defined WIN32
    #include <windows.h>
    #define TRETURN DWORD WINAPI
    #define THREADEND return 0;
#else
    #include <pthread.h>
    #define TRETURN void*
    #define THREADEND return NULL;
#endif

TRETURN threadedsub(void* p)
{
	pair<uint, HGData*>* pp = (pair<uint, HGData*>*)p;
	pp->second->SolveLayer(pp->first);
    
	THREADEND;
}

long GenModelHG::Solve()
{

    printf("Solving... \n");
	HGData* d = static_cast<HGData*>(solverdata);

    d->incol = 300;
    d->numcol = 0;
    printf("\tInitialization\n");
	d->prim.Init(*this, d->extrac, d->lb, d->lbname, d->ub, d->ubname);
    
    printf("Init has (%d %d) [%d %d]\n", d->lb.size(), d->ub.size(), d->prim.lbindex, d->prim.ubindex);
    
	d->extcost.resize(d->hg.he, 0.0);
	
    printf("\tCreate first column\n");
	d->SolveHgMp();
    
    //printf("********* HG obj = %f | %f ***************\n", d->currobj.front(), d->realobj.front());

	
    d->coldeb = 0;
    
    printf("\tAdd first column\n");
	for(int i = 0; i < int(d->realobj.size()); i++)
		d->prim.AddVar(d->currsol[i], d->realobj[i], d->mippath[i]);
    d->numcol += d->realobj.size();
	d->prim.CreateProblem();
    
	printf("\tThe master problem has %d vars %d consts (%d %d) and %d nz [%d %d]\n", d->prim.gm.nc, d->prim.gm.nr, d->lb.size(), d->ub.size(), d->prim.gm.nz, d->prim.lbindex, d->prim.ubindex);
	
	d->prim.Solve(d->extcost, d->lb2he, d->ub2he, d->extrac, d->extrav);
	
    printf("First solve\n");
    
	uint iter = 1;
    
	//while(currobj > 10*numeric_limits<float>::epsilon())
	while(!d->currobj.empty())
	{
		d->SolveHgMp();
		//printf("********* HG obj = %f | %f ***************\n", d->currobj.front(), d->realobj.front());
		if(iter%1==0)
		{
			printf("Iter %d (%f / %f) [%f] [%d/%d] stat %d \n", iter, d->currobj.front(), d->realobj.front(), d->prim.gm.objval, int(d->currobj.size()), d->numcol, d->prim.gm.solstat);
		}
    
		for(int i = 0; i < int(d->realobj.size()); i++)
			d->prim.AddColAlt(d->currsol[i], d->realobj[i], d->mippath[i]);
        d->numcol += d->realobj.size();
        //printf("\tThe master problem has %d vars %d consts (%d %d) and %d nz [%d %d]\n", d->prim.gm.nc, d->prim.gm.nr, d->lb.size(), d->ub.size(), d->prim.gm.nz, d->prim.lbindex, d->prim.ubindex);
		d->extcost.assign(d->hg.he, 0.0);
		d->prim.Solve(d->extcost, d->lb2he, d->ub2he, d->extrac, d->extrav);
		
		iter++;
       // if(iter>15)
       //     break;
		//getchar();
	}
    objval = d->prim.gm.objval;
    
    printf("done.\n");
	   
	return 0;
}

long GenModelHG::SetSol()
{
	vars.sol.clear();
	vars.sol.resize(nc,0);
	vars.rc.clear();
	vars.rc.resize(nc,0);
	HGData* d = (HGData*)solverdata;
	

    double tmpobj = 0.0;
    double tmpobj2 = 0.0;
	for(uint i = d->coldeb; i < d->prim.gm.vars.n; i++)
	{
		if(fabs(d->prim.gm.vars.sol[i]) >= 1e-9)
		{
			for(uint j = 0; j < nc/*d->hg.he*/; j++)
			{
				if(fabs(d->colrep[i-d->coldeb][j]) >= 1e-9)
                {
                    vars.sol[j] += d->colrep[i-d->coldeb][j]*d->prim.gm.vars.sol[i];
                    tmpobj += -vars.obj[j]*d->colrep[i-d->coldeb][j]*d->prim.gm.vars.sol[i];
                    //printf("\t%d %d %f %f\n", i, j, d->prim.gm.vars.sol[i], d->hg.hedges[j].val);
                }
			}
            tmpobj2 += d->prim.gm.vars.obj[i]*d->prim.gm.vars.sol[i];
            //printf("%d %f %f\n", i, d->prim.gm.vars.sol[i], d->realobj[i-d->coldeb]);
		}
	}
    printf("tmpobj = %f %f (%d) %d\n", tmpobj, tmpobj2, d->prim.gm.vars.n, d->coldeb);
    for(int i = 0; i < nr; i++)//int(consts.size()); i++)
    {
        
        double total2 = 0.0;
        double total3 = 0.0;
        double total = 0.0;
        for(int j = 0; j < int(d->extrac[i].size()); j++)
            total2 += d->extrac[i][j].second*vars.sol[d->extrac[i][j].first];
        if(d->extrac[i].size() == 0)
        {
            total2 = numeric_limits<double>::infinity();
        }
        else
        {
            for(int j = 0; j < int(d->prim.gm.consts[d->prim.cons2extrac[i]].cols.size()); j++)
                //for(uint ii = d->coldeb; ii < d->prim.gm.vars.n; ii++)
                    //if(fabs(d->prim.gm.vars.sol[ii]) >= 1e-9)
                        total3 = d->prim.gm.consts[d->prim.cons2extrac[i]].coefs[j]*d->prim.gm.vars.sol[d->prim.gm.consts[d->prim.cons2extrac[i]].cols[j]];
        }
        
        for(int j = 0; j < int(consts[i].cols.size()); j++)
        {
            total += consts[i].coefs[j]*vars.sol[consts[i].cols[j]];
        }
        consts[i].slack = consts[i].lrhs-total;
        if(consts[i].slack < -1e-9 && consts[i].sense == 'L')
            printf("constraints %d (<=) violated %f > %f %f %f\n", i, total, consts[i].lrhs, total2, total3);
        else if(consts[i].slack > 1e-9 && consts[i].sense == 'G')
            printf("constraints %d (>=) violated %f < %f %f %f\n", i, total, consts[i].lrhs, total2, total3);
        if(fabs(consts[i].slack) > 1e-9 && consts[i].sense == 'E')
            printf("constraints %d (==) violated %f != %f %f %f\n", i, total, consts[i].lrhs, total2, total3);
        if(d->extrac[i].size() > 0)
            consts[i].dual = d->prim.gm.consts[d->prim.cons2extrac[i]].dual;
        else
            consts[i].dual = 0.0; // Voir valeur d'un chemin partant du noeud correspondant dans le hg
    }
        
	    
	return 0;
}


long GenModelHG::CreateModel()
{
	HGData* d = static_cast<HGData*>(solverdata);
    
	uint nvar=0;
	uint nconst = 0;
	uint mnz=0;
    
    vector<int> cneg;
    vector<int> varneg;
    vector<bool> extramap;
    vector<vector<pair<int,double> > > numneg;
    
    numneg.resize(nc);
    varneg.resize(nc,-1);
    cneg.resize(nr, 0);
    extramap.resize(nr, false);
    d->extrac.resize(nr);
    d->extrav.resize(nc);
    d->mipedges.resize(nc);
    
    d->hg.AddEdges(nc);
    d->hg.AddNodes(nr);
    
    d->objmult = 1.0;
    
    if(boolParam.count("inverseobj") > 0 && boolParam["inverseobj"] == true)
        d->objmult = -1.0;
    
    bool ismip = boolParam.count("mip") > 0 && boolParam["mip"] == true;
    
    printf("Create problem...\n");
    
	for(unsigned long i = 0; i < nc; i++)
	{
        //d->hg.hnodes[nconst].name = vars.name[i];
        d->hg.hedges[nvar].name = vars.name[i];
		d->hg.hedges[nvar].ub = vars.ub[i];
        d->hg.hedges[nvar].lb = vars.lb[i];
		d->hg.hedges[nvar].val = d->objmult*vars.obj[i];
        if(!ismip)
            d->mipedges[nvar] = 'C';
        else
            d->mipedges[nvar] = vars.type[i];
        nvar++;
        //nconst++;
    }
    
    printf("\tNodes added!\n");
    
    int countextra=0;
    
    for(unsigned long i = 0; i < nr; i++)
	{
        //if(consts[i].lrhs != 0 && consts[i].lrhs < 999999999)
        //    printf("0 extra added %c %f\n", consts[i].sense, consts[i].lrhs);
        for(unsigned long j = 0; j < int(consts[i].cols.size()); j++)
		{
            if((consts[i].lrhs != 0 || (consts[i].lrhs == 0 && consts[i].sense == 'L')) && consts[i].lrhs < 999999999)
            {
                //printf("0 extra added %d %c %f\n", consts[i].cols[j], consts[i].sense, consts[i].lrhs);
                //d->extrac[i].push_back(pair<int,double>(consts[i].cols[j], consts[i].coefs[j]));
                //d->extrav[consts[i].cols[j]].push_back(pair<int,double>(i, consts[i].coefs[j]));
                extramap[i] = true;
            }
            if(ismip && (vars.type[consts[i].cols[j]] == 'B' || vars.type[consts[i].cols[j]] == 'I') && consts[i].coefs[j] < 0)
                extramap[i] = true;
            if(consts[i].coefs[j] < 0)
                numneg[consts[i].cols[j]].push_back(pair<int,double>(i,consts[i].coefs[j]));
        }
        d->hg.hnodes[nconst].name = consts[i].name;
		//d->hg.hnodes[nconst].ub = vars.ub[i];
        //d->hg.hnodes[nconst].lb = vars.lrhs[i];
		d->hg.hnodes[nconst].val = consts[i].lrhs;
        nconst++;
    }
    
    for(unsigned long i = 0; i < nc; i++)
	{
        if(numneg[i].size() > 1)
            for(int j = 0; j < int(numneg[i].size()); j++)
                cneg[numneg[i][j].first]++;
    }
      
    for(unsigned long i = 0; i < nc; i++)
	{
        if(numneg[i].size() > 1)
        {
            map<int,int> isort;
            map<int,int>::iterator it;
            int count = 0;
            for(int j = 0; j < int(numneg[i].size()); j++)
            {
                if(count < int(numneg[i].size())-1 && extramap[numneg[i][j].first])
                {
                    printf("1 extra added\n");
                    count++;
                    d->extrac[extramap[numneg[i][j].first]].push_back(pair<int,double>(i,numneg[i][j].second));
                    d->extrav[i].push_back(pair<int,double>(numneg[i][j].first,numneg[i][j].second));
                }
                else
                    isort[cneg[numneg[i][j].first]] = j;//numneg[i][j].first;
            }
            it = isort.begin();
            varneg[i] = numneg[i][it->second].first;
            it++;
            for(; it != isort.end() && count < int(numneg[i].size())-1; it++)
            {
                printf("2 extra added\n");
                count++;
                d->extrac[numneg[i][it->second].first].push_back(pair<int,double>(i,numneg[i][it->second].second));
                d->extrav[i].push_back(pair<int,double>(numneg[i][it->second].first, numneg[i][it->second].second));
                extramap[numneg[i][it->second].first] = true;
            }
        }
    }
    
    printf("\textra constraints added\n");
    
    for(int i = 0; i < int(nr); i++)
	{
        if(extramap[i])
            countextra++;
        for(int j = 0; j < int(consts[i].cols.size()); j++)
        {
            if(extramap[i])
            {
                d->extrac[i].push_back(pair<int,double>(consts[i].cols[j],consts[i].coefs[j]));
                d->extrav[consts[i].cols[j]].push_back(pair<int,double>(i,consts[i].coefs[j]));
            }
            if(consts[i].coefs[j] >= 0 || numneg[consts[i].cols[j]].size() == 1 || varneg[consts[i].cols[j]] == i)
            {
                if(!ismip || (vars.type[consts[i].cols[j]] != 'B' && vars.type[consts[i].cols[j]] != 'I') || consts[i].coefs[j] >= 0)
                {
                    d->hg.PairNodeAndEdge(i, consts[i].cols[j], consts[i].coefs[j]);
                    mnz++;
                }
            }
		}
	}
    
    
    printf("\tthere is %d consts %d vars and %d extrac and %d nz\n", nr, nc, countextra, nz);
	
    printf("\thypergraph built with %d nodes %d edges and %d nz\n", d->hg.hn, d->hg.he, mnz);
    
    d->AddOrigin(*this);
    printf("\torigin added\n");
    
    printf("\tprocessing layers\n");
    d->CalcLayers();
    printf("\t%d hypegraph layers processed\n", d->nlayers.size());
    
    printf("done.\n");
    
	return 0;
}

int HGData::AddOrigin(GenModelHG& gm)
{
	//eorigin_index = hg.he;
	norigin_index = hg.hn;
	hg.AddNodes(1);
	hg.hnodes.back().name = "Origin";
	for(uint i = 0; i < hg.he; i++)
	{
		if(hg.hedges[i].din == 0)
		{
            //printf("orphelin edges %s\n", hg.hedges[i].name.c_str());
			hg.PairNodeAndEdge(norigin_index, i, -1.0);
		}
        
		if(hg.hedges[i].din > 1)
		{
			printf("there is more than one entry for hedge %d\n", i);
			return 1;
		}
	}
    
    for(uint i = 0; i < hg.hn-1; i++)
	{
		if(gm.consts[i].lrhs > 0 && hg.hnodes[i].dout > 0)
		{
            //printf("supply node %s\n", gm.consts[i].name.c_str());
			hg.AddEdges(1);
			hg.hedges.back().name = "ori_to_"+hg.hnodes[i].name;
			hg.hedges.back().val = 0.0;
			hg.hedges.back().ub = numeric_limits<double>::infinity();
            hg.hedges.back().lb = 0;
			hg.PairNodeAndEdge(norigin_index, hg.he-1, -1.0);
			hg.PairNodeAndEdge(i, hg.he-1, 1.0);
		}
	}
    
	he2lb.resize(hg.he);
    he2ub.resize(hg.he);
    
	for(uint i = 0; i < hg.he; i++)
	{
		if(hg.hedges[i].lb > 0)
		{
			he2lb[i] = (uint)(lb2he.size());
			lb2he.push_back(i);
			lb.push_back(hg.hedges[i].lb);
			lbname.push_back(hg.hedges[i].name+string("_lb"));
		}
		else
			he2lb[i] = -1;
        if(hg.hedges[i].ub !=  numeric_limits<double>::infinity())
		{
			he2ub[i] = (uint)(ub2he.size());
			ub2he.push_back(i);
			ub.push_back(hg.hedges[i].ub);
			ubname.push_back(hg.hedges[i].name+string("_ub"));
		}
		else
			he2ub[i] = -1;
	}
    
	return 0;
}

int HGData::CalcLayers()
{
	vector<bool> nvis;
	nvis.resize(hg.hn, false);
	vector<uint> nnext;
    
	list<uint> tocheck;
	for(uint i = 0; i < hg.hn; i++)
		tocheck.push_back(i);
    
	uint nodecount = 0;
    
	nlayers.push_back(vector<uint>());
	list<uint>::iterator it;
	for(it = tocheck.begin(); it != tocheck.end();)
	{
		if(hg.hnodes[*it].dout == 0)
		{
            //	printf("Layer -1: Adding %s (hnode) [%d %d %d]\n", hg.hnodes[*it].name.c_str(), hg.hnodes[*it].din,
            //			hg.hnodes[*it].dout, hg.hnodes[*it].dall);
			nlayers.back().push_back(*it);
			nvis[*it] = true;
			nodecount++;
			it = tocheck.erase(it);
		}
		else
			it++;
	}
    
    printf("\t\tthere is %d nodes to checked\n", tocheck.size());
    
	uint iter = 0;
	bool cycle = false;
    
	while(!tocheck.empty() && iter < hg.hn)
	{
		nlayers.push_back(vector<uint>());
        
		vector<bool> tempvis;
		tempvis.resize(hg.hn, true);
        
        
		for(it = tocheck.begin(); it != tocheck.end();)
		{
			uint in = *it;
            
			for(uint j = 0; j < hg.hnodes[in].dout; j++)
			{
				uint ie = hg.hnodes[in].all[hg.hnodes[in].out[j]];
				for(uint k = 0; k < hg.hedges[ie].dout; k++)
				{
					uint inn = hg.hedges[ie].all[hg.hedges[ie].out[k]];
					if(!nvis[inn])
					{
						tempvis[in] = false;
                        //	printf("	%s (%s %d) %s\n", hg.hnodes[in].name.c_str(), hg.hedges[ie].name.c_str(), id2proc[ie], hg.hnodes[inn].name.c_str());
					}
				}
			}
			if(tempvis[in])
			{
				//printf("Layer %d: Adding %s (hnode) [%d %d %d]\n", iter, hg.hnodes[in].name.c_str(), hg.hnodes[in].din,
				//					hg.hnodes[in].dout, hg.hnodes[in].dall);
				nnext.push_back(in);
				nlayers.back().push_back(in);
				it = tocheck.erase(it);
				nodecount++;
			}
			else
				it++;
		}
		if(nnext.empty())
		{
			printf("cycle true\n");
			cycle = true;
			break;
		}
        
		for(uint i = 0; i < (uint)(nnext.size()); i++)
			nvis[nnext[i]] = true;
		nnext.clear();
        
		iter++;
		//printf("\t\titer = %d (layers = %d)\n", iter, nlayers.size());
	}
    printf("\t\tcycle detection...", iter, nlayers.size());
	
	if(iter == hg.he || cycle)
	{
		printf("\nThere is at least one cycle in the hypergraph (%d,%d [%d,%d] %d)\n", hg.he, (uint)(tocheck.size()), nodecount, hg.hn, iter);
		printf("%s\n", hg.hnodes[tocheck.back()].name.c_str());
		return 1;
	}
	else
		printf(" done.\n\t\tLayers constructed\n");
    
	return 0;
}

//long GenModelHG::CreateModel(string filename, int type, string dn)
//{
//    ReadFromFile(static_cast<GenModel*>(this), filename, type);
//    if(boolParam.count("maximize") > 0 && boolParam["maximize"] == true)
//        boolParam["inverseobj"] = true;
//    SetNumbers();
//    CreateModel();
//    
//    return 0;
//}

long GenModelHG::Init(string name)
{
	if(solverdata == NULL)
		solverdata = new HGData();
	
	return 0;
}

long GenModelHG::Clean()
{
	if(solverdata != NULL)
		delete static_cast<HGData*>(solverdata);

	return 0;
}

int HGData::SolveLayer(uint p)
{
	uint i = curr_layer;
	for(uint j = p; j < (uint)(nlayers[i].size()); j+=NUMMP)
	{
		long in = nlayers[i][j];
		
		if(hg.hnodes[in].dout == 0)
			inval[in] = 0.0;
		for(uint k = 0; k < hg.hnodes[in].dout; k++)
		{
			double alpha = fabs(hg.hnodes[in].tau[hg.hnodes[in].out[k]]);
			uint ie = hg.hnodes[in].all[hg.hnodes[in].out[k]];
			double newval = (hg.hedges[ie].val+extcost[ie])/alpha;
			for(uint l = 0; l < (uint)(hg.hedges[ie].dout); l++)
			{
				uint jn = hg.hedges[ie].all[hg.hedges[ie].out[l]];
				double beta = fabs(hg.hedges[ie].tau[hg.hedges[ie].out[l]]);
				newval+=inval[jn]/alpha*beta;
			}
			if(newval > inval[in])
			{
				inval[in] = newval;
				insol[in] = ie;
			}
		}
	}
	return 0;
}

int HGData::SolveHgMp()
{
	//vector<uint> sol;
	//insol.clear();
	if(insol.empty())
		insol.resize(hg.hn, -1);//numeric_limits<unsigned int>::max());
	else
		insol.assign(hg.hn, -1);//numeric_limits<unsigned int>::max());
	//vector<double> val;
	if(inval.empty())
		inval.resize(hg.hn, -numeric_limits<double>::infinity());
	else
		inval.assign(hg.hn, -numeric_limits<double>::infinity());
    
	//printf("SolveHgMp\n");
    
	for(int i = 0; i < int(nlayers.size()); i++)
	{
		curr_layer = i;
        //		for(uint j = 0; j < NUMMP && j < (uint)(nlayers[i].size()) ; j++)
        //			SolveLayer(j);
		pair<uint, HGData*> valpair[NUMMP];
		for(uint j=0; j < NUMMP; j++)
			valpair[j] = pair<uint, HGData*>(j, this);
        
        
#if defined WIN64 || defined WIN32
        HANDLE  hThreadArray[NUMMP];
        DWORD   dwThreadIdArray[NUMMP];
        for(uint j=0; j < NUMMP; j++)
            hThreadArray[j] = CreateThread(NULL, 0, threadedsub, &valpair[j], 0, &dwThreadIdArray[j]);
        
        WaitForMultipleObjects(NUMMP, hThreadArray, TRUE, INFINITE);
#else
        pthread_t tsubs[NUMMP];
        for(uint j=0; j < NUMMP; j++)
            pthread_create( &tsubs[j], NULL, threadedsub, &valpair[j]);
        for(uint j=0; j < NUMMP; j++)
            pthread_join(tsubs[j], NULL);
#endif
	}
    
	currobj.clear();
	currsol.clear();
	realobj.clear();
    mippath.clear();
    
    //getchar();
    
	for(int k = 0; k < int(hg.hnodes[norigin_index].dout); k++)
	{
		int inext = hg.hnodes[norigin_index].all[hg.hnodes[norigin_index].out[k]];
        
		double vobj = hg.hedges[inext].val+extcost[inext];
        
        //printf("voisin %d = %f (%d, %d)\n", k, vobj, norigin_index, insol[norigin_index]);
        
        //	if(inext != insol[norigin_index])
		//	continue;
        
		for(int vn = 0; vn < int(hg.hedges[inext].dout); vn++)
			vobj += inval[hg.hedges[inext].all[hg.hedges[inext].out[vn]]]*hg.hedges[inext].tau[hg.hedges[inext].out[vn]];
        
		//printf("comp %f %f\n", vobj, inval[norigin_index]);
        
		if(vobj < 10*numeric_limits<float>::epsilon())// 1e-8)
			continue;
        
        //printf("voisin %d = %f (%d, %d)\n", k, vobj, norigin_index, insol[norigin_index]);
        
		insol[norigin_index] = inext;
        
		currsol.push_back(vector<double>());
        
		list<uint> curr;
		list<double> calpha;
		curr.push_back(norigin_index);
		calpha.push_back(1.0);
		uint level = 0;
        
		//currobj = inval[norigin_index];
		currobj.push_back(vobj);
		//currsol.assign((uint)(ub.size()), 0.0);
		currsol.back().resize((uint)(prim.ubindex+ub.size()), 0.0);
        
		//realobj = 0.0;
		realobj.push_back(0.0);
		colrep.push_back(vector<double>(hg.he, 0.0));
        if(mipedges[inext])
            mippath.push_back(true);
        else
            mippath.push_back(false);
        
		while(!curr.empty())
		{
            //printf("level = %d (%d %d) %d %d\n", level, curr.size(), calpha.size(), curr.front(), insol[curr.front()]);
			level++;
			list<uint>::iterator it;
			list<double>::iterator ita = calpha.begin();
			uint max = (uint)(curr.size());
			uint ii = 0;
			for(it = curr.begin(); ii < max; ii++)
			{
                //printf("\t%d/%d\n", ii,max);
                if(insol[*it] == -1)
                {
                    it = curr.erase(it);
                    ita = calpha.erase(ita);
                    //if(it != curr.end())
                    //    printf("\t*next it = %d %d next a = %f\n", *it, insol[*it], *ita);
                    continue;
                }
                
                //printf("%d/%d (a = %f) %f %f\n", insol[*it], hg.he, *ita, (hg.hedges[insol[*it]].val)*(*ita)/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]),
                //       (extcost[insol[*it]])*(*ita)/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]));
                
                for(int ij = 0; ij < int(extrav[insol[*it]].size()); ij++)
                {
                    //printf("%d Basic %d %d %d (%f %f %f)\n", realobj.size(), *it, level, ii, extcost[insol[*it]], hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]], extrav[insol[*it]][ij].second);
                    //printf("%d %d Adding %f to %d (var %d) %d\n", realobj.size(), ij, *ita*extrav[insol[*it]][ij].second/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]),
                    //       prim.cons2extrac[extrav[insol[*it]][ij].first], insol[*it], extrav[insol[*it]][ij].first);
                    currsol.back()[prim.cons2extrac[extrav[insol[*it]][ij].first]] += *ita*extrav[insol[*it]][ij].second/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]);
                }
                if(he2lb[insol[*it]] != -1)
					currsol.back()[prim.lbindex+he2lb[insol[*it]]] += *ita/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]);
				if(he2ub[insol[*it]] != -1)
                {
                    //printf("Basic %d %d (%f %f)\n", hg.hn, hg.he, extcost[insol[*it]], hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]);
                    //printf("Adding %f to %d (var %d)\n", *ita/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]), prim.ubindex+he2ub[insol[*it]], insol[*it]);
					currsol.back()[prim.ubindex+he2ub[insol[*it]]] += *ita/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]);
                }
              
                //currsol[he2ub[insol[*it]]] += *ita/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]);
				colrep.back()[insol[*it]] += *ita/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]);
                
				//realobj+=(hg.hedges[insol[*it]].val)
                //			*(*ita)/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]);
				realobj.back()+=(hg.hedges[insol[*it]].val)*(*ita)/fabs(hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]);
              
                //printf("\tdout = %d\n", hg.hedges[insol[*it]].dout);
				for(uint i = 0; i < hg.hedges[insol[*it]].dout; i++)
				{
					uint ni = hg.hedges[insol[*it]].all[hg.hedges[insol[*it]].out[i]];
					
                    //if(hg.hedges[insol[*it]].din == 0)
                    //{
                    //    printf("din = 0\n");
                    //}
                    
                    //printf("\t\t%d %d %d %f %f i=%d/%d %d %d %d\n", insol[*it], hg.hedges[insol[*it]].dout, hg.hedges[insol[*it]].din, hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].out[i]], hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]], i,hg.hedges[insol[*it]].dout, ii, max, curr.size());
                    
					calpha.push_back(fabs((*ita)*hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].out[i]]
                                          /hg.hedges[insol[*it]].tau[hg.hedges[insol[*it]].in[0]]));
                    curr.push_back(ni);
                    
				}
				it = curr.erase(it);
				ita = calpha.erase(ita);
                //if(it != curr.end())
                //    printf("\tnext it = %d %d next a = %f\n", *it, insol[*it], *ita);
			}
		}
        if(realobj.size() > incol)
            break;
	}
    
    //getchar();
    
	return 0;
}

int PrimalProblem::Init(GenModelHG& gmhg, vector<vector<pair<int,double> > > exc, vector<double>& lb, vector<string>& lbname, vector<double>& ub, vector<string>& ubname)
{
    cindex = 0;
    lbindex = 0;
    ubindex = 0;
    cons2extrac.resize(exc.size(), -1);
    for(int i = 0; i < (uint)(exc.size())/*int(gmhg.nr)*/; i++)
    {
        if(exc[i].size() > 0)
        {
            //printf("cons %d %d\n", i, exc[i].size());
            cons2extrac[i] = lbindex;
            gm.AddConst(gmhg.consts[i].name , gmhg.consts[i].lrhs, gmhg.consts[i].sense);
            lbindex++;
            ubindex++;
        }
    }    
    
    for(uint i = 0; i < (uint)(lb.size()); i++)
    {
		gm.AddConst(lbname[i], lb[i], 'G');
        ubindex++;
    }
    for(uint i = 0; i < (uint)(ub.size()); i++)
		gm.AddConst(ubname[i], ub[i], 'L');
        
	
	return 0;
}

int PrimalProblem::AddVar(vector<double>& sol, double val, char type)
{
	char tmp[256];
    
	snprintf(tmp, 256, "col_%ld", gm.vars.n);
	gm.AddVar(string(tmp), val, 0.0, numeric_limits<double>::infinity(), type);
    gm.SetNumbers();
	for(uint i = 0; i < (uint)(sol.size()); i++)
	{
		if(fabs(sol[i]) > 10*numeric_limits<double>::epsilon())
		{
            //printf("%d/%d %d/%d %f\n", i, gm.nr, gm.vars.n-1, gm.nr, sol[i]);
			gm.AddNz(i, gm.vars.n-1, sol[i]);
		}
	}
    
	return 0;
}

int PrimalProblem::AddColAlt(vector<double>& sol, double val, char type)
{
	char tmp[256];
	snprintf(tmp, 256, "col_%ld", gm.vars.n);
    
	vector<int> inds;
	vector<double> vals;
    
	int nz = 0;
    
	for(uint i = 0; i < (uint)(sol.size()); i++)
	{
		if(fabs(sol[i]) > 10*numeric_limits<double>::epsilon())
		{
			inds.push_back(i);
			vals.push_back(sol[i]);
			nz++;
		}
	}
    
    ((GenModelCplex*)(&gm))->AddSolverCol(inds, vals, val, 0.0, numeric_limits<double>::infinity(), string(tmp), type);
    
	gm.SetNumbers();
    
	return 0;
}

int PrimalProblem::Solve(vector<double>& dual, vector<int>& lb2he, vector<int>& ub2he, vector<vector<pair<int,double> > >& exc, vector<vector<pair<int,double> > >& exv)
{
	((GenModelCplex*)(&gm))->Solve();
	((GenModelCplex*)(&gm))->SetSol();
	
    int cc = 0;
    for(uint i = 0; i < (uint)(exc.size()); i++)
    {
        if(exc[i].size() > 0)
        {
            for(int j = 0; j < int(exc[i].size()); j++)
            {
                //if(-gm.consts[cindex+cc].dual != 0)
                //    printf("dual edge extra %d = %f (vars %d)\n", cc, -exc[i][j].second*gm.consts[cindex+cc].dual, exc[i][j].first);
                if(fabs(-exc[i][j].second*gm.consts[cindex+cc].dual) > 1000000.0)
                    printf("%d/%d strange %f %f %f\n", cindex+cc, ubindex, exc[i][j].second, -gm.consts[cindex+cc].dual, -exc[i][j].second*gm.consts[cindex+cc].dual);
                    
                dual[exc[i][j].first] += -exc[i][j].second*gm.consts[cindex+cc].dual;
            }
            cc++;
        }
    }
    
    for(uint i = 0; i < (uint)(lb2he.size()); i++)
		dual[lb2he[i]] = gm.consts[lbindex+i].dual;
    
	for(uint i = 0; i < (uint)(ub2he.size()); i++)
    {
        if(-gm.consts[ubindex+i].dual != 0)
            printf("dual edge ub %d %d = %f\n", ubindex+i, ub2he[i], -gm.consts[ubindex+i].dual);
		dual[ub2he[i]] = -gm.consts[ubindex+i].dual;
    }
    
	return 0;
}

int PrimalProblem::CreateProblem()
{
	gm.SetNumbers();
	gm.boolParam["maximize"] = true;
	gm.boolParam["screenoff"] = true;
	gm.strParam["algo"] = "concurrent";
	
    ((GenModelCplex*)(&gm))->Init("PrimalProb");
    ((GenModelCplex*)(&gm))->CreateModel();
	
	return 0;
}

long GenModelHG::AddSolverRow(vector<int>& ind, vector<double>& val, double rhs, char sense, string name)
{
    return 0;
}

long GenModelHG::AddCut(int* cols, double* vals, int nz, double rhs, char sense, const char* name)
{
    return 0;
}

long GenModelHG::AddSolverCol(vector<int>& ind, vector<double>& val, double obj, double lb, double ub, string name, char type)
{
	return 0;
}

long GenModelHG::AddCol(int* newi, double* newcol, int nz, double obj, double lb, double ub, const char* name, char type)
{
	return 0;
}

long HGData::Reset()
{
	return 0;
}

HGData::HGData()
{
	Reset();
}

HGData::~HGData()
{
	Delete();
}
long HGData::Delete()
{
	return 0;
}
