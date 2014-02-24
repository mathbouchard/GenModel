/*
 * ProblemReaderOsi.h
 *
 *  Created on: 2012-10-02
 *      Author: mbouchard
 */

#ifndef PROBLEMREADER_H_
#define PROBLEMREADER_H_

#if defined WIN64 || defined WIN32
	#ifndef snprintf
		#define snprintf sprintf_s
	#endif
#endif

#include "GenModel.h"
#include "OsiClpSolverInterface.hpp"
//#include "OsiGlpkSolverInterface.hpp"
#include <string>

using namespace std;

int ReadFromFile(GenModel* pgm, string fn, int type=0, string dn = "");
int ReadFromObject(GenModel* pgm, OsiSolverInterface* pm);

#endif /* PROBLEMREADER_H_ */
