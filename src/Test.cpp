/***************************************************************************
 *  BaseApp.cpp
 *  An example application
 *
 *  January 8 11:32 2014
 *  Copyright  2014  Mathieu Bouchard
 *  mathbouchard@gmail.com
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "BaseInterface.h"

using namespace std;

int main(int argc, char** argv)
{
    string valo = "default";
    if(argc > 1)
        valo = string(argv[1]);
    string vali = string("interface-")+valo;
    
#if defined BASELIB_MODULE
    printf("Testing module object\n");
    printf("---------------------\n");
    BaseLib bl;
    bl.Run();
    bl.Set(valo);
    bl.Run();

    printf("\nTesting module interface\n");
    printf("---------------------\n");
    int token = _CreateModule();
    _Run(token);
    _Set(vali.c_str(), token);
    _Run(token);
#else
    printf("BaseLib Module not present\n");
#endif
    
    return 0;
}
