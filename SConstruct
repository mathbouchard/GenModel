# Tell SCons to create our build files in the 'build' directory

AddOption("--cplex", action="store_true", dest="cplex", default=False)
AddOption("--osi", "--cbc", "--coinor", action="store_true", dest="osi", default=False)
AddOption("--glpk", action="store_true", dest="glpk", default=False)
AddOption("--grb", "--gurobi", action="store_true", dest="grb", default=False)
AddOption("--scip", action="store_true", dest="scip", default=False)
AddOption("--hg", "--hypergraph", action="store_true", dest="hg", default=False)

from platform import system
print("System is "+ system());
VariantDir('build-release', 'src', duplicate=0)
SConscript('build-release/SConscript')
VariantDir('build-debug', 'src', duplicate=0)
SConscript('build-debug/SConscript')
