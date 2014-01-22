from platform import system
from os import listdir
from collections import OrderedDict

def validateLibrary(libraryName):
    """ Search in the libpath of the current platform if the library is present """
    for folder in libpath[system()]:
        for file in listdir(folder):
            if (libraryName + "." in file):
                return True
    return False

cpppath = {
    "Linux": [r"/opt/ibm/ILOG/CPLEX_Studio125/cplex/include",
        r"/home/mbouchard/source/gurobi510/linux64/include",
        r"/usr/local/include/"],
    "Windows": [r"C:\ILOG\CPLEX_Studio125\cplex\include",
        r"C:\gurobi550\win64\include",
        r"D:\Source\glpk-4.52\src"]
    }
        
libpath = {
    "Linux": [r"/opt/ibm/ILOG/CPLEX_Studio125/cplex/bin/x86-64_sles10_4.1/",
        r"/home/mbouchard/source/gurobi510/linux64/lib",
        r"/usr/local/lib/"],
    "Windows": [r"C:\ILOG\CPLEX_Studio125\cplex\lib\x64_windows_vs2010\stat_mta",
        r"C:\gurobi550\win64\lib",
        r"D:\Source\glpk-4.52\w64"]
    }

cppdefines = []
if (system() ==     "Windows"):
   cppdefines.append("WIN64")

libs = []

# Validate that the libs are present
solverLibs = {
    "Cplex": {
        "Linux": ["cplex125"],
        "Windows": ["cplex125"]
        },
    "Gurobi": {
        "Linux": ["gurobi_c++", "gurobi51"],
        "Windows": ["gurobi_c++md2012", "gurobi55"]
        },
    "HG": {
        "Linux": ["cplex125"],
        "Windows": ["cplex125"]
        },
    "Glpk": {
        "Linux": ["glpk"],
        "Windows": ["glpk_4_52"]
        }
    }

missingSolvers = []

for solver in solverLibs:
    allLibPresents = True
    for lib in solverLibs[solver][system()]:
        if (not validateLibrary(lib)):
            allLibPresents = False
    if (allLibPresents):
        libs += solverLibs[solver][system()];
        cppdefines.append((solver + "_SOLVER").upper())
    else:
        missingSolvers.append(solver)

libs = list(OrderedDict.fromkeys(libs))

cxxflags = {
    "Linux": "-O2 -m64 -Wall -fmessage-length=0",
    "Windows": "/O2 /MD /EHsc"
    }
    
# Create an environmnet
env = Environment(TARGET_ARCH = "x86_64",
    CPPPATH = cpppath[system()],
    LIBPATH = libpath[system()],
    LIBS = libs,
    CXXFLAGS = cxxflags[system()],
    CPPDEFINES = cppdefines
    )

# Change the build folder output
env.VariantDir("build", "src", duplicate = 0)

files = Glob("build/*.cpp")
if (system() == "Windows"):
    files += ["build/GenModelDLL.def"]
    target = "GenModelDLL"
elif (system() == "Linux"):
    target = "GenModelLib"

# Remove the cpp files for the missing solvers
filesToRemove = []
if (not env.GetOption("clean")):
    for file in files:
        for solver in missingSolvers:
            if solver in str(file):
                filesToRemove.append(file)

for file in filesToRemove:
    print(file)
    files.remove(file)
    

# Build the library
env.SharedLibrary(target = "lib/" + target, source = files)

if (missingSolvers):
    print("***************************************************************************")
    print("****** The following solvers will not compiled into the final library *****")
    print(missingSolvers)
    print("***************************************************************************")
