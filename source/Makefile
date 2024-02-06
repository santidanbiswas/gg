################################################################################
#  File :          Makefile (GNU LINUX makefile)                               #
################################################################################
#!/bin/sh

#local source dir:
SRCDIR= .
#Directory path of source:
BASEDIR=..

#Executable extension:
EXT=.X
#Program name:
PRO=TATA
#executable target:
NAME=$(PRO)$(EXT)

#Fortran-90 compilers and flags:
#FF = f90
FF = ifc
OPTIONS = -o
CFLAGS = -c -g

#Destination of the final binary file:
DEST=$(NAME)

#Extensions:
MOD = mod
.SUFFIXES: .f90 .o .$(MOD)

#SHELL commands:
RM= /bin/rm -f
MV= /bin/mv -f

#Objects for $(NAME):
OBJS = moduleGlobal.o             \
       moduleInput.o              \
       moduleRand.o               \
       modulePotMinima.o          \
       moduleSites.o              \
       moduleCalculation.o        \
       tata.o              

$(NAME): $(OBJS)  
	@echo ""
	@echo "Linking object files [ $^ ] ===>>> [ $@ ]"
	$(FF) $(OPTIONS) $(DEST) $(OBJS) 
#
%.o : %.f90
	@echo ""
	@echo "Linking object files [ $^ ] ===>>> [ $@ ]"
	@echo "Outdated dependencies : [ $? ]"
	$(FF) $(CFLAGS) $<
#
%.$(MOD):
	@echo ""
	@echo "Module file [ $@ ] considered updated through dependence on : [ $^ ]"
	touch  $@ 
#
clean:
	$(RM) $(OBJS) *.$(MOD) 
#
# Dependencies for creating object files (USE statements)
classInput.o           : moduleGlobal.$(MOD)
moduleMain.o           : moduleGlobal.$(MOD)
moduleSites.o          : moduleInput.$(MOD)  
modulePotMinima.o      : moduleInput.$(MOD)
moduleCalculation.o    : moduleSites.$(MOD)

# Location of modules in source files:
moduleGlobal.$(MOD)          : moduleGlobal.o 
moduleInput.$(MOD)           : moduleInput.o 
moduleRand.$(MOD)            : moduleRand.o
moduleSites.$(MOD)           : moduleSites.o
modulePotMinima.$(MOD)       : modulePotMinima.o
moduleCalculation.$(MOD)     : moduleCalculation.o

