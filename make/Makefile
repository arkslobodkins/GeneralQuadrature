# Arkadijs Slobodkins
# SMU Mathematics
# February 2021
#

CC = gcc
CFLAGS=-g
LIBSLA = -llapack
LIBSBLAS= -lblas
MATH= -lm
AR = ar
RANLIB = ranlib
VPATH= ../src:../src/Gauss_Lib:../results
OBJS =  GaussCore.o Hermite.o Jacobi.o Laguerre.o 

all : Driver.x 
Driver.x: Driver.o TestIntegral.o InsertionSort.o IntegralsOfBasisFunctions.o LeastSquaresNewton.o GaussTensor.o GetJacobian.o GetFunction.o NodeElimination.o  Phi.o  LegendrePoly.o BasisIndices.o Output.o JacobiPoly.o AddDimension.o InDomain.o ComputeDomain.o ImplementDomain.o SetParams.o SetDomain.o StoreCopy.o libgaussq.a
	$(CC)  $(CFLAGS)  -o $@ $^ $(LIBSLA) $(LIBSBLAS) $(MATH)
	
libgaussq.a: $(OBJS)
	$(AR) rv libgaussq.a $(OBJS)
	$(RANLIB) libgaussq.a
	
clean :
	\rm -f *.o  

realclean : clean
	\rm -f *.x *~ *.a

SUBDIR_ROOTS := ../results
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
GARBAGE_PATTERNS := *.txt *.dat
GARBAGE := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(GARBAGE_PATTERNS)))
txtclean: 
	rm -rf $(GARBAGE)

