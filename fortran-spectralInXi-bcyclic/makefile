EXTRA_COMPILE_FLAGS =
EXTRA_LINK_FLAGS = 

FC = ftn
FLINKER = ftn

TARGET = mmc

OBJ_FILES = computeB.o createGrids.o diagnostics.o indices.o mmc.o populateMatrix.o populateRHS.o stel_constants.o stel_kinds.o uniformDiffMatrices.o variables.o cyclicr_MPI_real.o

.PHONY: all clean

all: $(TARGET)

#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

$(TARGET): $(OBJ_FILES)
	$(FLINKER) -o $(TARGET) $(OBJ_FILES) $(EXTRA_LINK_FLAGS) $(PETSC_LIB)

clean::
	rm -f *.o *.mod *.a *.dat *~ $(TARGET)

# Dependencies:
computeB.o: variables.o
createGrids.o: stel_kinds.o variables.o
cyclicr_MPI_real.o: stel_constants.o stel_kinds.o 
diagnostics.o: indices.o stel_kinds.o variables.o
indices.o: variables.o
mmc.o: cyclicr_MPI_real.o variables.o
populateMatrix.o: indices.o sparsify.o stel_kinds.o variables.o
populateRHS.o: indices.o stel_kinds.o variables.o
stel_constants.o: stel_kinds.o
stel_kinds.o: 
uniformDiffMatrices.o:	
variables.o: stel_kinds.o

test_make:
	@echo FC is $(FC)
	@echo FLINKER is $(FLINKER)

