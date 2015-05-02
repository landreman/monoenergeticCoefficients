EXTRA_COMPILE_FLAGS =
EXTRA_LINK_FLAGS = 

#FLINKER = ftn

TARGET = mmc

OBJ_FILES = computeB.o createGrids.o diagnostics.o indices.o mmc.o populateMatrix.o populateRHS.o preallocateMatrix.o readInput.o sparsify.o uniformDiffMatrices.o variables.o

.PHONY: all clean

all: $(TARGET)

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

$(TARGET): $(OBJ_FILES)
	$(FLINKER) -o $(TARGET) $(OBJ_FILES) $(EXTRA_LINK_FLAGS) $(PETSC_LIB)

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

clean::
	rm -f *.o *.mod *.a *.dat *~ $(TARGET)

# Dependencies:
computeB.o: variables.o
createGrids.o: variables.o
diagnostics.o: indices.o variables.o
indices.o: variables.o
mmc.o: variables.o
populateMatrix.o: indices.o sparsify.o variables.o
populateRHS.o: indices.o variables.o
preallocateMatrix.o: variables.o
readInput.o: variables.o
sparsify.o:
uniformDiffMatrices.o:	
variables.o:

test_make:
	@echo FC is $(FC)
	@echo FLINKER is $(FLINKER)

