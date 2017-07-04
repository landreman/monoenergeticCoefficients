#FLINKER = ftn

TARGET = mmc

OBJ_FILES = apply_preconditioner.o computeB.o createGrids.o diagnostics.o indices.o mmc.o populateMatrix.o populateRHS.o preallocateMatrix.o readInput.o sparsify.o uniformDiffMatrices.o variables.o

.PHONY: all clean

all: $(TARGET)

PETSC_DIR=/opt/local/lib/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

EXTRA_COMPILE_FLAGS = -I${PETSC_DIR}/include -ffree-line-length-none -g -O0
EXTRA_LINK_FLAGS = 

$(TARGET): $(OBJ_FILES)
	$(FLINKER) -o $(TARGET) $(OBJ_FILES) $(EXTRA_LINK_FLAGS) $(PETSC_LIB)

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

clean::
	rm -f *.o *.mod *.a *.dat *~ $(TARGET)

# Dependencies:
apply_preconditioner.o: indices.o variables.o
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
	@echo PETSC_DIR is $(PETSC_DIR)
	@echo PETSC_ARCH is $(PETSC_ARCH)
	@echo FC is $(FC)
	@echo FLINKER is $(FLINKER)

