#FLINKER = ftn

TARGET = mmc

OBJ_FILES = \
	chooseFourierModes.o \
	computeB.o \
	diagnostics.o \
	dmrgrnkMod.o \
	FourierConvolutionMatrixMod.o \
	FourierDifferentiationMatrices.o \
	FourierTransformMod.o \
	indices.o \
	initFourier.o \
	mmc.o \
	populateMatrix.o \
	populateRHS.o \
	preallocateMatrix.o \
	readInput.o \
	sparsify.o \
	variables.o

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
chooseFourierModes.o: dmrgrnkMod.o FourierTransformMod.o variables.o
computeB.o: variables.o
diagnostics.o: FourierConvolutionMatrixMod.o FourierTransformMod.o indices.o variables.o
dmrgrnkMod.o:
FourierConvolutionMatrixMod.o: variables.o
FourierDifferentiationMatrices.o:
FourierTransformMod.o: variables.o
indices.o: variables.o
initFourier.o: variables.o
mmc.o: variables.o
populateMatrix.o: FourierConvolutionMatrixMod.o FourierTransformMod.o indices.o sparsify.o variables.o
populateRHS.o: FourierConvolutionMatrixMod.o FourierTransformMod.o indices.o variables.o
preallocateMatrix.o: variables.o
readInput.o: variables.o
sparsify.o:
variables.o:

test_make:
	@echo FC is $(FC)
	@echo FLINKER is $(FLINKER)

