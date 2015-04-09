EXTRA_COMPILE_FLAGS =
EXTRA_LINK_FLAGS = 

#FLINKER = ftn

TARGET = mmc

OBJ_FILES = geometry.o mmc.o populateMatrix.o populateRHS.o variables.o

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
geometry.o: variables.o
mmc.o: variables.o
populateMatrix.o: geometry.o variables.o
populateRHS.o: geometry.o variables.o
variables.o:

test_make:
	@echo FC is $(FC)
	@echo FLINKER is $(FLINKER)

