
#--------------------------------------General------------------------------------------
CXX      := g++  
#SRCDIR  := .
BUILDDIR := build
TARGET  := hall_conductivity


#--------------------------------------Sources and header files------------------------------------------
SRCEXT  := cpp
#SOURCES := $(shell ls -1 $(SRCDIR)/*.$(SRCEXT))
SOURCES := hall_conductivity.cpp
#OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(patsubst %,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

#--------------------------------------Compiler settings------------------------------------------

CFLAGS +=  -fopenmp  -std=c++11  #General compiler flags
LFLAGS += -O4 -fopenmp 
RUNFLAGS := -march=native   -O4   # Compiler flags for quick compile and run
DBFLAGS := -Wall -Wextra -g3  # Compiler flags for debugging
PROFFLAGS := -O2 -pg  # Compiler flags for profiling
OPTIMIZEFLAGS := -O4  # Compiler flags for optimal speed
LIB := -llapack -lblas # Specify Libraries
INC :=  -I /usr/include/eigen3 -I /usr/include # Additional include paths


#--------------------------------------Targets ------------------------------------------
$(TARGET): $(OBJECTS)
	@echo " Linking..."
#	@mkdir -p bin
	@echo " $(CXX) $(LFLAGS)  $^ -o $(TARGET) $(LIB)"; $(CXX) $(LFLAGS) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: %.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CFLAGS) $(INC) -c  $@ $< " ; $(CXX) $(CFLAGS) $(INC) -c -o  $@ $< 

run:    CFLAGS += $(RUNFLAGS)
run:    $(TARGET)

debug:  CFLAGS += $(DBFLAGS)
debug:  $(TARGET)

prof:   CFLAGS += $(PROFFLAGS)
prof: 	LFLAGS += -pg
prof:   $(TARGET)

optimize: CFLAGS += $(OPTIMIZEFLAGS)
optimize: $(SOURCES) $(HEADERS)
	$(FC) $(FFLAGS) $(INC) $(SOURCES) -o $(TARGET) $(LIB)

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean
