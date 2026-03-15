######### Written by thchiu 2013/02/03 ##########
CXX      = g++
LINK     = g++
CXXFLAGS = -O3 -Wall -DOUTPUT_FILE
#CXXFLAGS = -O3 -Wall -DPROFILE -DOUTPUT_FILE -DDEBUG
INCPATH  =
LIBS     = 
DEL_FILE = rm -f

####### Files

HEADERS = beta.h buildpathcttable.h count.h data_store.h general_utils.h KLdist.h math_utils.h output.h readdata.h sample_density_generation.h sampling.h SISfunctions.h stl.h timer.h type.h timer_ed520.h tree.h header.h density_tree_forest.h 

SOURCES = main.cpp count.cpp KLdist.cpp sampling.cpp SISfunctions.cpp timer_ed520.cpp tree.cpp density_tree_forest.cpp


OBJECTS = main.o count.o KLdist.o sampling.o SISfunctions.o timer_ed520.o tree.o density_tree_forest.o

TARGET  = gibbs

####### Implicit rules

.SUFFIXES: .o .cpp .cc .cxx .cu

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cu.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

####### Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS) $(INCPATH)

$(OBJECTS): $(HEADERS) Makefile

clean:
	@$(DEL_FILE) $(OBJECTS) $(TARGET)

