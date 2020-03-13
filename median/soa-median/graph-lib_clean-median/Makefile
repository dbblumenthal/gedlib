IDIR = ./include

## LSAPE directory
LSAPE_DIR=../lsape-for-graphlib/include

## Eigen library
EIGEN_DIR=../eigen3
OMP_DIR=/usr/local/lib

CXXFLAGS = -I$(IDIR) -I$(LSAPE_DIR) -I../lsape_sinkhorn -I$(EIGEN_DIR) -Wall  -std=c++11 -g -w  -O3 # -O1 -Werror -fsanitize=address heap-use-after-free.cc

BINDIR = ./bin
TESTDIR = ./test_VLDBJ_final
ODIR = ./obj
SRCDIR = ./src

_DEPS = graph.h  utils.h SymbolicGraph.h GraphEditDistance.h ConstantGraphEditDistance.h Dataset.h MultiGed.h BipartiteGraphEditDistance.h BipartiteGraphEditDistanceMulti.h RandomWalksGraphEditDistance.h RandomWalksGraphEditDistanceMulti.h IPFPGraphEditDistance.h  MultistartRefinementGraphEditDistance.h IPFPZetaGraphEditDistance.h  GNCCPGraphEditDistance.h CMUCostFunction.h CMUGraph.h  CMUDataset.h LetterCostFunction.h LetterGraph.h LetterDataset.h IBDDataset.h IBDCostFunction.h IBDGraph.h IBDMedianLabel.h WebCostFunction.h  WebGraph.h  WebDataset.h MedianGraph.h ConstantMedianLabel.h WebMedianLabel.h CMUMedianLabel.h LetterMedianLabel.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_SRCDEPS = utils.cpp  SymbolicGraph.cpp ConstantGraphEditDistance.cpp RandomWalksGraphEditDistance.cpp RandomWalksGraphEditDistanceMulti.cpp  CMUCostFunction.cpp  CMUGraph.cpp  CMUDataset.cpp  LetterCostFunction.cpp  LetterGraph.cpp LetterDataset.cpp IBDDataset.cpp IBDGraph.cpp WebCostFunction.cpp  WebGraph.cpp  WebDataset.cpp ConstantMedianLabel.cpp CMUMedianLabel.cpp LetterMedianLabel.cpp
DEPS_SRC += $(patsubst %,$(SRCDIR)/%,$(_DEPS_SRC))

_OBJ = utils.o SymbolicGraph.o ConstantGraphEditDistance.o RandomWalksGraphEditDistance.o CMUCostFunction.o CMUGraph.o CMUDataset.o LetterGraph.o LetterCostFunction.o LetterDataset.o IBDDataset.o IBDGraph.o WebCostFunction.o  WebGraph.o  WebDataset.o WebMedianLabel.o ConstantMedianLabel.o CMUMedianLabel.o LetterMedianLabel.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJ_QAP = utils.o QAPLibGraph.o QAPLibCostFunction.o QAPLibDataset.o
OBJ_QAP = $(patsubst %,$(ODIR)/%,$(_OBJ_QAP))

# all: $(BINDIR)/test_GraphEditDistance $(BINDIR)/contestGraphEditDistance
all:$(TESTDIR)/test_graph $(TESTDIR)/chemical-edit-distances $(TESTDIR)/benchmark

debug: CXXFLAGS += -DDEBUG -g
debug: $(TESTDIR)/chemical-edit-distances

with_times: CXXFLAGS += -D PRINT_TIMES
with_times: $(TESTDIR)/chemical-edit-distances

multithread: CXXFLAGS += -fopenmp
multithread: $(TESTDIR)/chemical-edit-distances

#uncomment next line on macOS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
multithread_with_times: CXXFLAGS += -fopenmp -D PRINT_TIMES
endif
ifeq ($(UNAME_S),Darwin)
multithread_with_times: CXXFLAGS += -Xpreprocessor -fopenmp -D PRINT_TIMES -L$(OMP_DIR) -lomp
endif
multithread_with_times: $(TESTDIR)/chemical-edit-distances

qaplib_mt: CXXFLAGS += -fopenmp
qaplib_mt: $(TESTDIR)/QAPlib

optim: CXXFLAGS += -O3
#optim: all


$(TESTDIR)/benchmark: $(DEPS) $(OBJ) $(TESTDIR)/benchmark.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml #-lomp

$(TESTDIR)/chemical-edit-distances: $(TESTDIR)/computeDistances.cpp $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml

$(TESTDIR)/QAPlib: $(TESTDIR)/QAPLib.cpp $(OBJ_QAP)
	$(CXX)  -D PRINT_TIMES -o $@ $^ $(CXXFLAGS) -ltinyxml #-lomp

$(TESTDIR)/test_graph: $(DEPS) $(OBJ) $(TESTDIR)/test_graph.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml #-lomp

$(BINDIR)/%: $(OBJ) $(SRCDIR)/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) #-lomp

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS) $(DEPS_SRC)
	$(CXX) -c -o $@ $< $(CXXFLAGS) #-lomp


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  $(BINDIR)/*
