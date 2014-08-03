.PHONY: clean
SAMTOOLDIR=../libs/samtools-0.1.19/
BOOST=../libs/boost_1_45_0_subset/
GTESTDIR=libs/gtest-1.7.0
GTEST=$(GTESTDIR)/include

SRCDIR=../src
INCLUDEDIR=../include
OBJDIR=obj
BINDIR=bin

CXX=g++
CPPFLAGS=#-DLOGDEBUG -DDEBUGREADS #-DNDEBUG -DMMTEST
CXXFLAGS=-I$(SAMTOOLDIR) -I$(BOOST) -I$(GTEST) -I$(INCLUDEDIR) -I$(SRCDIR) -Wno-deprecated -O3
LDFLAGS=-L$(SAMTOOLDIR) -L$(GTESTDIR) -lbam -lz -lpthread -lgtest

SOURCES := $(shell find $(SRCDIR) -name "*.cpp" ! -name "main.cpp")
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

all: pre-build main-build

pre-build:
	mkdir -p $(OBJDIR) $(BINDIR) $(OBJDIR)/haplotype_phaser

main-build: hapmuc_test

dependencies:
	cd libs;make

test: clean all
	bin/hapmuc_test data

hapmuc_test: $(OBJECTS) obj/main.o
	$(CXX) -o bin/$@ $(CXXFLAGS) $(CPPFLAGS) $(OBJECTS) obj/main.o $(LDFLAGS)

obj/main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR)/*
