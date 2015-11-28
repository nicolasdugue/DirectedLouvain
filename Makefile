#!/bin/bash

CC=g++

SRCDIR=src
HEADDIR=include
LIBDIR=obj
BINDIR=bin

CFLAGS= -ansi -O5 -Wall
LDFLAGS= -ansi -lm -Wall
EXEC=community convert hierarchy high_degree inverse reverse overlap

SRC= $(wildcard $(SRCDIR)/*.cpp)
OBJ1= $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/graph_binary.o) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/community.o)
OBJ2= $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/graph.o)
OBJ3 = $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/graph.o) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/graph_binary.o) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/community.o)

all: $(EXEC)

community : $(OBJ1) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/main_community.o)
	$(CC) -g -o $(BINDIR)/$@ $^ $(LDFLAGS)

convert : $(OBJ2) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/main_convert.o)
	$(CC) -g -o $(BINDIR)/$@ $^ $(LDFLAGS)

hierarchy : $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/main_hierarchy.o)
	$(CC) -g -o $(BINDIR)/$@ $^ $(LDFLAGS)

high_degree : $(OBJ3) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/main_high_degree.o)
	$(CC) -g -o $(BINDIR)/$@ $^ $(LDFLAGS)

inverse : $(OBJ3) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/main_compare.o)
	$(CC) -g -o $(BINDIR)/$@ $^ $(LDFLAGS)

reverse : $(OBJ3) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/main_reverse.o)
	$(CC) -g -o $(BINDIR)/$@ $^ $(LDFLAGS)

overlap : $(OBJ3) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/main_overlap.o)
	$(CC) -g -o $(BINDIR)/$@ $^ $(LDFLAGS)

##########################################
# Generic rules
##########################################

$(LIBDIR)/%.o: $(SRCDIR)/%%.cpp $(HEADDIR)/%.h
	$(CC) -g -o $@ -c $< $(CFLAGS)

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) -g -o $@ -c $< $(CFLAGS)

clean:
	rm -f $(LIBDIR)/*.o $(LIBDIR)/*~ $(SRCDIR)/*~ 
