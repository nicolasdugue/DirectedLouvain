#!/bin/bash

CC=g++

SRCDIR=src
HEADDIR=include
LIBDIR=obj
BINDIR=bin

CFLAGS := --std=c++11 -Wall -Wextra -pedantic -ggdb -Wno-unused-parameter -Wno-return-type -Wno-variadic-macros
LDFLAGS= -lm
EXEC=community hierarchy

SRC= $(wildcard $(SRCDIR)/*.cpp)
OBJ1= $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/graph.o) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/community.o)

all: $(EXEC)
Debug: CFLAGS += -DDEBUG -g
Debug: LDFLAGS += -DDEBUG -g
Debug: $(EXEC)

community : $(OBJ1) $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/main_community.o)
	[ -d $(BINDIR) ] || mkdir -p $(BINDIR)
	$(CC)  -o $(BINDIR)/$@ $^ $(LDFLAGS)

hierarchy : $(SRC:$(SRCDIR)/%.cpp=$(LIBDIR)/hierarchy.o)
	$(CC)  -o $(BINDIR)/$@ $^ $(LDFLAGS)

##########################################
# Generic rules
##########################################

$(LIBDIR)/%.o: $(SRCDIR)/%%.cpp $(HEADDIR)/%.hpp
	$(CC)  -o $@ -c $< $(CFLAGS)

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp
	[ -d $(LIBDIR) ] || mkdir -p $(LIBDIR)
	$(CC)  -o $@ -c $< $(CFLAGS)

clean:
	rm -f $(BINDIR)/* $(LIBDIR)/*.o $(LIBDIR)/*~ $(SRCDIR)/*~
