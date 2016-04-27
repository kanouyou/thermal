# Makefile for thermal simulation
# Author: 03.24.2016
# Y. Yang

TARGET  = SimFDM.exe
CXX     = g++
MAINSRC = testfdm.cpp
LIBSRC  = IFdmException.cpp IFdmRadiation.cpp\
          IFdmAluminium.cpp  IFdmCopper.cpp \
          IFdmNbTi.cpp IFdmKapton.cpp \
          IFdmHeatContainer.cpp

SRC     = $(MAINSRC) $(LIBSRC)
OBJS    = $(SRC:.cpp=.o)

CXXLIBS   =
CXXFLAGS  = -O3 -Wall
ROOTLIBS  = `root-config --evelibs`
ROOTFLAGS = `root-config --cflags`

CXXLIBS  += $(ROOTLIBS)
CXXFLAGS += $(ROOTFLAGS)

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXLIBS) $^ -o $@

.o:.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS)

.PHONY: install
install:
	mkdir -p bin
	cp -p $(TARGET) bin

