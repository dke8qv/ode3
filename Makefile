P56xxLIBS =  $(PWD)/p56xxlib

ifeq ($(P56xxLIBS),)
$(error P56xxLIBS environment variable is missing, source <path>/p56xxlib/setup.sh)
endif

#ROOTCFLAGS = $(shell root-config --cflags)
#ROOTLIBS   = $(shell root-config --libs)
#ROOTGLIBS  = $(shell root-config --glibs)
#ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
#CXXFLAGS  += $(ROOTCFLAGS) -I$(P56xxLIBS) -Wall -O3
#LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) -Wl,-rpath,$(P56xxLIBS) -L$(P56xxLIBS) -lP56xx
#GXX	   = g++ $(CXXFLAGS)

#SRCS = $(wildcard *.cpp)
#OBJ = $(SRC:.cpp=.o)
#EXES = $(SRCS:%.cpp=%)
#dep = $(OBJ:.o=.d)  # one dependency file for each source

SUBDIRS = p56xxlib src

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS): 
	make -C $@

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean ;\
	done

