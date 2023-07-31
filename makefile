TARGETS = summary_analysis
FLAGS = -g0 -O3 -D__USE_XOPEN2K8

DEPEND = include/*.h

all: $(TARGETS)

info:
	@echo "ROOTSYS: $(ROOTSYS)"
	@echo "targets: $(TARGETS)"
	@echo "flags: $(FLAGS)"

clean:
	rm -rf MC

$(TARGETS) : %: %.cc
	g++ -std=c++11 -Wall -Woverloaded-virtual -Wconversion $(FLAGS) -pthread -m64 -I/usr/include/root -I$(ROOTSYS)/include -L$(ROOTSYS)/lib `root-config --libs` -lMinuit -I. -ldl -pthread $< -o $@
