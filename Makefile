ROOTCXXFLAGS := $(shell root-config --cflags --glibs --ldflags)
LIBS := -lMinuit -lRooFit -lRooFitCore

.SECONDEXPANSION:
.PHONY: all clean


analysis: $$@.cpp
	$(CXX) -g -O3 $< -o $@.exe $(ROOTCXXFLAGS)

temp_analysis: $$@.cpp
	$(CXX) -g $< -o $@.exe $(ROOTCXXFLAGS) $(LIBS)


all: analysis temp_analysis

clean: 
	@rm -f *.exe 