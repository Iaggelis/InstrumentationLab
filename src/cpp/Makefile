ROOTCXXFLAGS := $(shell root-config --cflags --glibs --ldflags)
LIBS := -lMinuit -lRooFit -lRooFitCore


.SECONDEXPANSION:
.PHONY: all clean


analysis: $$@.cpp
	$(CXX) -g -O3 $< -o $@.exe $(ROOTCXXFLAGS)

temp_analysis: $$@.cpp
	$(CXX) -g $< -o $@.exe $(ROOTCXXFLAGS) $(LIBS)

cleaning_data: $$@.cpp
	$(CXX) -g $< -o $@.exe $(ROOTCXXFLAGS) 

root_convert_mt: $$@.cpp
	$(CXX) -g $< -o $@.exe -pthread $(ROOTCXXFLAGS) 



all: clean analysis temp_analysis cleaning_data root_convert_mt

clean: 
	@rm -f *.exe 