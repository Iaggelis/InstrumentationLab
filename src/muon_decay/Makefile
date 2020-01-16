ROOTCXXFLAGS := $(shell root-config --cflags --glibs --ldflags)
LIBS := -lMinuit -lRooFit -lRooFitCore


.SECONDEXPANSION:
.PHONY: all clean



: $$@.cpp
	$(CXX) -g $< -o $@.exe $(ROOTCXXFLAGS) $(LIBS)


all:

clean:
	@rm -f *.exe
