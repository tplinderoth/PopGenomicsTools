CXX = g++
CXXFLAGS = -O3 -Wall
BIN = fstWindow hetWindow ihsWindow dxyWindow

all: $(BIN)

fstWindow: fstWindow.cpp
	$(CXX) $(CXXFLAGS) fstWindow.cpp -o fstWindow
hetWindow: hetWindow.cpp
	$(CXX) $(CXXFLAGS) hetWindow.cpp -o hetWindow
ihsWindow: ihsWindow.cpp
	$(CXX) $(CXXFLAGS) ihsWindow.cpp -o ihsWindow
dxyWindow: dxyWindow.cpp
	$(CXX) $(CXXFLAGS) dxyWindow.cpp -o dxyWindow -lz -lboost_iostreams

clean:
	rm -f $(BIN) *.o *.d

.PHONY: clean all
