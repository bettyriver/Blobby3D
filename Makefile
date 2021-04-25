CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG
LIBS = -ldnest4 -lpthread -lfftw3

default:
	make noexamples -C $(DNEST4_PATH)/DNest4/code
	$(CXX) -I $(DNEST4_PATH) $(CXXFLAGS) -c src/*.cpp
	$(CXX) -pthread -L $(DNEST4_PATH)/DNest4/code -o Blobby3D *.o $(LIBS)
	rm *.o

nolib:
	$(CXX) -I $(DNEST4_PATH) $(CXXFLAGS) -c src/*.cpp
	$(CXX) -pthread -L $(DNEST4_PATH)/DNest4/code -o Blobby3D *.o $(LIBS)
	rm *.o
clean:
	rm -f *.o
	rm -f Blobby3D
	
