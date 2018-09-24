CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG
LIBS = -ldnest4 -lpthread -lfftw3

# For Artemis
# DNEST_PATH = /home/mvar0005
# BOOSTINC = /usr/include

# For MacBook Air
DNEST_PATH = /Users/mathewvaridel
BOOSTINC = /usr/include

default:
	make noexamples -C $(DNEST_PATH)/DNest4/code
	$(CXX) -I $(DNEST_PATH) $(CXXFLAGS) -c *.cpp -I $(BOOSTINC) $(CXXFLAGS)
	$(CXX) -pthread -L $(DNEST_PATH)/DNest4/code -o Blobby3D *.o $(LIBS)
	rm *.o

nolib:
	$(CXX) -I $(DNEST_PATH) $(CXXFLAGS) -c *.cpp
	$(CXX) -pthread -L $(DNEST_PATH)/DNest4/code -o Blobby3D *.o $(LIBS)
	rm *.o
