CC=g++
CFLAGS=-c -Wall -std=c++11
LDFLAGS=
SOURCES=main.cpp prognosis/Runge_good.cpp prognosis/GRK8.cpp prognosis/DataConverter.cpp prognosis/Influence/AdditionalFunction.cpp prognosis/Influence/InfluenceAtmosphere.cpp prognosis/Influence/InfluenceEarthRotation.cpp prognosis/Influence/InfluenceEGM96.cpp prognosis/Influence/InfluenceForce.cpp prognosis/Influence/InfluenceNutationEarth.cpp prognosis/Influence/InfluencePlanet.cpp prognosis/Influence/InfluencePoleEarth.cpp prognosis/Influence/InfluenceSun.cpp prognosis/Influence/Influence_Test.cpp prognosis/Influence/InfluenceTime.cpp prognosis/Influence/InfluiencePrecessionEarth.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o prognosis/*.o prognosis/Influence/*.o out
