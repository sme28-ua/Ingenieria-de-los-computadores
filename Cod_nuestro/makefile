CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11

TARGET = ejecutable_secuencial

all: $(TARGET)

$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o $(TARGET)

run: $(TARGET)
	./$(TARGET) 200 1000

clean:
	rm -f $(TARGET)