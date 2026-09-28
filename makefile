# Compilador y opciones por defecto
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11

# Nombre del ejecutable
TARGET = ejecutable_secuencial

all: $(TARGET)

$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o $(TARGET)

# Comando Make Run obligatorio según las especificaciones de la práctica
run: $(TARGET)
	./$(TARGET) 200 1000

clean:
	rm -f $(TARGET)