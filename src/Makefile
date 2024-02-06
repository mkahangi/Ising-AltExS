CXX = g++
TARGET = EoS

OBJDIR = ./build
SRCDIR = .
INCLUDEDIR = .

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
HEADERS = $(wildcard $(SRCDIR)/*.hpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
DEPFILES = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.d,$(SOURCES))

LIBS = -lgcc -lgomp
CXXFLAGS = -I$(INCLUDEDIR) -Wall -Wextra -g3 -O3 -std=gnu++20 -MMD -MP -fopenmp 
LDFLAGS = -g3 -fopenmp
MAKEFLAGS += -j$(nproc)
RM = rm -f

.PHONY: all clean run redo build rerun format
.SECONDARY: $(DEPFILES)

all: $(TARGET)

$(TARGET): $(OBJECTS) | build
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p $(OBJDIR)

clean:
	$(RM) $(TARGET) $(OBJECTS) $(DEPFILES)

run: build $(TARGET)
	./$(TARGET) ../input/parameterfile


redo: clean build
	$(MAKE)

rerun: clean build run

format:
	clang-format -style llvm -i $(SOURCES) $(HEADERS)

# do not include dependency files for these targets
ifneq ($(MAKECMDGOALS),clean, format)
-include $(DEPFILES)
endif
