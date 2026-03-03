CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra
HTSLIB_CFLAGS := $(shell command -v pkg-config >/dev/null 2>&1 && pkg-config --cflags htslib)
HTSLIB_LIBS := $(shell command -v pkg-config >/dev/null 2>&1 && pkg-config --libs htslib)
CPPFLAGS = -Iinclude $(HTSLIB_CFLAGS)
TARGET = PanelRecon
SRC = PanelRecon.cpp src/PanelIndex.cpp src/PanelFind.cpp src/utils.cpp src/FastqReader.cpp src/RefFasta.cpp
HEADERS = $(wildcard include/*.hpp) $(wildcard include/*.h)

ifeq ($(strip $(HTSLIB_LIBS)),)
LDLIBS = -lhts -lz
else
LDLIBS = $(HTSLIB_LIBS) -lz
endif

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDLIBS)

clean:
	rm -f $(TARGET) PanelRecon_index PanelRecon_find PanelRecon
