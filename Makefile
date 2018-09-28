SRC_DIR := src
BUILD_DIR := build
TARGET_DIR := bin

GENDATA_SRCS := $(SRC_DIR)/gen_data.cpp $(SRC_DIR)/model.cpp $(SRC_DIR)/reaction_info.cpp
GENDATA_OBJS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(GENDATA_SRCS:.cpp=.o))

IP_SRCS := $(SRC_DIR)/hydrogen_ip.cpp $(SRC_DIR)/compute.cpp $(SRC_DIR)/reaction_info.cpp $(SRC_DIR)/likelihood.cpp $(SRC_DIR)/model.cpp $(SRC_DIR)/qoi.cpp
IP_OBJS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(IP_SRCS:.cpp=.o))

QUESO_DIR := /home/rebecca/bin/queso
ANTIOCH_DIR := /home/rebecca/bin/antioch
COMBUSTION_DIR := /home/rebecca/repos/combustion

INC_PATHS := \
	-I. \
	-I$(QUESO_DIR)/include \
  -I$(ANTIOCH_DIR)/include \
  -I$(COMBUSTION_DIR)/include \

LIBS := \
	-L$(QUESO_DIR)/lib -lqueso \
	-L$(ANTIOCH_DIR)/lib -lantioch \
	-lboost_program_options \
	-lgsl -lgslcblas\

# CXX = mpic++
CXX = mpicxx.openmpi
CXXFLAGS += -O0 -g -Wall -c -std=c++11

default: hydrogen_ip

.SUFFIXES: .o .cpp

all:	 gen_data hydrogen_ip

clean:
	rm -f *~
	rm -f $(BUILD_DIR)/*.o
	rm -f $(TARGET_DIR)/gen_data $(TARGET_DIR)/hydrogen_ip

gen_data: $(GENDATA_OBJS)
	$(CXX) $^ -o $(TARGET_DIR)/gen_data $(LIBS)

hydrogen_ip: $(IP_OBJS)
	$(CXX) $^ -o $(TARGET_DIR)/hydrogen_ip $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(INC_PATHS) $(CXXFLAGS) -o $@ $<
