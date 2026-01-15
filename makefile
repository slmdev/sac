CXX := g++
CXXFLAGS := -Wall -O3 -std=c++20 -march=native
LDFLAGS  := -static -s

# OS detection
ifeq ($(OS),Windows_NT)
    EXE_EXT := .exe
    RM := del /Q /F
    RMDIR := rmdir /S /Q
    # Function to create directory if it doesn't exist
    MKDIR_P = if not exist $(subst /,\,$(1)) mkdir $(subst /,\,$(1))
else
    EXE_EXT :=
    RM := rm -f
    RMDIR := rm -rf
    MKDIR_P = mkdir -p $(1)
endif

SRC_DIR := src
OBJ_DIR := bin
TARGET  := sac$(EXE_EXT)


SRC_DIRS := src src/common src/file src/libsac src/model src/pred src/opt
SRCS := $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

$(info SRCS = $(SRCS))
$(info OBJS = $(OBJS))

all: $(TARGET)

# Link object files
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@

# Compile .cpp -> .o
$(OBJ_DIR)/%.o: src/%.cpp
	@$(call MKDIR_P,$(dir $@))
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET)
	$(RMDIR) $(OBJ_DIR)

.PHONY: all clean