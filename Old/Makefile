TARGET = a.out

HDRS = include

SRC = src/main.cpp

# Testing

TARGET_TEST = test.out

HDRS_TEST = include \
			test/include

SRC_TEST = test/src/main.cpp \
		   src/integer.cpp \
		   test/src/integerTest.cpp

.PHONY: all main test clean

all: main test

main: $(SRCS)
		$(CXX) -std=c++17 -O2 -fopenmp -Wall -Wextra $(addprefix -I,$(HDRS)) -o $(TARGET) $(CXXFLAGS) $(SRC)

#test: $(SRC_TEST) \
		$(CXX) -std=c++17 -Wall -Wextra $(addprefix -I,$(HDRS_TEST)) -o $(TARGET_TEST) $(CXXFLAGS) $(SRC_TEST)

clean:
	rm -rf $(TARGET) $(TARGET_TEST)