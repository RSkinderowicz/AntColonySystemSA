CXX = g++

# Change the following to debug to build a debug version
mode = release

ifeq ($(mode),release)
	CXXFLAGS = -pipe -std=c++11 -Wall -pedantic -DNDEBUG -O3 -mtune=native -march=native 
else
mode = debug
	CXXFLAGS = -g -pipe -std=c++11 -Wall -pedantic -O0 -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self  -Wmissing-include-dirs -Woverloaded-virtual -Wredundant-decls  -Wsign-conversion -Wsign-promo  -Wstrict-overflow=5 -Wundef -Werror -Wno-unused
endif

LDFLAGS =
BUILDDIR = obj
TARGET = acs
SRCDIR = src
SOURCES = main.cpp\
		  tsp_instance.cpp\
		  standard_pheromone_memory.cpp\
		  rng.cpp\
		  lru_pheromone_memory.cpp\
		  docopt.cpp\
		  k_opt.cpp\
		  sop_instance.cpp\
		  ant_sop.cpp\
		  sop_local_search.cpp\
		  cmd_line_args.cpp\
		  log.cpp\
		  stop_condition.cpp\
		  fs_utils.cpp


OBJS = $(SOURCES:.cpp=.o)

$(info $$OBJS is [${OBJS}])

OUT_OBJS = $(addprefix $(BUILDDIR)/,$(OBJS))

DEPS = $(SOURCES:%.cpp=$(BUILDDIR)/%.depends)

$(warning $(DEPS))

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OUT_OBJS)
	$(CXX) $(CXXFLAGS) $(OUT_OBJS) $(LDFLAGS) -o $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR)/%.depends: $(SRCDIR)/%.cpp
	@mkdir -p depends
	$(CXX) -MF"$@" -MG -MM -MP  -MT"$(<F:%.cpp=$(BUILDDIR)/%.o)" $(CXXFLAGS) $< > $@

clean:
	rm -f $(OUT_OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
