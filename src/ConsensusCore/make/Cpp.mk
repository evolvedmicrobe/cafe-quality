include make/Defs.mk

CONFIG_MAKE  = make/Config.mk
ifneq ("$(wildcard $(CONFIG_MAKE))","")
include $(CONFIG_MAKE)
else
$(error Run ./configure first!)
endif


CXX_OBJS        := $(addprefix $(OBJDIR)/, $(CXX_SRCS:.cpp=.o))
SWIG_SRCS       := $(shell find src/SWIG/ -name "*.i")

$(CXX_LIB): $(OBJDIR) $(CXX_OBJS)
	ar crs $(CXX_LIB) $(CXX_OBJS)
	touch $(CXX_LIB)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(CXX_OBJS): $(OBJDIR)/%.o: %.cpp
	$(CXX_STRICT) $(COVERAGE) -c $< -o $@
	$(CXX) -MT $(OBJDIR)/$*.o -MM $< > $(OBJDIR)/$*.d

-include $(CXX_OBJS:.o=.d)

tests: $(CXX_LIB)
	make -f make/Tests.mk run-tests

test: tests
check: tests

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)
