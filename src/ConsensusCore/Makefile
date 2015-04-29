#
# Basic makefile. Use this for building the C++ library and the SWIG
# bindings.
#
# Notes:
#  - can compile with clang++ using:
#      % make GXX=clang++
#  - 32-bit build:
#      % make MACHINE=-m32
#  - debug build:
#      % make DEBUG=1
#
include make/Defs.mk

all: lib

lib:
	${MAKE} -f make/Cpp.mk

#
# SWIG targets
#

csharp: lib
	make -f make/CSharp.mk

	
before-xbuild: csharp
	-cp build/CSharp/libConsensusCore.so ../../lib/libConsensusCore.so

perftest: lib
	clang++ -O3 -DNDEBUG -g -msse3 $(INCLUDES) -fno-omit-frame-pointer -std=c++11  $(OBJDIR)/libConsensusCore.a src/C++/Quiver/MatrixTester.cpp -o $(OBJDIR)/perftest
#
# Clean targets
#

clean-csharp:
	-rm -rf $(CSHARP_BUILD_DIR)
	-rm -rf bin obj

clean-cxx:
	-rm -rf $(BUILD_ROOT)/C++

clean: clean-cxx clean-csharp
	-rm -rf $(BUILD_ROOT)

#
# Test targets
#

test-csharp:
	@make -f make/CSharp.mk test-csharp

test: lib
	@make -f make/Tests.mk

mattest: lib
	@make -f make/Cpp.mk mattest

check: test
tests: test


#
# Targets used by PBI internal build
#
.PHONY: all lib clean-cxx clean test tests check  \
	csharp clean-csharp \
	test-csharp  
