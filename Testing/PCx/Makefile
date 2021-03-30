# Makefile.  Generated from Makefile.in by configure.

export PCx_ARCH = linux
export WSSMP_LIB = 
export NG_LIB = ../Ng-Peyton/cholesky.a
export TARGETDIR = .

export CC ?= gcc
export CFLAGS ?= -Wall -O3 -DDECIMAL_WIDTH=6 
export FC ?= gfortran
export FCFLAGS ?= -Wall -O3 

# default
PCx:: NgPeyton
	@echo "***************************************************";
	@echo "Building PCx with the Ng Peyton Cholesky solver....";
	@echo "***************************************************";
	$(MAKE) $(MFLAGS) -C SRC PCx_NgPeyton

NgPeyton::
	@echo "***************************************************";
	@echo "Building the Ng Peyton Cholesky solver....";
	@echo "***************************************************";
	$(MAKE) $(MFLAGS) -C Ng-Peyton

PCx_wssmp::
	@echo "***************************************************";
	@echo "Building PCx with the IBM WSSMP Cholesky solver....";
	@echo "***************************************************";
	$(MAKE) $(MFLAGS) -C SRC PCx_wssmp

PCx_mysolver::
	@echo "***********************************************";
	@echo "Building PCx with the user-supplied solver     ";
	@echo "Requires file ./SRC/mysolver.c and             ";
	@echo "    library ./mysolver/libmysolver.a           ";
	@echo "***********************************************";
	$(MAKE) $(MFLAGS) -C SRC PCx_mysolver

mex:: NgPeyton
	@echo "***********************************************";
	@echo "The Matlab Interface is only available with    ";
	@echo "    the Ng-Peyton solver.                      ";
	@echo "***********************************************";
	$(MAKE) $(MFLAGS) -C mex

doc::
	-mkdir DOC/_doxygen
	cd SRC; doxygen
	$(MAKE) $(MFLAGS) -C DOC html

clean::
	rm -f PCx
	rm -f *.log
	rm -f Ng-Peyton/*.[ao]
	rm -f SRC/*.[ao]
	rm -f mex/*.mex*
	rm -f mex/*.[ao]
	rm -f F2C/*.[ao]
