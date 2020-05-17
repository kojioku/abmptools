FC            = gfortran
# CFLAGS        = -O4 -Wall -I/usr/local/include
# DEST          = /usr/local/bin
# LDFLAGS       = -L/usr/local/lib
# LIBS          = -lhoge -lm
# OBJS          = main.o sub.o
SRCDIR        = ./abmptools/f90/src
BINDIR        = ./abmptools/f90/bin
# OPT           = -static-libgfortran -static-libgcc
OPT           = -static
OPT2          = -shared -fPIC
SUB           = ./abmptools/f90/src/readifiepiedalib.f90
PROGRAM       =  $(BINDIR)/readifiepiedalib


all: $(PROGRAM)

install: $(PROGRAM)
	python setup.py install

# for f90
$(BINDIR)/readifiepiedalib: $(SRCDIR)/readifiepiedalib.f90
	$(FC) $(SRCDIR)/readifiepiedalib.f90 $(OPT2) -o $(BINDIR)/readifiepiedalib.so

# clean
clean:; rm -f $(PROGRAM)

