#                                                  # solaris/cc
#ELF_CC      = cc 
#ELF_CFLAGS  = -Xa -v -G -Kpic

#                                                  # linux/gcc
ELF_CC      = gcc
ELF_CFLAGS  =  -O2 -W -Wall -shared -fpic

U = new_pdf

SLANG_INC = $(CXC)/$(ARCH)/include
SLANG_LIB = -L/nfs/cxc/a1/i686/lib -lslang
LIBS = $(SLANG_LIB) -Wl,-R/nfs/cxc/a1/i686/lib

INCS = -I../src -I$(SLANG_INC)

all:	$(U).so

$(U).so:	$(U).c
	$(ELF_CC) $(ELF_CFLAGS) $(INCS) $(LIBS) $(U).c -o $(U).so

clean:
	rm -f *.o *.so .##*
	
distclean:
	make clean
	
	
