#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# BR copied from http://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html

## BR definitions.
MYGSLDIR = /usr/local

# define the C compiler to use
CC = gcc

# define any compile-time flags
##CFLAGS = -Wall -g -DREALFLOAT
## nevermind, make the reals doubles!!
CFLAGS = -Wall -g -fopenmp 

# define any directories containing header files other than /usr/include
#
MYINCLUDE = -I$(MYGSLDIR)/include

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -lm -lgsl -lgslcblas

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
MYLIB = -L$(MYGSLDIR)/lib

# define the C source files
# SRCS = cosmoutils.c misc.c readradecz.c xibin.c header.c readnbody.c xi.c xigrid.c
SRCS = header.c cosmoutils.c misc.c angupweight.c readradecz.c readnbody.c xigrid.c xibin.c computesepnbody.c computesepradecz.c countpairs.c output.c xi.c


# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
# define the executable file 
MAIN = xi

all:    $(MAIN)
	@echo  Compiled xi!

OBJS = $(SRCS:.c=.o)

header.o: header.h
misc.o: header.o 
angupweight.o: header.o misc.o 
cosmoutils.o: header.o misc.o 
computesepnbody.o: header.o misc.o 
computesepradecz.o: header.o misc.o  
readradecz.o: header.o misc.o 
readnbody.o: header.o misc.o 
xibin.o: header.o misc.o 
xigrid.o: header.o misc.o 
output.o: header.o misc.o 
countpairs.o: header.o misc.o xigrid.o computesepnbody.o computesepradecz.o angupweight.o
xi.o: header.o misc.o xibin.o readnbody.o readradecz.o angupweight.o countpairs.o output.o



#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)
