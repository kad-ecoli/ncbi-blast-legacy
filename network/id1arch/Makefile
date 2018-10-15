include $(NCBI)/ncbi.mk


INCLUDES = -I . -I$(NCBI_INCDIR) 
CFLAGS=$(CCOPT) $(INCLUDES) -g #-DDEBUG_ASN_PRINT

LDFLAGS = -L. -L$(NCBI_ALTLIB)  -lncbiid1 -lnetcli -lncbitool -lblastcompadj -lncbiobj -lncbi -lm

CC = gcc -g

BINARIES = idfetch

##
## some things to make
##

all : $(BINARIES)

idfetch : idfetch.o 
	$(CC) -o $@ idfetch.o  $(LDFLAGS)

idfetch.pure : idfetch.o 
	purify $(CC) -o $@ idfetch.o $(LDFLAGS)

test : test.o
	$(CC) -o $@ test.o $(LDFLAGS)

seqidtest : seqidtest.o
	$(CC) -o $@ seqidtest.o $(LDFLAGS)
sample : sample.o
	$(CC) -o $@ sample.o $(LDFLAGS)
gi2hash : gi2hash.o id_hash.o
	$(CC) -o $@ gi2hash.o id_hash.o $(LDFLAGS)
generate : id1gen.c

library : libncbiid1.a

libncbiid1.a : id1gen.o id1arch.o accid1.o
	rm -f $@
	ar -q $@ id1gen.o id1arch.o accid1.o


#id1gen.c : id1.l id1.h all.l id1map.h
#	rm -f id1generr id1genout
#	 $(CODEGEN) -i id1map.h -l id1.l,all.l -o id1gen -d . > id1genout 2>id1generr

#id1.l : id1.asn
#	rm -f id1.l*
#	rm -f id1.h id1stat.h
#	$(ASNTOOL) -m id1.asn -l id1.h
#	rm -f id1.h
#	$(ASNTOOL) -m id1.asn -o id1.h
#	mv id1.l* id1.l

#all.l : 
#	rm -f all.l*
#	$(ASNTOOL) -m $(NCBI_ASNALL) -l all.h
#	cp all.l* all.l
#
#id1.h : id1.asn
#	$(ASNTOOL) -m id1.asn -o id1.h



clean :
	- rm -f *.o
