PROGRAM	= asymcorrlda
CC	= gcc
CFLAGS	= -O3
SRCS	= asymcorrlda.c learn.c likelihood.c hyper.c writer.c feature.c imatrix.c dmatrix.c util.c
OBJS	= $(SRCS:.c=.o)
HEADERS	= $(SRCS:.c=.h)
LDFLAGS	= -lm -lgsl -lgslcblas -L/usr/local/lib/
INCLUDE	= -I/usr/local/include/
VERSION	= 0.1
PKGNAME	= asymcorrldacgs-$(VERSION)
DISTDIR	= ../dist
DISTFILES	= $(SRCS) $(HEADERS) Makefile

all: depend $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

depend:
	@$(CC) $(INCLUDE) -MM $(SRCS) > .depend
clean:
	@rm -f .depend $(OBJS)
pkg:
	@[ -d $(PKGNAME) ] || mkdir $(PKGNAME)
	@cp -p $(DISTFILES) $(PKGNAME)
	@tar czvf $(DISTDIR)/$(PKGNAME).tar.gz $(PKGNAME)
	@rm -r $(PKGNAME)

-include .depend
