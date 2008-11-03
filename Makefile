# Makefile 
#
#
OPT_CFLAGS = -O2 -g2 -Wall

#DATAIO_DEPS = -I/usr/local/include
DATAIO_LIBS =  -ldataio

CFLAGS	= ${OPT_CFLAGS} ${DATAIO_DEPS} 
ALL_LIBS = ${DATAIO_LIBS} -lm

PROGRAM = pcmx
SRCS = pcmx.c rateconv.c
HDRS = 
OBJS = $(SRCS:.c=.o)

TVER := $(shell grep BUILD version.h | cut -d' ' -f3)
BUILDVER := $(shell expr 1 + ${TVER})

PREFIX=/usr/local
MODE=655
OWNER=root
GROUP=root
INSTALL=/usr/bin/install -b -D -o ${OWNER} -g ${GROUP}

all: ${PROGRAM}

${PROGRAM} : ${OBJS} 
	$(CC) ${OPT_CFLAGS} -o ${PROGRAM} $(OBJS) ${ALL_LIBS}

clean:
	rm -f core $(OBJS) aplot

spotless: clean
	rm -f .depend

dep depend:
	$(CPP) -M $(CFLAGS) $(SRCS) > .depend

install: ${PROGRAM}
	${INSTALL} -m ${MODE} ${PROGRAM} ${PREFIX}/bin/${PROGRAM}

release:
	-rm -f version.h
	-echo "#define APLOT_BUILD "$(BUILDVER) > version.h
		-echo "#define APLOT_VERSION "`date +%Y%m%d` >> version.h

tar: 
	make spotless
	make release
	tar -C ../ -cplf ../${PROGRAM}_`date +%Y%m%d`.tar ${PROGRAM} --exclude ${PROGRAM}/misc
	bzip2 -9 ../${PROGRAM}_`date +%Y%m%d`.tar
	mv -f ../${PROGRAM}_`date +%Y%m%d`.tar.bz2 /var/www/html/software/

#
# include a dependency file if one exists
#
ifeq (.depend,$(wildcard .depend))
include .depend
endif

