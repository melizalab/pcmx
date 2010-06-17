# Makefile 
#
#
CFLAGS	= -O2 -g2 -Wall
LIBS =  -lm

PROGRAM = pcmx
SRCS = pcmx.c rateconv.c pcmfix.c pcmio.c pcmraw.c pcmseq2.c pcmseq2_read.c pcmseq2_write.c pcmwav.c expidx.c
HDRS = 
OBJS = $(SRCS:.c=.o)

CC = gcc
PREFIX=/usr/local
MODE=655
OWNER=root
GROUP=root
INSTALL=/usr/bin/install -b -D -o ${OWNER} -g ${GROUP}

all: ${PROGRAM}

${PROGRAM} : ${OBJS} 
	$(CC) ${OPT_CFLAGS} -o ${PROGRAM} $(OBJS) ${LIBS}

clean:
	rm -f core $(OBJS) pcmx
