CC = mpicc
CFLAGS = -O3 -Wall -Iinclude $(shell pkg-config --cflags PETSc)
LDFLAGS = $(shell pkg-config --libs PETSc) -llapack -lblas -lm

TARGET = test_dlk
SRCS = src/dmdc_lqr.c src/kspdlk.c tests/test_dlk.c
OBJS = $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
