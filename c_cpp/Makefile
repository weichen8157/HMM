.PHONY: all clean
CC=gcc
CFLAGS+=-w
LDFLAGS+=-lm     # link to math library

TRAIN_TARGET=train
TEST_TARGET=test

all: train test

train: train.c hmm.h
	$(CC) $(LDFLAGS) $(CFLAGS) $< -o $@

test: test.c hmm.h
	$(CC) $(LDFLAGS) $(CFLAGS) $< -o $@

clean:
	$(RM) $(TEST_TARGET) $(TRAIN_TARGET)   # type make clean to remove the compiled file
