CC = gcc
CFLAGS = -O2 -Wall -Wdouble-promotion -fno-math-errno
SOURCES = navtoolbox.c test.c linalg.c
OUTPUT = navtest

all: $(OUTPUT)

$(OUTPUT): $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -o $(OUTPUT) -lm

run_test: $(OUTPUT)
	./$(OUTPUT)
	@RETVAL=$$?; \
	if [ $$RETVAL -ne 0 ]; then \
		echo "Test failed with return value $$RETVAL"; \
		exit 1; \
	else \
		echo "Test passed"; \
	fi

.PHONY: analyze
analyze:
	scan-build $(CC) -fsanitize=address $(CFLAGS) $(SOURCES) -o $(OUTPUT)

clean:
	rm -f $(OUTPUT)
