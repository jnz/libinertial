CC = gcc
CFLAGS = -O2 -Wall -Wdouble-promotion -fno-math-errno
SOURCES = navtoolbox.c test.c linalg.c
HEADERS = navtoolbox.h linalg.h
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

format: $(SOURCES) $(HEADERS)
	clang-format -i $(SOURCES)
	clang-format -i $(HEADERS)

clean:
	rm -f $(OUTPUT)
