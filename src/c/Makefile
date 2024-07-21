CC = gcc
CFLAGS = -O2 -Wall -Wdouble-promotion -fno-math-errno
SOURCES = navtoolbox.c test.c linalg.c
HEADERS = navtoolbox.h linalg.h
OUTPUT = navtest
LDFLAGS = -lm

all: $(OUTPUT)

$(OUTPUT): $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -o $(OUTPUT) $(LDFLAGS)

run_test: $(OUTPUT)
	./$(OUTPUT)
	@RETVAL=$$?; \
	if [ $$RETVAL -ne 0 ]; then \
		echo "Test failed with return value $$RETVAL"; \
		exit 1; \
	else \
		echo "Test passed"; \
	fi

.PHONY: analyze tidy format cppcheck clean
analyze:
	scan-build -v clang -fsanitize=address $(CFLAGS) $(SOURCES) -o $(OUTPUT) $(LDFLAGS)

tidy:
	clang-tidy $(SOURCES) --

format: $(SOURCES) $(HEADERS)
	clang-format -i $(SOURCES)
	clang-format -i $(HEADERS)

cppcheck:
	cppcheck --enable=all --inconclusive --error-exitcode=1 $(SOURCES)

clean:
	rm -f $(OUTPUT)