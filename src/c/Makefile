CC = gcc
CFLAGS = -O2 -Wall -Wdouble-promotion -fno-math-errno -I/usr/include/eigen3
SOURCES = navtoolbox.c test.c linalg.c blasmini.c benchmark/benchmark.cpp kalman_takasu.c kalman_udu.c
HEADERS = navtoolbox.h linalg.h blasmini.h benchmark/benchmark.h kalman_takasu.h kalman_udu.h
OUTPUT = navtest
LDFLAGS = -lm -lstdc++

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

asm: $(SOURCES)
	$(CC) $(CFLAGS) -fverbose-asm -S $(SOURCES)

.PHONY: profile prof lto analyze tidy format cppcheck clean
analyze:
	scan-build -v clang -fsanitize=address $(CFLAGS) $(SOURCES) -o $(OUTPUT) $(LDFLAGS)

profile: prof
	./$(OUTPUT)
	gprof $(OUTPUT) gmon.out > analysis.txt

lto: CFLAGS += -flto
lto: LDFLAGS += -flto
lto:
	$(CC) $(CFLAGS) $(SOURCES) -o $(OUTPUT) $(LDFLAGS)

prof: CFLAGS += -pg
prof: LDFLAGS += -pg
prof:
	$(CC) $(CFLAGS) $(SOURCES) -o $(OUTPUT) $(LDFLAGS)

tidy:
	clang-tidy $(SOURCES) --

format: $(SOURCES) $(HEADERS)
	clang-format -i $(SOURCES)
	clang-format -i $(HEADERS)

cppcheck:
	cppcheck --enable=all --inconclusive --error-exitcode=1 $(SOURCES)

clean:
	rm -f $(OUTPUT) *.s gmon.out analysis.txt
