all: prep
	+cmake --build build --config Release

prep:
	+cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

clean:
	rm -rf build

.PHONY: all prep clean

