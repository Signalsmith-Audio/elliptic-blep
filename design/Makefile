main: out/main
	@cd out && ./main

out/coeffs.h: design.py
	mkdir -p out
	python3 design.py | tee out/coeffs.h
	
out/main: out/coeffs.h *.cpp ../*.h
	g++ -std=c++11 -Ofast \
		-Wall -Wextra -Wfatal-errors -Wpedantic -pedantic-errors \
		main.cpp -o out/main

clean:
	rm -rf out
