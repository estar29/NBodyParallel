CXXFLAGS=-O3

nbody: nbody.cpp
	g++ -O3 nbody.cpp -o nbody

solar.out: nbody
	date
	./nbody "planet" 200 5000000 10000 > solar.out # maybe a minutes
	date 

# Test case 1
# Same solar case as the sequential code to check if parallel implementation
# is correct.
solar.out: nbody
	date
	./nbody "planet" 200 5000000 10000 > solar_par.out # maybe a minutes
	date

solar.pdf: solar.out
	date
	python plot.py solar.out solar.pdf 1000
	date

# Test case 2
# 100 particles where time step size = 1, and number of steps = 10000.
hun_part.out: nbody
	date
	./nbody 100 1 10000 100 > hun_part.out
	date

# Test case 3
# 1000 particles where time step size = 1, and number of steps = 10000.
thou_part.out: nbody
	date
	./nbody 1000 1 10000 100 > thou_part.out
	date
