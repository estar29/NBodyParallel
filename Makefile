nbody_par: nbody-parallel.cpp
	   g++ nbody-parallel.cpp -o nbody_par

# Test case with a simple 3 particle system.
test.out: nbody_par
	       ./nbody_par 3 1 100 5 > test.out

test.pdf: test.out
	  python3 plot.py test.out test.pdf 1000
