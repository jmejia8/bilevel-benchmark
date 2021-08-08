# Bilevel Optimization Test Problems

Test Function Suite for Bilevel Optimization implemented in C/C++.

Wrappers available in [Julia](https://github.com/jmejia8/BilevelBenchmark.jl) and [MATLAB](https://github.com/jmejia8/bilevel-benchmark-matlab). 

## Benchmarks

### Single-objective Bilevel Optimization

- PMM: Five test function with pseudo-feasible solutions.
- SMD: Twelve test functions, constrained and unconstrained problems.
- TP: Ten test function from classical literature.

## Build Project

If you are using a Linux-based operative system, you can build using the make command.

```
make
```

That command will generate the <ins>shared libraries</ins> corresponding  to each test suite:
`pmm.so`, `smd.so`, `tp.so` . That can be useful for users working on Python, R, Julia, MATLAB, etc. 


Also you can run the following command (on linux) to build a shared library including all functions (`make` also create this library).

```
gcc blb18_op.c -lm -O2 -march=native  -o blb18_op.so
```

## API References

See `test.c` for examples.

### PMM

```C
/*
 * D_ul: upper level dimension
 * D_ll: lower level dimension
 * x: upper level decision vector
 * y: lower level decision vector
 * fnum: function number, i.e., PMM<fnum>
 * Save in F, G the function and constraints values.
*/
void PMM_leader(int D_ul, int D_ll, double *x, double *y, double *F, double *G, int fnum)
```

```C
/*
 * D_ul: upper level dimension
 * D_ll: lower level dimension
 * x: upper level decision vector
 * y: lower level decision vector
 * fnum: function number, i.e., PMM<fnum>
 * Save in f, g the function and constraints values.
*/
void PMM_follower(int D_ul, int D_ll, double *x, double *y, double *f, double *g, int fnum)
```


```C
/*
 * D_ul: upper level dimension
 * D_ll: lower level dimension
 * x: upper level decision vector
 * fnum: function number, i.e., PMM<fnum>
 * Save in y the lower level optimal solution, i.e., y in argmin f(x, y)
*/
void PMM_Psi(int D_ul, int D_ll, int k, double *x, double *y, int fnum){
```


## Contributing

Please, be free to send me your PR, issue or any comment about this repo.
