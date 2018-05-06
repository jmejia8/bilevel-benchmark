# A Benchmark for Bilevel Optimization

Test Function Suite for Bilevel Optimization

## Example

Build test:
```
gcc test.c -lm -Ofast -march=native  -o test && time  ./test
```

If you want to call this test function suit from Python, R, Julia, MATLAB, etc. 

Build using (linux):
```
gcc blb18_op.c -lm -Ofast -march=native  -o blb18_op.so
```
