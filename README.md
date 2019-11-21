# STARPS
A simulation tool for cultivation of bacterial populations according to a Wright-Fisher model of inheritance. Developed by Petter Lindgren and Jon Ahlinder

## Compilation
To compile the C code, the [Gnu Scientific Library](https://www.gnu.org/software/gsl/) is needed. 
```
gcc -o WrightFisher_call WrightFisher_call.c -Wall -I/usr/include -lm -lgsl -lgslcblas
```

