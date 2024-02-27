# Day2 exercises

## Profiling

### a. Investigate the performance of the matmult.py script

- I used cProfile to check this script
- this script does not contain any functions, but only loops so there is nothing really to check with line_profiler
- overall speed: 3.453 seconds

### b. Investigate the performance of the euler72.py script

- I used both cProfile and line_profiler to check this script
- line_profiler suggests that line 52 take most time and I would try to optimize this one first

### c. Improve the performance of the matmult.py script

- chatGPT suggested to use numpy to improve the performance
- cProfiler: 0.843 seconds vs 3.453 seconds