CPU set up:
1. Single-threaded process requires 1 thread on any node: --ntasks=1 --cpus=per-task=1
2. Multi-threaded process requires 32 threads on any node: --ntasks=1 --cpus-per-task=32
3. Multi-threaded process requires 32 threads on a node with exclusivity (use only if undoubtedly needed): --ntasks=1 --cpus-per-task=32 --exclusive
4. 8 concurrent but isolated tasks, each task with 1 thread (no multithreading): --ntasks=8 --nodes=1
5. 8 concurrent but isolated tasks, each task with 32 threads (multithreading only within each single task): --ntasks=8 --nodes=1 --cpus-per-task=32
