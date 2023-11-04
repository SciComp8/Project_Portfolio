CPU set up:
  - Single-threaded process requires 1 thread on any node: --ntasks=1 --cpus=per-task=1
  - Multi-threaded process requires 32 threads on any node: --ntasks=1 --cpus-per-task=32
  - Multi-threaded process requires 32 threads on a node with exclusivity (use only if undoubtedly needed): --ntasks=1 --cpus-per-task=32 --exclusive
  - 8 concurrent but isolated tasks, each task with 1 thread (no multithreading): --ntasks=8 --nodes=1
  - 8 concurrent but isolated tasks, each task with 32 threads (multithreading only within each single task): --ntasks=8 --nodes=1 --cpus-per-task=32
