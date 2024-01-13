# Show the info of current jobs
squeue --help
squeue -u username -l
squeue -u username --states=COMPLETED
squeue -u username --states=RUNNING
squeue -u username --states=PENDING
squeue -u username --states=SUSPENDED
squeue -u username -p=panda

# Show the info of past jobs
sacct 
sacct -j 11114265
sacct -r panda

# Cancel/pause/resume jobs
scancel 11114268
scontrol hold 11114269
scontrol release 11114269
