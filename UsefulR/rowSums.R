# Remove any probes that have any samples showing detection p value >= 0.01
probe.keep <- rowSums(det.p.filter < 0.01) == ncol(norm.set) 
