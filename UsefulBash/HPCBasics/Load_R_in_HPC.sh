ssh aphrodite
ssh curie.pbtech
spack find r
# -- linux-centos7-broadwell / gcc@8.2.0 --------------------------
# r@3.6.3  r@4.0.0  r@4.0.2  r@4.0.3

# -- linux-centos7-cascadelake / gcc@9.2.0 ------------------------
# r@3.6.1

# -- linux-centos7-sandybridge / gcc@8.2.0 ------------------------
# r@3.6.1  r@4.0.3  r@4.0.3

# -- linux-centos7-x86_64 / gcc@4.8.5 -----------------------------
# r@4.2.2

# -- linux-centos7-x86_64 / gcc@8.2.0 -----------------------------
# r@4.1.1  r@4.1.1  r@4.1.3  r@4.2.0  r@4.2.2

# -- linux-centos7-zen / gcc@8.2.0 --------------------------------
# r@4.0.3

# -- linux-centos7-zen / gcc@9.2.0 --------------------------------
# r@3.6.1

spack load r@4.2.2%gcc@8.2.0
which R
# /pbtech_mounts/software/spack/centos7/opt/spack/linux-centos7-x86_64/gcc-8.2.0/r-4.2.2-5vatd3o562m3jqy2o5apixkg22gqjepj/bin/R

R --version
# R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"

R
