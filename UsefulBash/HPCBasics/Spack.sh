ssh username@aphrodite.med.cornell.edu
ssh curie.pbtech

# Show the information of installed packages
spack find # x
spack find | less -R # Useful
spack --color always find | less -R # Useful
# SPACE: forward | b: backward | g: beginning | G: end | q: quit
spack find python # x spack find py
# -- linux-centos7-x86_64 / gcc@4.8.5 -----------------------------
# python@3.5.0  python@3.8.12
# -- linux-centos7-x86_64 / gcc@8.2.0 -----------------------------
# python@3.8.12  python@3.8.12  python@3.9.12  python@3.9.15
# x86_64
echo $PATH
vim ~/.bashrc
spack_target_find () {
    spack find | sed -r "s/\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]//g" | grep -E "${1}|-- linux" | sed 's/-- linux/\n-- linux/g'
}
spack_target_find py
# -- linux-centos7-zen / gcc@9.2.0 --------------------------------
# py-cycler@0.10.0
# py-cython@0.29.14
# py-kiwisolver@1.1.0
# py-six@1.12.0
# python@3.7.6
spack find -v -l python@3.7.6 # Varied hashes indicate the differences between installations of the same-version python

# Show the information of available packages (packages.spack.io)
spack list
spack list | grep "^py-"
spack list | grep "python"
spack list | grep "bowtie"

# Load/unload the python package
spack load python@3.9.15%gcc@8.2.0 # version number: 3.9.15, compiler: gcc@8.2.0
spack find --loaded # View the loaded packages
which python
# /pbtech_mounts/software/spack/centos7/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-3.9.15-t253edbrqg3xsv3olujo5utfeb6ctade/bin/python
python -V
python help
spack unload python@3.9.15%gcc@8.2.0

spack load python@3.8.12%gcc@8.2.0
# ==> Error: python@3.8.12%gcc@8.2.0 matches multiple packages.
#   Matching packages:
#     rern7pd python@3.8.12%gcc@8.2.0 arch=linux-centos7-haswell
#     4mo5tkb python@3.8.12%gcc@8.2.0 arch=linux-centos7-haswell
#     kbqv5ik python@3.8.12%gcc@8.2.0 arch=linux-centos7-x86_64
#     qefl54h python@3.8.12%gcc@8.2.0 arch=linux-centos7-x86_64
#   Use a more specific spec (e.g., prepend '/' to the hash).
spack load /kbqv5ik
spack unload /kbqv5ik

spack load bowtie2@2.3.5.1%gcc@8.2.0
# ==> Error: bowtie2@2.3.5.1%gcc@8.2.0 matches multiple packages.
#   Matching packages:
#     tlj7iop bowtie2@2.3.5.1%gcc@8.2.0 arch=linux-centos7-broadwell
#     ztcq4ql bowtie2@2.3.5.1%gcc@8.2.0 arch=linux-centos7-sandybridge
#   Use a more specific spec (e.g., prepend '/' to the hash).
spack load /ztcq4ql
bowtie2

spack unload # Unload all packages
spack find --loaded # ==> 0 loaded packages

