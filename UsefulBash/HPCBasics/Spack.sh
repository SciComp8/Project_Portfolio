ssh username@aphrodite.med.cornell.edu
ssh curie.pbtech

# Show the information of installed packages
spack find 
spack find python # x spack find py
spack list
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

# Load/unload the python package
spack load python@3.9.15%gcc@8.2.0
which python
# /pbtech_mounts/software/spack/centos7/opt/spack/linux-centos7-x86_64/gcc-8.2.0/python-3.9.15-t253edbrqg3xsv3olujo5utfeb6ctade/bin/python
python -V
python help

spack unload python@3.9.15%gcc@8.2.0
