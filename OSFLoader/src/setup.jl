# Setup described here
# https://github.com/JuliaPy/PyCall.jl#python-virtual-environments

@info "Setting up PyCall"

# Step 1: activate environment .osf_pyenv
#$ source osf_pytorch/.osf_py_env/bin/python
#$ pip install -r requirements.txt

# Step 2:
ENV["PYCALL_JL_RUNTIME_PYTHON"] = "/home/simon/research/repos/osf-pytorch/.osf_pyenv/bin/python"
ENV["PYTHON"] = "/home/simon/research/repos/osf-pytorch/.osf_pyenv/bin/python"

# Step 3:
using Pkg
Pkg.build("PyCall")
using PyCall

# Step 4:
# pushfirst!(pyimport("sys")."path", "")
global OSF_PATH = joinpath("/home/simon/research/repos/osf-pytorch")
pushfirst!(pyimport("sys")."path", OSF_PATH)
