
for i in range(40):
    with open("slurm_tmp.sh", "r") as f:
        filestring = f.read()
        filestring = filestring.replace("replace_render_start", str(i))
        filestring = filestring.replace("replace_render_end", str(i+1))
    with open("slurm_comp_1_nerf_without_S_{}.sh".format(str(i)), "w") as fout:
        fout.write(filestring)
