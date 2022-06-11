import glob
dir_list = list(glob.glob("/viscam/u/rhgao/datasets/ObjectFiles/modelhaven/*/VisionNet_osf_data/"))

#f = open("YCB_render_test.sh", "w")
#f = open("GoogleScannedObjects_render_test.sh", "w")
#f = open("modelhaven_render_val.sh", "w")

for dir in dir_list:
    obj_name = dir.split("/")[-3]
    f = open("{}.txt".format(obj_name), "w")
    f.write("expname: {}\n".format(obj_name))
    f.write("basedir: ./logs/modelhaven\n")
    f.write("datadir: {}\n".format(dir[:-1]))
    f.write("dataset_type: osf\n\n")
    f.write("no_batching: True\n\n")
    f.write("use_viewdirs: True\n")
    f.write("use_lightdirs: True\n")
    f.write("white_bkgd: True\n")
    f.write("lrate: 0.0005\n\n")
    f.write("N_samples: 64\n")
    f.write("N_importance: 128\n")
    f.write("N_rand: 1024\n\n")
    f.write("half_res: False\n\n")
    f.write("i_video: 10000000000000000\n")
    f.write("i_testset: 10000000000000000\n")
    f.write("scaled_sigmoid: True\n")
    f.close()


#for dir in dir_list:
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}google_16k/textured.obj --mtl_file {}google_16k/textured.mtl --results_path {}google_16k/VisionNet_data/train\n".format(dir, dir, dir))
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}google_16k/textured.obj --mtl_file {}google_16k/textured.mtl --results_path {}google_16k/VisionNet_data/val\n".format(dir, dir, dir))
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}google_16k/textured.obj --mtl_file {}google_16k/textured.mtl --results_path {}google_16k/VisionNet_data/test\n".format(dir, dir, dir))
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}meshes/model.obj --mtl_file {}meshes/model.mtl --results_path {}VisionNet_data/train\n".format(dir, dir, dir))
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}meshes/model.obj --mtl_file {}meshes/model.mtl --results_path {}VisionNet_data/val\n".format(dir, dir, dir))
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}meshes/model.obj --mtl_file {}meshes/model.mtl --results_path {}VisionNet_data/test\n".format(dir, dir, dir))
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}*.obj --mtl_file {}*.mtl --results_path {}VisionNet_data/train\n".format(dir, dir, dir))
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}*.obj --mtl_file {}*.mtl --results_path {}VisionNet_data/val\n".format(dir, dir, dir))
    #f.write("blender -b -noaudio -P gpu.py -E CYCLES --python 360_view.py -- --mesh_file {}*.obj --mtl_file {}*.mtl --results_path {}VisionNet_data/test\n".format(dir, dir, dir))
