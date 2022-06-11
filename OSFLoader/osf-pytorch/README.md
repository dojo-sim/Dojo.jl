# Learning Object-Centric Neural Scattering Functions for Free-Viewpoint Relighting and Scene Composition
Learning Object-Centric Neural Scattering Functions for Free-Viewpoint Relighting and Scene Composition.
Michelle Guo, Koven Yu, Alireza Fathi, Yen-Yu Chang, Eric Ryan Chan, Ruohan Gao, Thomas Funkhouser, Jiajun Wu ([arXiv](https://arxiv.org/abs/2012.08503) 2020)

Pytorch implementation of learning object-centric neural scattering functions (OSF) based on the [NeRF](https://github.com/yenchenlin/nerf-pytorch) codebase.


## Setup

```
git clone https://github.com/yuyuchang/osf-pytorch.git
cd osf-pytorch
pip install -r requirements.txt
bash download_data.sh #### TODO
```

## Querying a density function from a trained model
A demo for extracting density from a trained model can be found by running:
```python
python extract_density.py --config=configs/osf/bunny_trans/bunny_trans.txt
```

## Training OSFs

```
# Train OSFs on checker with black background
python run_osf.py --config=configs/osf/checkers/rand_light/checkers_black.txt

# Train OSFs on objects from the ObjectFolder dataset
python run_osf.py --config=configs/osf/48/48.txt
python run_osf.py --config=configs/osf/71/71.txt

# Train OSFs on translucent objects
python run_osf.py --config=configs/osf/bunny_trans/bunny_trans.txt
python run_osf.py --config=configs/osf/bunny_trans/cup_trans.txt

# Train OSFs on real objects captured by a cellphone
python run_osf.py --config=configs/osf/bluesoap_first20/bluesoap_first20.txt
python run_osf.py --config=configs/osf/halfsoap_first20/halfsoap_first20.txt
```

We provide backwards compatibility with training NeRFs. The following 
configuration files train NeRFs on various scenes.

```
python run_osf.py --config=configs/nerf_synthetic/bunny_trans/bunny_trans.txt
python run_osf.py --config=configs/nerf_synthetic/cup_trans/cup_trans.txt
```

## Testing OSFs
```
python run_osf.py --config=configs/nerf_synthetic/bunny_trans/bunny_trans.txt --render_only --render_test
```

# TODO: Composing OSFs

After training individual OSFs, you can compose them together into arbitrary
scene arrangements at test time. This example composes the checkers background
and the ObjectFolder objects together into a scene.

```
python run_osf.py --config=configs/checkers_comp_scene_1_0.txt
```

## Configuration Files

We maintain backwards compatibility with the structure of original NeRF 
configuration files. If you'd like to turn on various flags that OSF uses, you 
need to specify them for each object in `object_params` (list of object 
configurations), like the following example. See the `configs` subdirectory for
more examples of configurations to use.

```
object_params: [
    {
        exp_dir: Path to the experiment directory containing a pretrained model for this object,
        intersect_bbox: Whether to use ray-bbox intersections to compute near/far sampling bounds,
        box_center: [-0.05, 0, 0.015],
        box_dims: [0.25, 0.4, 0.2],
        translation: [0.1, 0.15, 0.0875],
        rotation: [0, 0, 90],
        translation_delta: [-0.01, 0, 0],
        lightdirs_method: <lightdirs_method>
    },
    ...
]
```

Options for `lightdirs_method` include:

- `viewdirs`: Use viewing direction of each ray as the incoming light direction
for each sample along the ray. This is useful for datasets where your light
coincides with the camera (e.g., using a phone + flash), such as in the NRF
dataset.

TODO: Do you still need `args.use_lightdirs`? What if different objects have different `use_lightdirs` and `lightdirs_method`?

For composition, set the data directory to be the dataset that you want to use
scene parameters for (?). For instance, we use the XXX(Koven's compose scenes) dataset as the
datadir for `checkers_comp_scene_1_0` configuration file.

Translation / rotation delta is for rendering moving objects.

## TODO: Datasets

Datasets are downloaded to your `./data` folder, and are organized under folders
with the same name as the `dataset_type` used to load the datasets in 
`run_osf.py`.

## Citation

```
@article{guo2020osf,
    title={Object-Centric Neural Scene Rendering},
    author={Guo, Michelle and Fathi, Alireza and Wu, Jiajun and Funkhouser, Thomas},
    journal={arXiv preprint arXiv:2012.08503},
    year={2020}
}
```

We especially thank the authors of the following papers for sharing the NRF dataset with us. If you use the NRF dataset, please cite the following papers:
```
@inproceedings{bi2020deep,
  title={Deep reflectance volumes: Relightable reconstructions from multi-view photometric images},
  author={Bi, Sai and Xu, Zexiang and Sunkavalli, Kalyan and Ha{\v{s}}an, Milo{\v{s}} and Hold-Geoffroy, Yannick and Kriegman, David and Ramamoorthi, Ravi},
  booktitle={Computer Vision--ECCV 2020: 16th European Conference, Glasgow, UK, August 23--28, 2020, Proceedings, Part III 16},
  pages={294--311},
  year={2020},
  organization={Springer}
}

@article{bi2020neural,
  title={Neural reflectance fields for appearance acquisition},
  author={Bi, Sai and Xu, Zexiang and Srinivasan, Pratul and Mildenhall, Ben and Sunkavalli, Kalyan and Ha{\v{s}}an, Milo{\v{s}} and Hold-Geoffroy, Yannick and Kriegman, David and Ramamoorthi, Ravi},
  journal={arXiv preprint arXiv:2008.03824},
  year={2020}
}
```

## Rendering an OSF

Download pretrained model by running

```
download_example_weights.sh
```

Then, run: 

```
python run_nerf.py --config configs/nrf/pony_fwd_render_test.txt
```

You can use `render_only` flag to render out a camera path.
You can use `render_only` and `render_test` to render the test set.

## Training an OSF

To train a nerf on the lego dataset, run
```
python run_nerf.py --config configs/lego.txt
```

To train a light-conditioned nerf:
```
python run_nerf.py --config configs/lego_lightdir.txt
```
(todo: actually feed in light dirs instead of views as light dirs)

Add the flag `--expname_ts` to automatically append `<date>_<time>` to the experiment
name.

## Paper Experiments

Reproducing experiments requires knowing the config file to run, which points to the
appropriate dataset. However, we also need to know how the dataset was generated. Given
a dataset path, there needs to be a mapping to the bash script used to generate the 
dataset. The alternative is to look through every single script and look for the script
that has the dataset name in its args.

## Google Research code

https://github.com/google-research/google-research/blob/master/osf
