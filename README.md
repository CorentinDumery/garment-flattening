
# GG WP: General Garment Woven Parameterization

<table>
    <tr align="middle" >
          <td width="70%" align="left">
          <i>"If your pair of jeans isn't even isotropic, why should its parameterization be?"</i> </br></br>

Previous UV mapping methods (ABF++, LSCM, ARAP, BFF, SCAF, AutoCuts, OptCuts, ...) are rotation invariant and minimize triangle distortion without any assumption on the material being flattened. These methods are unfit for anisotropic materials which stretch unequally along different axes.  

This repository provides an implementation of the anisotropic parameterization described in <a href="https://igl.ethz.ch/projects/computational-patternmaking/computational-pattern-making-paper.pdf">Computational Pattern Making from 3D Garment Models</a> by Pietroni et al.
It's original intent is to accurately flatten woven textiles, for which the thread structure induces
anisotropy in their ability to stretch.</td>
          <td width="30%"><img style="float: right;" src="images/woven_viz.gif" margin="35px"></td>
    </tr>
</table>

## Testing

```
git clone https://github.com/CorentinDumery/garment-flattening
cd garment-flattening/
git submodule update --init --recursive
mkdir build && cd build
cmake ..
make -j woven_param       #build library
make -j woven_viz         #build test app
./woven_viz [path/to/input/mesh.obj]
```

Then, the main features of our parameterization can be visualized with the following commands:
```
make -j woven_viz && ./woven_viz ../data/jumpsuit_multipose/front1.obj
make -j multiple_poses && ./multiple_poses
make -j reflec_param && ./reflec_param
```

## Features

* Shearing: we allow threads to shear, if this helps reduce stretch on the grain axes.

![teaser](images/both_semispheres.png)

* Vertical alignment: patterns are aligned with the vertical axis in 3D.

![teaser](images/align_viz.png)

* Reflectability: opposite sides or seams are constrainted to be a reflection. This makes sewing significantly easier.

![teaser](images/reflec_illus.png)

* Multiple poses: given multiple target 3D meshes representing different poses, we produce a 2D pattern that best fits all targets.

## Adding to your project

```
git add submodule https://github.com/CorentinDumery/garment-flattening
git submodule update --init --recursive
```

Then, in your `CMakeLists.txt`, add:
```
add_subdirectory([path_to_garment-flattening]/garment-flattening)
...
target_link_libraries([target] woven_param)
```

## Citing

Thank you for reading! If this repository is useful to you, feel free to reach out and/or cite our paper:

```
@article{ComputationalPatternmaking:2022,
  title = {Computational Pattern Making from {3D} Garment Models},
  author = {Nico Pietroni and Corentin Dumery and Raphael Falque and Mark Liu and Teresa Vidal-Calleja and Olga Sorkine-Hornung},
  journal = {ACM Transactions on Graphics},
  volume = {41},
  number = {4},
  pages = {157:1–14},
  year = {2022},
  publisher = {ACM}
}
```

![animals_figure](images/animals.png)

## Acknowledgments

Some of the models used in this repository are adapted from
[Generating Datasets of 3D Garments with Sewing Patterns](https://zenodo.org/record/5267549#.YhepENso_mF) by Maria Korosteleva and Sung-Hee Lee.
