
# GG WP: General Garment Woven Parameterization

<table>
    <tr align="middle" >
          <td width="70%" align="left">
          <i>"If your pair of jeans isn't even isotropic, why should its parameterization be?"</i> </br></br>

Previous UV mapping methods (ABF++, LSCM, ARAP, BFF, SCAF, AutoCuts, OptCuts, ...) are rotation invariant and minimize triangle distortion without any assumption on the material being flattened. These methods are unfit for anisotropic materials which stretch differently along different axes.  

This repository provides an implementation of the anisotropic parameterization described in TODO by TODO.
It's original intent is to accurately flatten woven textiles, for which the thread structure induces
anisotropy in their ability to stretch.</td>
          <td width="30%"><img style="float: right;" src="images/woven_viz.gif" margin="35px"></td>
    </tr>
</table>

## Testing

```
git clone https://github.com/corentinDumery/TODO
git submodule update --init --recursive
mkdir build && cd build
cmake ..
make -j woven_param       #build library
make -j woven_viz         #build test app
./woven_viz [path/to/input/mesh.obj]
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
git add submodule https://github.com/corentinDumery/TODO
git submodule update --init --recursive
``` 

Then, in your `CMakeLists.txt`, add: 
```
add_subdirectory([path_to_woven_param]/woven-param)
...
target_include_directories([target] PUBLIC [path_to_woven_param]/include)
target_link_libraries([target] woven_param)
```

## Citing

Thank you for reading! If this repository is useful to you, feel free to reach out and/or cite our paper:

TODO

![teaser](images/animals.png) 

## Acknowledgments

Some of the models used in this repository are adapted from the dataset associated with https://github.com/maria-korosteleva/Garment-Pattern-Generator