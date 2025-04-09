# Interface-Maker

[![DOI](https://zenodo.org/badge/DOI/xxx.xxx/xxx.xxx.svg)](https://doi.org/xxx.xxx/xxx.xxx)

## Overview

An application for generating customizable slabs and interfaces for first-principles simulations.

This code is based on the A. Zur et al. paper: "Lattice match: An application to heteroepitaxy, Journal of applied physics 55(2) (1984) 378-386".

## Get Started

* 🌐 [Try it Online](https://interface-maker.streamlit.app/)

* ⭐️ Please STAR this repository if you find it helpful :)

* For any questions or suggestions, please contact us (gliu4@wpi.edu; yzhong@wpi.edu). 

## Workflow

![Workflow](res/image.png)

## Usage Description

1. **Upload Your Structures**: In the `Upload Your Structures` section, you can upload the conventional cell of the lower and upper systems in VASP format (.vasp/.poscar). After successfully uploading the files, you will see the structure of the lower and upper systems displayed below the section. 

2. **Define Miller Indices**: In the `Define Miller Indices` section, you can define the maximum Miller indices of h, k, l for lower and upper slabs. If you are interested in specific Miller indices, you can assign them by checking the `Assign Specific Miller Indices for Lower and Upper Systems` checkbox. This will allow you to input the Miller indices directly.

3. **Set Lattice Match Parameters**: In the `Set Lattice Match Parameters` section, you can set the tolerance for the misfit of lattice vectors (in %) and angle (in degree). The default values are 5% for the misfit of lattice vectors and 5 degrees for the angle. You can adjust these values according to your needs.

4. **Customize the Interface Geometry**: In the `Customize the Interface Geometry` section, you can customize the interface geometry including the minimum thickness of the slab, slab vacuum, interface gap, area range for the matched interfaces, and the shape filter option. Note that checking the `Shape Filter` option will only keep the square-like interfaces, which will help to reduce the number of interfaces generated.

5. **Generate Interfaces**: After setting all the parameters, click the `Generate Interfaces` button. The program will generate all possible interfaces based on the parameters you set. The interface generation details will be displayed in the `Interface Generation Details` section.

6. **Download the Generated Interfaces**: After the interface generation is complete, you can download the generated interfaces in a zip file. The zip file will contain the generated interface structures in VASP format (.vasp/.poscar) and profiles in .txt and .csv formats. The profiles include the atoms, Miller indices, area, lattice matching misfit, and transformation matrix.

## Cite Us
If you use this code in your research, please cite our paper:

```
@article{
  title={A Computational Framework for Interface Design: From First Principles to Machine Learning, with a Case Study on LaCoO3/La2NiO4},
  author={G. Liu and S. Yang and Y. Zhong},
  journal={},
  volume={},
  number={},
  pages={},
  year={2025}
}
```