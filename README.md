# blenderPolygons

Welcome to the blenderPolygons wiki!

## What does it do ?

Detect polygons from unordered splines like in architectural 2d drawings within a single click. Fix precision issues filling small gaps between lines encountred when importing. Once detected, you are able to pick polygons with a select tool and extract curves from your selection.

## Why ?

Cleaning dxf files or redrawing walls in worst cases is often a painfull task. This script does automagically fix most issues encountred in such cases (not well formed lines), and allow you to simply pick needed polygons and get a clean curve of polygons boundarys.

## How does it do ?

Under the hood, i do write from scratch the cleaning code and intersection detection to handle precision and optimization issuees occuring on medium data sets using a spacial tree (rtree). Polygonization is only a wrapper arround shapely as it is realy powerfull and quick.

## Dependency

In order to work this script depends on
- shapely

available as python wheels.

## Setup under windows:

The most easiest way to install shapely Python Binding is to use the packages available here rtree shapely. Choose the one that match the version of Python bundle with Blender eg: cp35-cp35m versions for blender 2.78.

You need pip to install *.whl package files. The Python installation bundle with Blender do not include pip but include distutils. So it's possible to install pip with a standard console:

blender_install_folder\2.7x\python\bin\python.exe -m ensurepip

And then install the wheel:

blender_install_folder\2.7x\python\bin\python.exe -m pip install xxx-cp35-cm35m-win_amd64.whl

To test the install open Blender Python console and type:

from shapely import geometry

