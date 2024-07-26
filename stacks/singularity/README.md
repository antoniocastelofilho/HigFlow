## First

You must have singularity already installed in your computer. 

See singularity developer page for instruction on how to install singularity: https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps

## Second

Build the singularity image

Bind where the bibliotecas directory (inside Higflow main directory) and scripts are  with the container

```
singularity build --bind PATH_TO_BIBLIOTECAS:/bibliotecas,$PWD/scripts:/scripts higflow-image-ubuntu22.sif higflow_image.def
```

Ex:

```
singularity build --bind $USER/HigFlow/bibliotecas:/bibliotecas,$PWD/scripts:/scripts higflow-image-ubuntu22.sif higflow_image.def
```