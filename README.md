# DenseCRF 3D

For nifti medical images.

Based on:
- mia and itkio libraries https://github.com/biomedia-mira/oak2
- Philipp Krähenbühl and Vladlen Koltun, "Efficient inference in fully connected crfs with gaussian edge potentials", in NIPS, 2011.
- https://github.com/deepmedic/dense3dCrf


## Docker


build

```
docker build -t densecrf .
```


run

```
export DATA_DIR=/mnt/data
docker run -v $DATA_DIR:/data -it --rm densecrf dense3DCrfInference --help
```

```
# /data/test.cfg
image = /data/vol.nii.gz
probmap = /data/probmap.nii.gz
output = /data/output
posW = 2
biW = 1
numberOfForegroundClasses = 1
biXYZStds = 5 5 5
posXYZStds = 3 3 3
biModStds = 30
maxIterations = 5
```

```
docker run -v $DATA_DIR:/data -it --rm densecrf dense3DCrfInference --unary --probs --config /data/test.cfg
```

```
docker run -v $DATA_DIR:/data -it --rm densecrf dense3DCrfInference \
--probs --unary --image /data/vol.nii.gz \
--probmap /data/probmap.nii.gz \
--output /data/output \
--numberOfForegroundClasses 1 \
--maxIterations 5 \
--biXYZStds "3 3 3" \
--posXYZStds "4 4 4" 
```