# 3D Reconstruction of Human Faces: a Laplacian Approach 
Research project supported by FAPESP from January 2024 to December 2024.

**Developer:** João Pedro Viguini Tolentino Taufner Correa

**Supervisor:** João do Espírito Santo Batista Neto

## About this project [in development]
The reconstruction of meshes from a reduced set of points constitutes the core of this work. The objective is to reconstruct a human face (including the nose, mouth, eyes, and other elements) from a simplified version, such as a caricature represented by simplified curves of the face. This reduced information is referred to as robust features. These features will be extracted from the three-dimensional mesh of a human face and will be used to reconstruct the 3D surface through the Laplacian operator.

## Getting started

- There is a **_python _notebook__** containing each step of the algorithms in this repository: https://colab.research.google.com/drive/1lm-FyJcLh4ZcomAdfGPfnJpKkfKnvTTL#scrollTo=Q1EQNd8o0u2D

- Make sure to install all the necessary dependencies listed in the **requirements.txt** file.

## Summary
1. Point cloud extraction
2. Robust features
3. Surface reconstruction
4. Results

## 1. Point cloud extraction
The acquisition of the three-dimensional point cloud of the human face is essential for surface reconstruction. For this task, the MediaPipe tool was used, which employs machine learning models to detect landmarks and facial expressions from images or image sequences (videos).

These models are capable of operating on both single images and continuous image streams, providing real-time results. The input data for point cloud extraction was an image (PNG, JPG, or JPEG) of a human face in a frontal position.

[point_cloud_image]

- MediaPipe is capable of detecting 468 points.
- This task will generate a .obj file containing all the vertices and faces.

Each face determines the connection of three vertices. They are calculated by applying the Delaunay Triangulation algorithm.

[delaunay triangulation image]

There are several face files that can define different topologies for the face. In this project, the same face file (**faces.obj**) was used to maintain the consistency of all reconstructions. It is worth noting that a different topology would result in a different outcome, even for the same set of points.

## 2. Robust features 
There are two main features that you can extract from a face: parabolic curves and _ridges_. These are known as **robust features**.
These are features that doesn't change when the surface is deformed.

[parabolic image]

[ridges image]




### Parabolic curves

### Ridges

## 3. Surface reconstruction

### Choosing the anchor points


## 4. Results

## 5. 3D characters


## Acknowledgements


## Bibliography

