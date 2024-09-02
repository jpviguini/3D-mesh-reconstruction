# 3D Reconstruction of Human Faces: a Laplacian approach 
Research project supported by FAPESP from January 2024 to December 2024.

**Developer:** João Pedro Viguini Tolentino Taufner Correa

**Supervisor:** João do Espírito Santo Batista Neto

## About this project [in development...]
The reconstruction of meshes from a reduced set of points constitutes the core of this work. The objective is to reconstruct a human face (including the nose, mouth, eyes, and other elements) from a simplified version, such as a caricature represented by simplified curves of the face. This reduced information is referred to as robust features. These features will be extracted from the three-dimensional mesh of a human face and will be used to reconstruct the 3D surface through the Laplacian operator.

## Getting started

- I'm currently implementing a **_python _notebook__** containing each step of the algorithms in this repository: https://colab.research.google.com/drive/1lm-FyJcLh4ZcomAdfGPfnJpKkfKnvTTL#scrollTo=Q1EQNd8o0u2D

- Make sure to install all the necessary dependencies listed in the **requirements.txt** file.

**Disclaimer:** This GitHub page is still in **development** phase...

## Summary
0. Introduction
1. Point cloud extraction
2. Robust features
3. Surface reconstruction
4. Results

## 0. Introduction
The main goal of this repository is to provide hands-on experience to apply the surface reconstruction technique explored in this research project. **To get more details about the project, please check the full article attached on this repository (fully in portuguese)**.

## 1. Point cloud extraction
The acquisition of the three-dimensional point cloud of the human face is essential for surface reconstruction. For this task, the MediaPipe tool was used, which employs machine learning models to detect landmarks and facial expressions from images or image sequences (videos).

These models are capable of operating on both single images and continuous image streams, providing real-time results. The input data for point cloud extraction was an image (PNG, JPG, or JPEG) of a human face in a frontal position.


<img src="https://github.com/user-attachments/assets/90dd2d55-44f9-44f7-936c-c141b72aca7b" width="300">

- MediaPipe is capable of detecting 468 points.
- This task will generate a .obj file containing all the vertices and faces.

Each face determines the connection of three vertices, calculated by applying the Delaunay Triangulation algorithm.

Several face files can define different topologies for the face. In this project, the same face file (faces.obj) was used to maintain consistency across all reconstructions. It is worth noting that a different topology would result in a different outcome, even for the same set of points.

## 2. Robust features 
Two main features can be extracted from a face: parabolic curves and _ridges_. These are known as robust features. These features do not change when the surface is deformed.

<img src="https://github.com/user-attachments/assets/39788517-0bf7-47b1-84cb-1d9d7e2334c1" width="300">

<img src="https://github.com/user-attachments/assets/c64859aa-3801-4ab4-9b74-f7fbc3a87bca" width="300">

The MATLAB code that extracts these features was originally developed by Danilo Marques during his doctorate:
- **main24atu2.m:** Extracts parabolic curves.
- **main.m**: Extracts _ridges_.

To execute these scripts, you should first provide the .obj file generated in the previous session.

## 3. Surface reconstruction
After obtaining the anchor points (parabolic curves or _ridges_), the Laplace operator is applied to reconstruct the surface from these points. You can also choose how many (%) anchor points to use. The anchor points are selected using the `linspace` function from the NumPy library, which allows us to select linearly spaced points, providing a well-distributed spatial arrangement along the human face.

### Why Laplace?
The approach is based on the Laplace operator and differential representations among the vertices of a given neighborhood in a point cloud mesh. In contrast to the traditional representation by global Cartesian coordinates, the differential representation of a surface reveals information about its local shape, as well as the size and orientation of local details. In addition to providing information that results in a reconstruction with better detail preservation, it is a linear system, which makes it computationally efficient.


## 4. Results
In development...


## Acknowledgements
I would like to thank professors João do E.S. Batista Neto, Antônio Castelo and Farid Tari for their guidance and support throughout the development of this project. I also want to thank Matheus Paiva for his support.

## Bibliography
- ANGAROLA, M. P. Curvature Estimation Using Machine Learning Algorithms. 2024. Last access: 19-08-2024. Available on: ⟨https://github.com/MatheusPaivaa/CurvatureML ⟩.

- BRUCE, J. W.; GIBLIN, P. J.; TARI, F. Ridges, crets an sub-parabolic lines of evolving
surfaces. International Journal of Computer Vision, 1996.

- FAKHOURY, A. Image curves reconstruction by means of robust features. 2021. Last access: 18-08-2024. Available on: ⟨https://github.com/andrefakhoury/image-curve-reconstruction ⟩.

- IZUMIYA, S. et al. Differential Geometry From A Singularity Theory Viewpoint.
Singapura: World Scientific, 2015. 139 p.

- LUGARESI, C. et al. Mediapipe: A framework for perceiving and processing reality. In:
Third Workshop on Computer Vision for AR/VR at IEEE Computer Vision and Pattern
Recognition (CVPR) 2019. [s.n.], 2019. Available on ⟨https://mixedreality.cs.cornell.edu/s/NewTitle May1 MediaPipe CVPR CV4ARVR Workshop 2019.pdf ⟩.

- SORKINE, O. Differential representations for mesh processing. Computer Graphics
Forum, European Association for Computer Graphics, v. 25, n. 4, p. 789–807, 2006.
YANG, S. et al. Wider face: A face detection benchmark. In: IEEE Conference on
Computer Vision and Pattern Recognition (CVPR). [S.l.: s.n.], 2016
