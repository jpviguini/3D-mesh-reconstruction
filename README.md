# 3D Reconstruction of Human Faces: a Laplacian approach 🧮
<p align="center">
  <a href="#introduction">Introduction</a> • 
  <a href="#point_cloud">Point cloud extraction</a> • 
  <a href="#robust_features">Robust features</a> • 
  <a href="#mesh_update">Updating the mesh</a> •
  <a href="#reconstruction">Surface reconstruction</a> •
  <a href="#results">Results</a>
</p>

---

Research project supported by FAPESP from January 2024 to December 2024.

**Developer:** João Pedro Viguini Tolentino Taufner Correa

**Supervisor:** João do Espírito Santo Batista Neto

## About this project 
The reconstruction of meshes from a reduced set of points constitutes the core of this work. The objective is to reconstruct a human face (including the nose, mouth, eyes, and other elements) from a simplified version, such as a caricature represented by simplified curves of the face. This reduced information is referred to as robust features. These features will be extracted from the three-dimensional mesh of a human face and will be used to reconstruct the 3D surface through the Laplacian operator.

**Note (full article)**: To get more details about the project, please check the full article attached on this repository (fully in portuguese) - `final_report.pdf`.

## Getting started
Follow the steps below to set up the environment and get started with the project:

### Prerequisites
Ensure you have the following software installed:
- Python 3.8 or higher
- MATLAB or Octave (required for robust feature extraction)
- Git

### Installation

1. **Clone the repository**
   ```
   git clone https://github.com/jpviguini/3D-mesh-reconstruction.git
   cd 3D-mesh-reconstruction
   ```
2. **Install the required Python dependencies**
   ```
   pip install -r requirements.txt
   ```

### Usage

1. **Point Cloud Extraction**  
   Execute `gerarObjs.py` section to generate the `.obj` file (original mesh) from a static image.
   
   **Note:** There's an example (from `rosto1.png`) inside `pasta_destino` folder - in case you don't want to run this step.

3. **Feature Extraction**  
   Run the MATLAB scripts:  
   - `main24atu2.m` for parabolic curves.  
   - `main.m` for ridges.  
   These scripts should be executed on the `.obj` files to generate the robust features.

  **Note:** If you don't want to run these MATLAB scripts you can simply use the `.obj` files that already contains robust features in the folder `robust_features`. These curves are from the `rosto1.png` file, so if you want to add texture to it, you should use this image to capture the right colors.

3. **Surface Reconstruction**  
   Execute the reconstruction algorithm (`reconstrucao.py`) to generate the final 3D mesh. You will get a `recovered.obj` file.

4. **Visualization**
   Execute `plotar.py` with the `.obj` files (original and recovered meshes) to see the original mesh and the reconstructed mesh side-by-side.

5. **Texture (Bonus)**
   If you want to add texture to the mesh, simply execute `texture.py` with the right paths to the .obj files and the original image.


## <div id="introduction"></div>0. Introduction
The main goal of this repository is to provide hands-on experience to apply the surface reconstruction technique explored in this research project. Full article: `final_report.pdf`

## <div id="point_cloud"></div>1. Point cloud extraction
The acquisition of the three-dimensional point cloud of the human face is essential for surface reconstruction. For this task, the MediaPipe tool was used, which employs machine learning models to detect landmarks and facial expressions from images or image sequences (videos).

These models are capable of operating on both single images and continuous image streams, providing real-time results. The input data for point cloud extraction was an image (PNG, JPG, or JPEG) of a human face in a frontal position.



<img src="https://github.com/user-attachments/assets/f92f1ea7-8e33-4f06-adb1-211df30be488" width="300">

- MediaPipe is capable of detecting 468 points.
- This task will generate a .obj file containing all the vertices and faces.


In this original mesh, each face determines the connection of three vertices. The set of vertices that make up a mesh can be connected in different ways, resulting in distinct topologies for the same object, such as a face. _Mediapipe_, for example, always uses the same connections between the vertices, ensuring a consistent topology. However, it is important to note that different topologies produce different results, even if the set of vertices remains unchanged, as the way the vertices are connected alters the properties of the resulting mesh.

## <div id="robust_features"></div>2. Robust features
Two main features can be extracted from a face: parabolic curves and _ridges_ (blue and red)). These are known as **robust features**. These features do not change when the surface is deformed.

<img src="https://github.com/user-attachments/assets/eb37081d-0c9e-414a-af7c-b50fe4854c1b" width="300">


<img src="https://github.com/user-attachments/assets/4599240d-0e99-4bb6-8aa5-714147792965" width="300">


The MATLAB code that extracts these features was originally developed by Danilo Marques during his doctorate:
- **main24atu2.m:** Extracts parabolic curves.
- **main.m**: Extracts _ridges_.

To execute these scripts, you should first provide the .obj file generated by _Mediapipe_ in the previous session.

## <div id="mesh_update"></div>3. Updating the mesh
To apply the reconstruction algorithm, it is essential to update the original mesh to include the vertices resulting from the intersections with the curves. The new mesh generates an updated adjacency matrix, which contains crucial information for the anchor points in the surface reconstruction. Since the curve points necessarily intersect the edges of the original mesh, composed of triangular faces, the following intersection scenarios occur:

- **0 intersections:** No points from the curve intersect the face of the original mesh. In this case, the face remains unchanged.

- **1 intersection:** Only one point from the curve intersects one of the face's edges. The original triangle is divided into two new triangles, with the intersection point connected to the vertex opposite the intersected edge.

- **2 intersections:** Two points intersect the face. In this case, the two points are connected, creating one triangle and one quadrilateral. This is the most common scenario.

- **3 intersections:** Three points intersect the face, each on a different edge. The face is divided into four new triangles, connecting the intersection points to the original vertices of the face.

- **4 or more intersections:** Although rare, if there are four or more intersections in a single face, the original triangle is maintained for simplification. However, subsequent refinement is recommended for proper handling.

These scenarios ensure appropriate subdivision of triangular faces. For a more refined mesh, it is necessary to extend the handling to the division of other polygonal shapes, such as quadrilaterals and pentagons.

With the updated mesh, it is possible to extract the new adjacency matrix, which is fundamental for capturing the relationships between the curve vertices and the mesh vertices. This matrix is essential for surface reconstruction, as it provides the necessary information for the reconstruction algorithm, ensuring better integration of the curve vertices with the new mesh.

The image below shows the updated mesh after inserting the vertices from the parabolic curve into the original mesh from _Mediapipe_ and texturizing it with the colors from the original image.

<img src="https://github.com/user-attachments/assets/f4fe3867-995d-4520-b8fe-707cbde1e3e3" width="400">




## <div id="reconstruction"></div>4. Surface reconstruction
From the selected anchor points on the updated mesh, we apply the reconstruction algorithm to obtain the reconstructed mesh. This approach has previously been used for curves in a two-dimensional context (FAKHOURY, 2021). Similarly, we will use it for curves in three dimensions.

Assuming that the OBJ files containing the vertices of the parabolic curves or ridges have already been computed, the algorithm can be summarized in the following steps:

- Select the desired number of anchor points;

- Compute the adjacency matrix;

- Generate the weights for each vertex (default for all vertices);

- Calculate the differential coordinates;

- Generate the Laplacian matrix;

- Retrieve the Cartesian coordinates from the differential coordinates and obtain the reconstructed mesh.


### Why Laplace?
The approach is based on the Laplace operator and differential representations among the vertices of a given neighborhood in a point cloud mesh. In contrast to the traditional representation by global Cartesian coordinates, the differential representation of a surface reveals information about its local shape, as well as the size and orientation of local details. In addition to providing information that results in a reconstruction with better detail preservation, it is a linear system, which makes it computationally efficient.


## <div id="results"></div>5. Results

### Using parabolic curves

![reconstruct_parabolic_en](https://github.com/user-attachments/assets/48aad452-38d9-4ac8-92b4-0346502bc93d)

### Using _ridges_

![reconstruction_ridges_en](https://github.com/user-attachments/assets/a124fefe-0969-4c16-a0cc-4fb7370d8b59)



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
