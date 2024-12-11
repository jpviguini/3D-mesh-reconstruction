import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def read_obj(file_path):
    vertices = []
    faces = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('v '):  # Vertices
                vertex = list(map(float, line.strip().split()[1:]))
                vertices.append(vertex)
            elif line.startswith('f '):  # Faces
                face = [int(part.split('/')[0]) - 1 for part in line.strip().split()[1:] if part]  
                faces.append(face)

    return np.array(vertices), faces

def plot_mesh(ax, vertices, faces, title='Mesh'):
 
    poly3d = [[vertices[vertex] for vertex in face] for face in faces]
  
    mesh = Poly3DCollection(poly3d, alpha=0.5, edgecolor='k', facecolor='cyan')
    ax.add_collection3d(mesh)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)

    x_limits = [vertices[:, 0].min() - 0.1, vertices[:, 0].max() + 0.1]
    y_limits = [vertices[:, 1].min() - 0.1, vertices[:, 1].max() + 0.1]
    z_limits = [vertices[:, 2].min() - 0.1, vertices[:, 2].max() + 0.1]

    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_zlim(z_limits)

def main():
    
    file_path1 = 'pasta_destino/face_rosto1.obj' 
    file_path2 = 'recovered.obj'  


    vertices1, faces1 = read_obj(file_path1)
    vertices2, faces2 = read_obj(file_path2)


    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')


    plot_mesh(ax1, vertices1, faces1, title='Malha original (Mediapipe)')
    plot_mesh(ax2, vertices2, faces2, title='Reconstrução')

    plt.show()

if __name__ == "__main__":
    main()
