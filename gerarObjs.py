import os
import cv2
import mediapipe as mp

# Extrai a nuvem de pontos usando o Face Landmark Detection (MediaPipe)
def extract_point_cloud(frame):
    mp_face_mesh = mp.solutions.face_mesh
    face_mesh = mp_face_mesh.FaceMesh()

    rgb_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    points = []
    results = face_mesh.process(rgb_frame)
    if results.multi_face_landmarks:
        for face_landmarks in results.multi_face_landmarks:
            for i, landmark in enumerate(face_landmarks.landmark, start=0):
                x, y, z = landmark.x, landmark.y, landmark.z
                points.append((x, y, z))
    return points

# Extrai as faces que definem as conexões de cada vértice
def extract_faces(points, filename):
    faces = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('f '):
                face = [int(x) for x in line.strip().split()[1:]]
                faces.append(face)
    return faces

# Escreve os vértices e as faces no arquivo OBJ
def write_obj(vertices, faces, filename, image_path):
    with open(filename, 'w') as f:
        f.write(f"# {image_path}\n")
        for vertex in vertices:
            f.write(f"v {vertex[0]} {vertex[1]} {vertex[2]}\n")
        for face in faces:
            f.write(f"f {' '.join([str(idx+1) for idx in face])}\n")



def process_images(image_folder, pasta_destino):
    caminhos_imagens = [os.path.join(image_folder, arquivo) for arquivo in os.listdir(image_folder) if arquivo.endswith(('.png', '.jpg', '.jpeg'))]

    for idx, image_path in enumerate(caminhos_imagens, start=1):
        frame = cv2.imread(image_path)
        points = extract_point_cloud(frame)  
        filename = os.path.join(pasta_destino, f"face_{image_path.split('/')[1].split('.')[0]}.obj")  
        faces_path = './faces.obj' 
        faces = extract_faces(points, faces_path)  
        write_obj(points, faces, filename, image_path) 

process_images("imgs", "pasta_destino")
