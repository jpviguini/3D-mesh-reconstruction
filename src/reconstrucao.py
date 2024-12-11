import numpy as np
import matplotlib.pyplot as plt
from adjacencia import (
    ler_obj_vertices, 
    ler_obj_vertices_e_faces, 
    gerar_nova_malha, 
)
from scipy.spatial import KDTree

original_mesh_path = 'pasta_destino/face_rosto1.obj'
curve_path = 'robust_features/ridges_rosto1.obj'

### Algoritmo de reconstrução
def gen_default_weights(A):
    """
    Gera os pesos normalizados e os graus a partir da matriz de adjacência.
    """
    c = np.sum(A, axis=1)
    c = np.array(c)
    # Evitar divisão por zero
    c[c == 0] = 1
    w = A / c[:, np.newaxis]
    return w, c


def gen_delta_coords(v, w, c):
    """
    Calcula as coordenadas delta como L V, onde L é a matriz laplaciana.
    """
    # L = D - A, onde D é a matriz diagonal de graus
    # delta = L V = D V - A V
    delta = np.diag(c) @ v - w @ v * c[:, np.newaxis]
    return delta


def gen_laplacian(A):
    """
    Gera a matriz laplaciana a partir da matriz de adjacência.
    """
    return np.diag(np.sum(A, axis=1)) - A  # L = D - A


def recover_coords(L, delta, anchor_id, anchor_val):
    """
    Recupera as coordenadas dos vértices resolvendo o sistema linear L V = delta,
    com condições de âncora adicionadas.
    """
    n = L.shape[0]
    k, ndims = anchor_val.shape

    # Adicionar linhas de âncoras à matriz laplaciana
    add_L = np.zeros((k, n))
    add_L[np.arange(k), anchor_id] = 1
    L_aug = np.vstack((L, add_L))

    # Adicionar valores de âncoras às coordenadas delta
    delta_aug = np.vstack((delta, anchor_val))

    # Resolver o sistema para cada dimensão
    v = np.zeros((n, ndims))
    for i in range(ndims):
        v[:,i] = np.linalg.lstsq(L_aug, delta_aug[:,i], rcond=None)[0]
    return v


def get_error(original, reconstructed):
    """
    Calcula o erro de reconstrução como a norma Euclidiana entre as malhas original e reconstruída.
    """
    return np.linalg.norm(original - reconstructed)


def save_as_obj(vertices, faces, file_path):
    """
    Salva a malha reconstruída em um arquivo OBJ.
    """
    with open(file_path, 'w') as file:
        for vertex in vertices:
            file.write(f"v {vertex[0]} {vertex[1]} {vertex[2]}\n")
        for face in faces:
            face_str = ' '.join(str(idx + 1) for idx in face)  # +1 para formato OBJ
            file.write(f"f {face_str}\n")

## Encontra os vértices de fronteira da malha
def find_boundary_vertices(faces):
    """
    Identifica os vértices de fronteira na malha para isolarmos da reconstrução (a borda se mantém fixa).
    """

    edge_count = {}
    for face in faces:
        # Considera todas as arestas de cada face
        edges = [
            tuple(sorted((face[0], face[1]))),
            tuple(sorted((face[1], face[2]))),
            tuple(sorted((face[2], face[0]))),
        ]
        for edge in edges:
            if edge in edge_count:
                edge_count[edge] += 1
            else:
                edge_count[edge] = 1

    # Arestas que aparecem apenas uma vez são de fronteira
    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]

    # Vértices de fronteira são os extremos dessas arestas
    boundary_vertices = set()

    for edge in boundary_edges:
        boundary_vertices.update(edge)
    return list(boundary_vertices)


def main():
    # Porcentagem de pontos âncoras a serem selecionados (entre 0 e 100)
    porcentagem_ancoras = 100

    # Carregar os vértices e faces da malha original 
    vertices_malha, faces_malha = ler_obj_vertices_e_faces(original_mesh_path)                             

    # Carregar os vértices das curvas
    vertices_curva = ler_obj_vertices(curve_path)    
    curve_vertices = np.array(vertices_curva)

    # Gerar a nova malha incluindo os vértices da curva parabólica
    novos_vertices, novas_faces, matriz_adjacencia = gerar_nova_malha(vertices_curva, vertices_malha, faces_malha)

    # Convertendo a nova malha para um array numpy
    v = np.array(novos_vertices)
    n = len(v)

    # KD-Tree com os novos vértices (para encontrarmos quais vértices da curva estão na nova malha)
    kdtree = KDTree(v)

    # Definir a tolerância para a correspondência
    tolerance = 1e-1 

    # Inicializar listas para IDs e valores dos âncoras (nesse caso são os vértices da curva)
    anchor_id = []
    anchor_val = []

    # Iterar sobre os vértices da curva e adicionar como âncora apenas os encontrados na nova malha
    for vertex in curve_vertices:
        distance, idx = kdtree.query(vertex)
       
        if distance < tolerance:
            anchor_id.append(idx)
            anchor_val.append(v[idx])
        else:
            print(f"Vértice da curva {vertex} não encontrado na malha dentro da tolerância.")

    # Identificar os vértices de fronteira e adicionar como âncoras 
    boundary_vertex_indices = find_boundary_vertices(faces_malha)
    boundary_vertices = v[boundary_vertex_indices]
    anchor_id.extend(boundary_vertex_indices)
    anchor_val.extend(boundary_vertices)

    # Converter as listas dos âncoras para arrays numpy
    anchor_id = np.array(anchor_id)
    anchor_val = np.array(anchor_val)

    print("\n")
    print("Número total de âncoras encontradas na nova malha (curva + fronteira):", len(anchor_id))

    if len(anchor_id) == 0:
        print("Erro: Nenhum vértice da curva ou da fronteira foi encontrado na malha.")
        return

    # Selecionar uma porcentagem dos âncoras de forma linearmente espaçada
    num_ancoras_total = len(anchor_id)
    num_ancoras_selecionadas = max(1, int(np.ceil((porcentagem_ancoras / 100) * num_ancoras_total)))
    
    # Gera índices linearmente espaçados
    indices_selecionados = np.linspace(0, num_ancoras_total - 1, num=num_ancoras_selecionadas, dtype=int)
    
    # Seleciona os IDs e valores das âncoras
    anchor_id_selecionados = anchor_id[indices_selecionados]
    anchor_val_selecionados = anchor_val[indices_selecionados]

    print(f"Número de âncoras selecionadas ({porcentagem_ancoras}%):", len(anchor_id_selecionados))

    ### Aqui começa o algoritmo de reconstrução ###

    # Gerar os pesos e graus padrão
    w, c = gen_default_weights(matriz_adjacencia)

    # Inicializar delta como zeros (já que delta = L V e V será recalculado)
    delta = np.zeros((n, 3), dtype=np.float64)
    
    # Gerar a matriz laplaciana com base na matriz de adjacência que criamos com a nova malha
    L = gen_laplacian(matriz_adjacencia)

    # Reconstruir os vértices 
    v_rec = recover_coords(L, delta, anchor_id_selecionados, anchor_val_selecionados)

    # Salvar a malha reconstruída
    save_as_obj(v_rec, novas_faces, 'recovered.obj')

    print("Número de pontos da curva:", len(vertices_curva))
    print("\n")
    print("Tamanho da malha original:", len(vertices_malha))
    print("\n")
    print("Tamanho da nova malha:", len(v))
    print("Tamanho da malha reconstruida:", len(v_rec))
    print("\n")

    # Calcular e imprimir o erro da reconstrução
    error = get_error(v, v_rec)
    print("Erro da reconstrução:", error)

if __name__ == "__main__":
    main()


