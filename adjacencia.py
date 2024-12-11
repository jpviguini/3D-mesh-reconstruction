import numpy as np

def ler_obj_vertices_e_faces(caminho_arquivo):
    """
    Lê um arquivo .obj e extrai os vértices e as faces.
    
    Args:
        caminho_arquivo: caminho para o arquivo .obj.
        
    Returns:
        vertices: lista de vértices (cada vértice é uma tupla de coordenadas).
        faces: lista de faces (cada face é uma lista de índices de vértices).
    """
    vertices = []
    faces = []
    
    with open(caminho_arquivo, 'r') as arquivo:
        for linha in arquivo:
            partes = linha.strip().split()
            if not partes:
                continue
            if partes[0] == 'v':
                # Lê um vértice (x, y, z)
                vertices.append(tuple(map(float, partes[1:4])))
            elif partes[0] == 'f':
                # Lê uma face (índices de vértices)
                # Ajusta índices de 1-base para 0-base
                indices = [int(p.split('/')[0]) - 1 for p in partes[1:]]
                faces.append(indices)
    
    return vertices, faces

def ler_obj_vertices(caminho_arquivo):
    """
    Lê um arquivo .obj e extrai apenas os vértices (para a curva), ignorando vértices com coordenadas (0, 0, 0).
    
    Args:
        caminho_arquivo: caminho para o arquivo .obj.
        
    Returns:
        vertices: lista de vértices (cada vértice é uma tupla de coordenadas).
    """
    vertices = []
    
    with open(caminho_arquivo, 'r') as arquivo:
        for linha in arquivo:
            partes = linha.strip().split()
            if not partes:
                continue
            if partes[0] == 'v':
                # Lê um vértice (x, y, z)
                coords = tuple(map(float, partes[1:4]))
                # Ignora vértices com coordenadas (0, 0, 0)
                if coords != (0.0, 0.0, 0.0):
                    vertices.append(coords)
    
    return vertices

def save_as_obj(vertices, faces, file_path):
    with open(file_path, 'w') as file:
        for vertex in vertices:
            file.write(f"v {vertex[0]} {vertex[1]} {vertex[2]}\n")
        for face in faces:
            face_str = ' '.join(str(idx + 1) for idx in face)  # soma 1 porque o .obj começa em 1, não em 0.
            file.write(f"f {face_str}\n")


def obter_arestas(face, vertices):
    """
    Retorna pares de vértices (arestas) de uma face.

    Args:
        face: Lista de índices dos vértices da face.
        vertices: Lista de coordenadas dos vértices.

    Returns:
        List[Tuple[np.array, np.array]]: Lista de arestas, onde cada aresta é um par de vértices.
    """
    num_vertices = len(face)
    arestas = []
    for i in range(num_vertices):
        v1 = vertices[face[i]]
        v2 = vertices[face[(i + 1) % num_vertices]]
        arestas.append((v1, v2))
    return arestas


def ponto_intersecta_aresta(ponto, aresta, tolerancia=1e-6):
    """
    Verifica se o ponto está em uma aresta específica (v1, v2).
    
    Args:
        ponto: coordenadas do ponto de interseção (np.array ou lista).
        aresta: tupla com dois vértices (v1, v2), onde v1 e v2 são np.array ou listas.
        tolerancia: tolerância para a verificação de colinearidade e de pertencimento ao segmento.
    
    Returns:
        bool: True se o ponto estiver na aresta, False caso contrário.
    """
    v1, v2 = map(np.array, aresta)
    ponto = np.array(ponto)
    
    # Vetores da aresta e do ponto em relação ao v1
    vetor_aresta = v2 - v1
    vetor_ponto = ponto - v1

    # Calcula a norma do vetor aresta para normalização
    norma_aresta = np.linalg.norm(vetor_aresta)
    if norma_aresta < tolerancia:
        return False  # Aresta de comprimento zero

    # Projeta o vetor ponto no vetor aresta para encontrar a posição escalar
    t = np.dot(vetor_ponto, vetor_aresta) / (norma_aresta ** 2)

    # Verifica se t está no intervalo [0, 1]
    if t < -tolerancia or t > 1 + tolerancia:
        return False

    # Calcula a distância perpendicular do ponto à aresta
    proj = v1 + t * vetor_aresta
    distancia = np.linalg.norm(proj - ponto)

    return distancia < tolerancia

def construir_matriz_adjacencia(faces, num_vertices):
    """
    Constrói uma matriz de adjacência como uma matriz 2D.

    Args:
        faces: lista de faces, onde cada face é uma lista de índices de vértices.
        num_vertices: número total de vértices na malha.

    Returns:
        np.array: matriz de adjacência de tamanho (num_vertices x num_vertices).
    """
    adjacencia = np.zeros((num_vertices, num_vertices), dtype=int)

    for face in faces:
        num_vertices_face = len(face)
        for i in range(num_vertices_face):
            v1 = face[i]
            v2 = face[(i + 1) % num_vertices_face]

            # Marca a adjacência entre v1 e v2
            adjacencia[v1, v2] = 1
            adjacencia[v2, v1] = 1  

    return adjacencia



## Gera a nova malha com os vértices da curva + vértices do mediapipe
def gerar_nova_malha(curve_intersections, mesh_vertices, mesh_faces):
    new_vertices = list(mesh_vertices)  # Os novos vértices iniciam com o os vértices do mediapipe
    new_faces = []
    intersection_point_indices = {}  #  Dicionário para mapear cada vértice da curva a um índice

    # Converte os vértices da curva para array numpy
    curve_points = [np.array(pt) for pt in curve_intersections]


    # Para cada face, adiciona os vértices da curva na nova malha (se intersectarem alguma aresta da face)
    for face in mesh_faces:
        arestas = obter_arestas(face, mesh_vertices)
        pontos_intersecao = []

        for idx_aresta, aresta in enumerate(arestas): # para cada aresta
            for pt in curve_points:  # para cada ponto da curva

                if ponto_intersecta_aresta(pt, aresta):  # verifica se o ponto está intersectando a aresta atual
                    pt_tuple = tuple(pt)

                    # Verifica se o ponto de intersecao ja foi registrado
                    if pt_tuple not in intersection_point_indices:                       
                        intersection_point_indices[pt_tuple] = len(new_vertices)
                        new_vertices.append(pt_tuple) # adiciona o ponto da curva aos vértices da nova malha

                    pontos_intersecao.append((idx_aresta, intersection_point_indices[pt_tuple]))
                    break  # Para ao encontrar uma interseção na aresta


        # Dicionário onde cada chave é um edge_idx (índice da aresta) e o valor é idx_p (índice do ponto de interseção).
        pontos_intersecao = list({edge_idx: idx_p for edge_idx, idx_p in pontos_intersecao}.items())

        # Divide a face com base no número de pontos de interseção nas arestas
        novos_faces = subdividir_face(face, pontos_intersecao, new_vertices)
        new_faces.extend(novos_faces)

    num_vertices = len(new_vertices)

    # Gera a matriz de adjacência com a nova malha criada
    nova_adjacencia = construir_matriz_adjacencia(new_faces, num_vertices)  

    return new_vertices, new_faces, nova_adjacencia


def subdividir_triangulo(triangle, pontos_intersecao, vertices):
    """
    Subdivide o triângulo com base nos pontos de interseção da curva

    Args:
        triangle: Lista de índices dos vértices do triângulo da malha original
        pontos_intersecao: Lista de tuplas (edge_idx, idx_ponto_intersecao).
        vertices: Lista atualizada de vértices

    Returns:
        novas_faces: Lista das novas faces resultantes da divisão do triângulo
    """
    novas_faces = []

    num_intersections = len(pontos_intersecao)
    v0, v1, v2 = triangle

    # Define os índices das arestas em um dicionário
    edge_vertices = {
        0: (v0, v1),
        1: (v1, v2),
        2: (v2, v0)
    }

    ## A curva não intersecta o triângulo ---> mantém o triângulo original                        
    if num_intersections == 0:
        novas_faces.append([v0, v1, v2])

    # Uma interseção --> divide o triângulo em dois triângulos
    elif num_intersections == 1:
        
        (edge_idx, idx_p) = pontos_intersecao[0]

        # Vértices da aresta insertectada
        va, vb = edge_vertices[edge_idx]

        # Vértice oposto da aresta intersectada
        vc = list(set([v0, v1, v2]) - set([va, vb]))[0]

        # Cria dois novos triângulos
        novas_faces.append([va, idx_p, vc])
        novas_faces.append([idx_p, vb, vc])
        #print("1 intersecao")

    ## Duas interseções (onde mais acontece!!) --> cria um triângulo e um quadrilátero
    elif num_intersections == 2:
        
        (edge_idx1, idx_p1), (edge_idx2, idx_p2) = pontos_intersecao
        
        # Identifica a aresta que não é intersectada
        non_intersected_edge_idx = list(set([0, 1, 2]) - set([edge_idx1, edge_idx2]))[0]

        # O vértice oposto a aresta não intersectada é o que não faz parte dessa aresta
        va, vb = edge_vertices[non_intersected_edge_idx]
        opposite_vertex = list(set([v0, v1, v2]) - set([va, vb]))[0]

        # Forma um triângulo com os pontos de interseção e o vértice oposto
        novas_faces.append([idx_p1, idx_p2, opposite_vertex])

        # Forma um quadrilátero
        # Usa os vertices da aresta intersectada que são diferentes do vértice oposto
        va1, vb1 = edge_vertices[edge_idx1]
        va2, vb2 = edge_vertices[edge_idx2]
        vertex_a = va1 if va1 != opposite_vertex else vb1
        vertex_b = va2 if va2 != opposite_vertex else vb2

        # Forma o quadrilátero
        quad_face = [idx_p1, vertex_a, vertex_b, idx_p2]
        novas_faces.append(quad_face)
    
    ## 3 Interseções ---> cria quatro triângulos
    elif num_intersections == 3:
        
        #print("3 interseções no triângulo")

        # Obter os pontos de interseção
        idx_p1 = pontos_intersecao[0][1]
        idx_p2 = pontos_intersecao[1][1]
        idx_p3 = pontos_intersecao[2][1]

        # Criar o triângulo interno
        novas_faces.append([idx_p1, idx_p2, idx_p3])

        # Criar três triângulos externos
        novas_faces.append([v0, idx_p1, idx_p3])
        novas_faces.append([v1, idx_p2, idx_p1])
        novas_faces.append([v2, idx_p3, idx_p2])

    ## 4 ou mais intersecoes --> muito improvável de acontecer, vamos manter o mesmo triangulo original
    else:
        print("4 ou mais intersecoes")
        novas_faces.append([v0, v1, v2])

    return novas_faces



## Função que divide a face dependendo do número de vértices da face e pontos de interseção.
def subdividir_face(face, pontos_intersecao, vertices):
    novas_faces = []

    num_vertices_face = len(face)
    num_intersections = len(pontos_intersecao)

    # Define os vértices das arestas
    edge_vertices = {}
    for i in range(num_vertices_face):
        edge_vertices[i] = (face[i], face[(i + 1) % num_vertices_face])

    if num_intersections == 0: ## Nao tem interseçao, mantem a mesma face
        novas_faces.append(face)
    
    ## Divide o triângulo
    elif num_vertices_face == 3:
        novas_faces = subdividir_triangulo(face, pontos_intersecao, vertices)

    ## Divide o quadrilátero (para uma malha mais refinada)
    elif num_vertices_face == 4:
        pass
    else:
        print(f"Faces com {num_vertices_face} vértices não são suportadas.")
        novas_faces.append(face)

    return novas_faces


def main():
    # Carrega a malha do mediapipe
    vertices_malha, faces_malha = ler_obj_vertices_e_faces('pasta_destino/face_rosto1.obj')

    # Carrega os vértices da curva
    vertices_curva = ler_obj_vertices('robust_features/parabolicas_rosto1.obj')

    # Gera a nova malha, que agora inclui os pontos da curva
    novos_vertices, novas_faces, matriz_adjacencia = gerar_nova_malha(vertices_curva, vertices_malha, faces_malha)
    #print(len(novos_vertices))

    # Salva a nova malha em um arquivo .obj
    save_as_obj(novos_vertices, novas_faces, 'nova_malha.obj')

if __name__ == "__main__":
    main()