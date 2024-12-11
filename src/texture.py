import numpy as np
from PIL import Image
import trimesh
from mayavi import mlab
from tvtk.api import tvtk

recovered_mesh = 'recovered.obj'
original_image = 'imgs/rosto1.jpg'

# Carrega a malha a partir do arquivo .obj
mesh = trimesh.load(recovered_mesh)

# Carrega a imagem original
imagem = Image.open(original_image)
largura_imagem, altura_imagem = imagem.size
dados_imagem = np.array(imagem)

# Obtém as coordenadas dos vértices
vertices = mesh.vertices.copy()

# Extrai as coordenadas X e Y normalizadas dos vértices
x_coords = vertices[:, 0]
y_coords = vertices[:, 1]

# Verifica se as coordenadas estão no intervalo [0, 1], caso contrário, normaliza
if x_coords.min() < 0 or x_coords.max() > 1 or y_coords.min() < 0 or y_coords.max() > 1:
    x_coords = (x_coords - x_coords.min()) / (x_coords.max() - x_coords.min())
    y_coords = (y_coords - y_coords.min()) / (y_coords.max() - y_coords.min())

# Multiplica pelas dimensões da imagem para obter coordenadas em pixels
u = x_coords * (largura_imagem - 1)
v = y_coords * (altura_imagem - 1)

# Converte para inteiros e garante que estejam dentro dos limites da imagem
pixel_x = np.clip(u.round().astype(int), 0, largura_imagem - 1)
pixel_y = np.clip(v.round().astype(int), 0, altura_imagem - 1)

# Obtém as cores dos pixels correspondentes aos vértices
vertex_colors = dados_imagem[pixel_y, pixel_x, :3]

# --- PLOTA A MALHA GERADA COM PHONG SHADING USANDO MAYAVI ---

# Extrai os vértices e faces da malha
faces = mesh.faces.copy()
points = vertices

# Cria um objeto PolyData do TVTK
poly_data = tvtk.PolyData(points=points, polys=faces)

# Atribui as cores aos vértices
vertex_colors_uint8 = vertex_colors.astype(np.uint8)
poly_data.point_data.scalars = vertex_colors_uint8
poly_data.point_data.scalars.name = 'colors'

# Configura a visualização com o Mayavi
mlab.figure(size=(800, 600))

# Cria a superfície com Phong shading
mesh_plot = mlab.pipeline.surface(poly_data)
mesh_plot.actor.property.interpolation = 'phong'  # Aplica o Phong shading

# Ajusta as propriedades de iluminação (opcional)
mesh_plot.actor.property.specular = 0.1
mesh_plot.actor.property.specular_power = 10

# Exibe a malha
mlab.show()
