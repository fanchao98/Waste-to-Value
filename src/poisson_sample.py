import sys
import open3d as o3d
import numpy as np
from pathlib import Path

target_mesh_name = sys.argv[1]

base_dir = Path(".") / "output" / target_mesh_name

mesh_path = base_dir / "D_S.off"
mesh = o3d.io.read_triangle_mesh(str(mesh_path))
mesh.compute_vertex_normals()

triangle_areas = o3d.geometry.TriangleMesh.get_surface_area(mesh)
print("D_S Mesh:", triangle_areas)

pcd = mesh.sample_points_poisson_disk(
    number_of_points=int(triangle_areas / 4),
    init_factor=10,
    use_triangle_normal=True
)

output_path = base_dir / "D_S.xyzn"
o3d.io.write_point_cloud(str(output_path), pcd)

mesh_path = base_dir / "D_A.off"
mesh = o3d.io.read_triangle_mesh(str(mesh_path))
mesh.compute_vertex_normals()

triangle_areas = o3d.geometry.TriangleMesh.get_surface_area(mesh)
print("D_A Mesh:", triangle_areas)

pcd = mesh.sample_points_poisson_disk(
    number_of_points=int(triangle_areas / 4),
    init_factor=10,
    use_triangle_normal=True
)

output_path = base_dir / "D_A.xyzn"
o3d.io.write_point_cloud(str(output_path), pcd)
