import gzip
import io
import math
import os
import sys
from pathlib import Path

import bpy
import ezdxf
import pybase64
import shapely

from shapely import Polygon, validation
from ezdxf.addons import odafc
from src.ezdxf.acis import api

dwg_file_name = "output.dwg"
dxf_file_name = "output.dxf"
obj_dir_name = "obj_files"
gltf_file_name = "output.glb"
addon_name_import_dxf = "io_import_dxf"

new_addons_path = "/blender-4.1.1-linux-x64/4.1/scripts/addons"
if new_addons_path not in sys.path:
    sys.path.append(new_addons_path)
new_python_path = "/blender-4.1.1-linux-x64/4.1/python/lib/python3.11/site-packages"
if new_python_path not in sys.path:
    sys.path.append(new_python_path)

import addon_utils


# def lambda_handler(event, context):
#     # Extract the data from the event
#     data = event['data']
#
#     # Call the main function
#     data_decoded = decode_and_unzip_base64_gzip(data, dwg_file_name)
#     result = convert_dwg_to_gltf(data_decoded)
#     result_encoded = encode_file_to_base64(result)
#
#     # Return a response
#     return {
#         'statusCode': 200,
#         'body': result_encoded
#     }


def encode_file_to_base64(file_path: str) -> bytes:
    with open(file_path, 'rb') as file:
        file_data = file.read()
    return pybase64.b64encode(file_data)


def decode_and_unzip_base64_gzip(encoded_data: bytes, output_filename: str):
    # Decode the base64 encoded data
    decoded_data = pybase64.b64decode(encoded_data)

    # Unzip the gzip data
    with gzip.GzipFile(fileobj=io.BytesIO(decoded_data)) as gzip_file:
        unzipped_data = gzip_file.read()

    # Save the unzipped data to a file in the current directory
    output_path = os.path.join(os.getcwd(), output_filename)
    with open(output_path, 'wb') as output_file:
        output_file.write(unzipped_data)
    return output_filename


def convert_dwg_to_dxf(filename: str):
    try:
        output_filename = Path(filename).with_suffix('.dxf').name

        ezdxf.options.read_file("./config.ini")
        ezdxf.options.update_cached_options()

        odafc.convert(source=filename, dest=output_filename, replace=True)
    except Exception as e:
        raise Exception(f"Failed to convert CAD data: {e}")


def convert_dwg_to_gltf(filename: str) -> str:
    # convert dwg to dxf
    # convert_dwg_to_dxf(filename=filename)

    # import dxf to blender
    blender_addon_import_dxf = "bl_ext.blender_org.import_autocad_dxf_format_dxf"
    addon_utils.enable(blender_addon_import_dxf)
    bpy.ops.import_scene.dxf(files=[{"name": dxf_file_name}])

    # import obj's
    extract_3dsolid_with_position_to_obj(dxf_file_name)
    obj_files = [f for f in os.listdir(obj_dir_name) if f.endswith('.obj')]
    for file in obj_files:
        bpy.ops.wm.obj_import(filepath=f"{obj_dir_name}/{file}", forward_axis='Y', up_axis='Z')

    # Export to gltf
    bpy.ops.export_scene.gltf(filepath=gltf_file_name)

    return gltf_file_name


# 메쉬의 구멍(holes)을 기준으로 면(faces)을 분할
def split_faces_by_holes(faces: list[tuple[int]], holes: list[tuple[int]], verts: list[tuple[float]]) -> list:
    # hole이 있는 face와 해당 hole을 매칭
    matching_pairs = [(face, hole) for face, hole in zip(faces, holes) if hole]
    # hole이 없는 face를 필터링
    faces_without_holes = [face for face, hole in zip(faces, holes) if not hole]

    # face와 hole의 인덱스를 좌표로 변환
    def convert_to_coordinates(pairs, vertices):
        converted_pairs = []
        for face, hole in pairs:
            face_coords = [vertices[idx] for idx in face]
            hole_coords = [vertices[idx] for idx in hole]
            converted_pairs.append((face_coords, hole_coords))
        return converted_pairs

    face_to_hole_as_verts = convert_to_coordinates(matching_pairs, verts)

    # 2d 좌표로 변환
    def remove_constant_axis(coords):
        if all(coord[0] == coords[0][0] for coord in coords):
            return [(coord[1], coord[2]) for coord in coords], 0, coords[0][0]
        elif all(coord[1] == coords[0][1] for coord in coords):
            return [(coord[0], coord[2]) for coord in coords], 1, coords[0][1]
        elif all(coord[2] == coords[0][2] for coord in coords):
            return [(coord[0], coord[1]) for coord in coords], 2, coords[0][2]
        return coords, -1, None  # In case no axis is constant, return original

    def add_constant_axis(coords, axis, value):
        if axis == 0:
            return [(value, coord[0], coord[1]) for coord in coords]
        elif axis == 1:
            return [(coord[0], value, coord[1]) for coord in coords]
        elif axis == 2:
            return [(coord[0], coord[1], value) for coord in coords]
        return coords  # In case no axis is constant, return original

    polygons = []
    triangle_3ds = []
    for face_coords, hole_coords in face_to_hole_as_verts:
        converted_2d_pairs = []
        face_2d, axis, value = remove_constant_axis(face_coords)
        hole_2d, _, _ = remove_constant_axis(hole_coords)
        converted_2d_pairs.append((face_2d, hole_2d))

        # create shapely polygon with faces in converted_2d_pairs
        for face, hole in converted_2d_pairs:
            polygon = Polygon(face, [hole])
            if not polygon.is_valid:
                print(f"Invalid polygon: {validation.explain_validity(polygon)}")
                continue
            polygons.append(polygon)
            triangles = [tri for tri in shapely.delaunay_triangles(polygon).geoms if tri.within(polygon)]

            # Revert 2D triangles back to 3D coordinates
            triangles_3d = []
            for tri in triangles:
                tri_coords = list(tri.exterior.coords)[:-1]  # Remove the closing coordinate
                tri_3d = add_constant_axis(tri_coords, axis, value)
                triangles_3d.append(tri_3d)

            triangle_3ds.extend(triangles_3d)

    # Create a mapping from vertex coordinates to their indices
    vert_to_index = {tuple(v): i for i, v in enumerate(verts)}

    # Convert triangles_3d to tuples of indices
    triangles_as_indices = []
    for tri in triangle_3ds:
        tri_indices = tuple(vert_to_index[tuple(coord)] for coord in tri)
        triangles_as_indices.append(tri_indices)

    all_faces = triangles_as_indices + faces_without_holes
    return all_faces


def extract_3dsolid_with_position_to_obj(dxf_path):
    """DXF에서 3DSOLID를 추출하고 위치(Translation) + 스케일(Scale) + 회전(Rotation) 적용 후 OBJ 변환"""
    try:
        doc = ezdxf.readfile(dxf_path)
    except Exception as e:
        print(f"❌ DXF 파일 로드 실패: {e}")
        return

    msp = doc.modelspace()
    solid_count = 0

    for entity in msp:
        if entity.dxftype() == "INSERT":
            # 🔹 블록 참조(INSERT)의 위치 및 변환 정보 가져오기
            block_name = entity.dxf.name
            translation = (entity.dxf.insert.x, entity.dxf.insert.y, entity.dxf.insert.z)
            scale = (entity.dxf.xscale, entity.dxf.yscale, entity.dxf.zscale)
            rotation = entity.dxf.rotation  # AutoCAD 기준, degrees

            # 🔹 블록 내부의 3DSOLID 객체 찾기
            try:
                block = doc.blocks.get(block_name)
                for solid in block.query("3DSOLID"):
                    bodies = []
                    try:
                        bodies = api.load_dxf(solid)  # ACIS 데이터 로드
                    except Exception as e:
                        print(f"⚠️ 3DSOLID 로드 실패 (Handle: {solid.dxf.handle}) - {e}")
                        continue

                    for body in bodies:
                        mesh_list = []
                        try:
                            mesh_list = api.mesh_from_body(body, merge_lumps=True)  # ACIS → 메쉬 변환
                        except Exception as e:
                            print(f"⚠️ ACIS -> 메쉬 변환 실패 (Handle: {solid.dxf.handle}) - {e}")
                            continue

                        for mesh in mesh_list:
                            mesh.normalize_faces()
                            verts = [tuple(v) for v in mesh.vertices]
                            faces = [tuple(face) for face in mesh.faces]
                            holes = mesh.holes

                            if not all(hole == () for hole in holes):
                                print("Hole Detected: solid count = ", solid_count)
                                faces = split_faces_by_holes(faces, holes, verts)

                            if verts and faces:
                                # 🔹 블록의 위치, 스케일, 회전 적용
                                transformed_verts = apply_transformation(verts, translation, scale, rotation)
                                obj_filename = f"3DSOLID_{solid_count}.obj"
                                save_obj(obj_filename, transformed_verts, faces)
                                solid_count += 1
            except KeyError:
                print(f"⚠️ 블록 '{block_name}'이 존재하지 않음.")

        elif entity.dxftype() == "3DSOLID":
            bodies = []
            try:
                bodies = api.load_dxf(entity)
            except Exception as e:
                print(f"⚠️ 3DSOLID 로드 실패 (Handle: {entity.dxf.handle}) - {e}")
                continue

            translation = (0, 0, 0)  # 기본값
            scale = (1, 1, 1)
            rotation = 0

            try:
                if hasattr(entity.dxf, "insert"):
                    translation = (entity.dxf.insert.x, entity.dxf.insert.y, entity.dxf.insert.z)
                elif hasattr(entity.dxf, "location"):
                    translation = (entity.dxf.location.x, entity.dxf.location.y, entity.dxf.location.z)
            except AttributeError:
                pass

            print(f"🔹 3DSOLID 변환 중... 위치: {translation}, 스케일: {scale}, 회전: {rotation}°")

            for body in bodies:
                mesh_list = []
                try:
                    mesh_list = api.mesh_from_body(body, merge_lumps=True)
                except Exception as e:
                    print(f"⚠️ ACIS -> 메쉬 변환 실패 (Handle: {entity.dxf.handle}) - {e}")
                    continue

                for mesh in mesh_list:
                    mesh.normalize_faces()
                    verts = [tuple(v) for v in mesh.vertices]
                    faces = [tuple(face) for face in mesh.faces]
                    holes = mesh.holes

                    if not all(hole == () for hole in holes):
                        print("Hole Detected: solid count = ", solid_count)
                        faces = split_faces_by_holes(faces, holes, verts)
                    if verts and faces:
                        transformed_verts = apply_transformation(verts, translation, scale, rotation)
                        obj_filename = f"3DSOLID_{solid_count}.obj"
                        save_obj(obj_filename, transformed_verts, faces)
                        solid_count += 1

    if solid_count > 0:
        print(f"✅ {solid_count} 개의 3DSOLID 객체가 OBJ로 변환되었습니다.")
    else:
        print("⚠️ 변환할 3DSOLID 객체가 없습니다.")


def save_obj(filename, vertices, faces):
    """OBJ 파일을 저장 (위치, 스케일, 회전 반영)"""
    obj_folder = "obj_files"
    file_path = os.path.join(obj_folder, filename)

    with open(file_path, "w") as obj_file:
        obj_file.write("# Exported OBJ from DXF 3DSOLID (Preserve Position, Scale, Rotation)\n")
        for v in vertices:
            obj_file.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for f in faces:
            obj_file.write("f " + " ".join(str(i + 1) for i in f) + "\n")


def apply_transformation(vertices, translation, scale, rotation):
    """위치, 스케일, 회전을 적용하여 원래 위치 유지"""
    transformed_vertices = []

    # 회전 변환 (AutoCAD는 Z축 회전을 사용)
    angle = math.radians(rotation)  # AutoCAD는 degree 단위 → radian 변환
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)

    for v in vertices:
        # 스케일 적용
        x, y, z = v[0] * scale[0], v[1] * scale[1], v[2] * scale[2]

        # Z축 회전 적용
        x_rot = x * cos_a - y * sin_a
        y_rot = x * sin_a + y * cos_a

        # 최종 위치 변환 적용
        transformed_vertices.append((x_rot + translation[0], y_rot + translation[1], z + translation[2]))

    return transformed_vertices


result = convert_dwg_to_gltf(dwg_file_name)
result_encoded = encode_file_to_base64(result)
