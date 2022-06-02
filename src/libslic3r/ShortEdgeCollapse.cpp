#include "ShortEdgeCollapse.hpp"

#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>

namespace Slic3r {

void its_short_edge_collpase(indexed_triangle_set &mesh, size_t target_triangle_count) {
    std::unordered_map<size_t, size_t> vertices_index_mapping; // this map has two usages:
            // in the first step, it contains mapping from original vertex index to final vertex index (which is different from original only for removed vertices)

    std::unordered_set<size_t> vertices_to_remove;
    std::unordered_set<size_t> faces_to_remove;
    std::vector<Vec3i> triangles_neighbors;
    std::vector<Vec3f> face_normals;
    std::vector<Vec3f> vertex_normals;
    std::vector<size_t> face_indices;
    std::random_device rd;
    std::mt19937 generator(rd());

    float decimation_ratio = 1.0f;
    float edge_size = 0.2f;
    size_t triangle_count = mesh.indices.size();

    while (triangle_count > target_triangle_count) {
        if (decimation_ratio < 0.4) {
            edge_size *= 1.5f;
        }
        if (decimation_ratio < 0.05) {
            edge_size *= 1.5f;
        }

        float max_edge_len_squared = edge_size * edge_size;
        triangles_neighbors = its_face_neighbors_par(mesh);
        face_normals = its_face_normals(mesh);

        vertex_normals.resize(mesh.vertices.size());
        face_indices.resize(mesh.indices.size());
        for (size_t face_idx = 0; face_idx < mesh.indices.size(); ++face_idx) {
            Vec3i t = mesh.indices[face_idx];
            Vec3f n = face_normals[face_idx];
            vertex_normals[t[0]] = n;
            vertex_normals[t[1]] = n;
            vertex_normals[t[2]] = n;

            face_indices[face_idx] = face_idx;
        }

        std::vector<float> min_vertex_score(mesh.vertices.size(), 1);
        for (size_t face_idx = 0; face_idx < mesh.indices.size(); ++face_idx) {
            Vec3i t = mesh.indices[face_idx];
            Vec3f n = face_normals[face_idx];
            min_vertex_score[t[0]] = std::min(min_vertex_score[t[0]], n.dot(vertex_normals[t[0]]));
            min_vertex_score[t[1]] = std::min(min_vertex_score[t[1]], n.dot(vertex_normals[t[1]]));
            min_vertex_score[t[2]] = std::min(min_vertex_score[t[2]], n.dot(vertex_normals[t[2]]));
        }

        for (size_t vertex_index = 0; vertex_index < mesh.vertices.size(); ++vertex_index) {
            vertices_index_mapping[vertex_index] = vertex_index;
        }

        std::shuffle(face_indices.begin(), face_indices.end(), generator);

        for (const size_t &face_idx : face_indices) {
            if (faces_to_remove.find(face_idx) != faces_to_remove.end()) {
                continue;
            }
            for (size_t edge_idx = 0; edge_idx < 3; ++edge_idx) {
                size_t vertex_index_keep = mesh.indices[face_idx][edge_idx];
                size_t vertex_index_remove = mesh.indices[face_idx][(edge_idx + 1) % 3];

                if ((mesh.vertices[vertex_index_keep] - mesh.vertices[vertex_index_remove]).squaredNorm()
                        > max_edge_len_squared) {
                    continue;
                }

                if (min_vertex_score[vertex_index_remove] < min_vertex_score[vertex_index_keep]) {
                    size_t tmp = vertex_index_keep;
                    vertex_index_keep = vertex_index_remove;
                    vertex_index_remove = tmp;
                }

                if (min_vertex_score[vertex_index_remove] < -1.0f) {
                    continue;
                }

                if (vertices_to_remove.find(vertex_index_keep) != vertices_to_remove.end()
                        || vertices_to_remove.find(vertex_index_remove) != vertices_to_remove.end()) {
                    break;
                }

                int neighbor_face_idx = triangles_neighbors[face_idx][edge_idx];
                if (neighbor_face_idx > 0 && faces_to_remove.find(neighbor_face_idx) != faces_to_remove.end()) {
                    continue;
                }

                faces_to_remove.insert(face_idx);
                faces_to_remove.insert(neighbor_face_idx);
                vertices_to_remove.insert(vertex_index_remove);
                vertices_index_mapping[vertex_index_remove] = vertices_index_mapping[vertex_index_keep];
                min_vertex_score[vertex_index_keep] = -2.0;
                break;
            }

        }

        //flatten the mapping
        for (auto &pair : vertices_index_mapping) {
            while (vertices_index_mapping[pair.second] != pair.second) {
                pair.second = vertices_index_mapping[pair.second];
            }
        }

        std::vector<Vec3f> new_vertices;
        new_vertices.reserve(mesh.vertices.size() - vertices_to_remove.size());
        for (size_t vertex_index = 0; vertex_index < mesh.vertices.size(); ++vertex_index) {
            if (vertices_to_remove.find(vertex_index) == vertices_to_remove.end()) { //not removed
                new_vertices.push_back(mesh.vertices[vertex_index]);
                assert(vertices_index_mapping[vertex_index] == vertex_index);
                vertices_index_mapping[vertex_index] = new_vertices.size() - 1;
            }
        }

        std::vector<Vec3i> new_indices;
        new_indices.reserve(mesh.indices.size() - faces_to_remove.size());
        for (size_t face_idx = 0; face_idx < mesh.indices.size(); ++face_idx) {
            if (faces_to_remove.find(face_idx) != faces_to_remove.end()) {
                continue; //skip removed triangles
            }
            Vec3i new_triangle;
            for (int t_vertex_idx = 0; t_vertex_idx < 3; ++t_vertex_idx) {
                size_t orig_index = mesh.indices[face_idx][t_vertex_idx];
                if (vertices_to_remove.find(orig_index) == vertices_to_remove.end()) { //this vertex was not removed
                    new_triangle[t_vertex_idx] = vertices_index_mapping[orig_index];
                } else { // this vertex was removed, so use the vertex mapping points to
                    new_triangle[t_vertex_idx] = vertices_index_mapping[vertices_index_mapping[orig_index]];
                }
            }
            if (new_triangle[0] == new_triangle[1] || new_triangle[1] == new_triangle[2]
                    || new_triangle[2] == new_triangle[0]) {
                continue; //skip degenerate
            }
            new_indices.push_back(new_triangle);
        }

        decimation_ratio = float(faces_to_remove.size()) / float(mesh.indices.size());

//        std::cout << " DECIMATION RATIO: " << decimation_ratio << std::endl;

        mesh.vertices = new_vertices;
        mesh.indices = new_indices;

        vertices_index_mapping.clear();
        vertices_to_remove.clear();
        faces_to_remove.clear();
        triangles_neighbors.clear();
        face_normals.clear();
        vertex_normals.clear();
        face_indices.clear();
        triangle_count = mesh.indices.size();
    }
}

}
