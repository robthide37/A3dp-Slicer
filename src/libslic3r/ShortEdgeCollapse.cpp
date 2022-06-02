#include "ShortEdgeCollapse.hpp"

#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>

namespace Slic3r {

void its_short_edge_collpase(indexed_triangle_set &mesh, size_t target_triangle_count) {
    std::unordered_map<size_t, size_t> vertices_index_mapping; // this map has two usages:
            // in the first step, it contains mapping from original vertex index to final vertex index (which is different from original only for removed vertices)
            // in the second step, it keeps the value for the removed vertices from the first step, but for those NOT removed, it contains mapping to the new position of the new vertices vector.
    std::unordered_set<size_t> vertices_to_remove;
    std::unordered_set<size_t> faces_to_remove;
    std::vector<Vec3i> triangles_neighbors;
    std::vector<Vec3f> face_normals;
    std::vector<Vec3f> vertex_normals; // Vertex normal in this algorithm is normal of any! triangle that contains the vertex

    std::vector<size_t> face_indices; //vector if indices, serves only the purpose of randomised traversal of the faces
    std::random_device rd;
    std::mt19937 generator(rd());

    float decimation_ratio = 1.0f; // decimation ratio updated in each iteration. it is number of removed triangles / number of all
    float edge_size = 0.2f; // Allowed collapsible edge size. Starts low, but is gradually increased
    size_t triangle_count = mesh.indices.size();

    std::vector<Vec3f> new_vertices;
    std::vector<Vec3i> new_indices;


    while (triangle_count > target_triangle_count) {
        // the following two switches try to keep the decimation ratio within 0.05 and 0.4 range (but nothing too clever :)
        if (decimation_ratio < 0.4) { // if decimation ratio is not too large, increase it
            edge_size *= 1.5f;
        }
        if (decimation_ratio < 0.05) { // if decimation ratio is very low, increase it once more
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
            //assign the normal of the current triangle to all its vertices. Do not care if its overridden, normal of any triangle will work
            // (so in the end, each vertex has assigned normal of one of its triangles with largest face index)
            vertex_normals[t[0]] = n;
            vertex_normals[t[1]] = n;
            vertex_normals[t[2]] = n;

            face_indices[face_idx] = face_idx; // just an index vector, will be later shuffled for randomized face traversal
        }

        std::vector<float> min_vertex_dot_product(mesh.vertices.size(), 1); // now compute vertices dot product - this is used during edge collapse,
        // to determine which vertex to remove and which to keep;  We try to keep the one with larger angle, because it defines the shape "more".
        // The min vertex dot product is lowest dot product of its normal (assigned in previous step) with the normals of faces around it.
        // the lower the dot product, the more we want to keep the vertex
        for (size_t face_idx = 0; face_idx < mesh.indices.size(); ++face_idx) {
            Vec3i t = mesh.indices[face_idx];
            Vec3f n = face_normals[face_idx];
            min_vertex_dot_product[t[0]] = std::min(min_vertex_dot_product[t[0]], n.dot(vertex_normals[t[0]]));
            min_vertex_dot_product[t[1]] = std::min(min_vertex_dot_product[t[1]], n.dot(vertex_normals[t[1]]));
            min_vertex_dot_product[t[2]] = std::min(min_vertex_dot_product[t[2]], n.dot(vertex_normals[t[2]]));
        }

        // prepare vertices mapping, intitalize with self index, no vertex has been removed yet
        for (size_t vertex_index = 0; vertex_index < mesh.vertices.size(); ++vertex_index) {
            vertices_index_mapping[vertex_index] = vertex_index;
        }

        //shuffle the faces and traverse in random order, this MASSIVELY improves the quality of the result
        std::shuffle(face_indices.begin(), face_indices.end(), generator);

        for (const size_t &face_idx : face_indices) {
            if (faces_to_remove.find(face_idx) != faces_to_remove.end()) {
            // if face already removed from previous collapses, skip (each collapse removes two triangles [at least] )
                continue;
            }
                // look at each edge if it is good candidate for collapse
            for (size_t edge_idx = 0; edge_idx < 3; ++edge_idx) {
                size_t vertex_index_keep = mesh.indices[face_idx][edge_idx];
                size_t vertex_index_remove = mesh.indices[face_idx][(edge_idx + 1) % 3];
                //check distance, skip long edges
                if ((mesh.vertices[vertex_index_keep] - mesh.vertices[vertex_index_remove]).squaredNorm()
                        > max_edge_len_squared) {
                    continue;
                }
                // swap indexes if vertex_index_keep has higher dot product (we want to keep low dot product vertices)
                if (min_vertex_dot_product[vertex_index_remove] < min_vertex_dot_product[vertex_index_keep]) {
                    size_t tmp = vertex_index_keep;
                    vertex_index_keep = vertex_index_remove;
                    vertex_index_remove = tmp;
                }

                // If vertex has already been part of a collapse, skip it (mark value -2.0 is set after collapse)
                if (min_vertex_dot_product[vertex_index_remove] < -1.1f) {
                    continue;
                }
                // if any of the collapsed vertices is already marked for removal, skip as well. (break, since no edge is collapsible in this case)
                if (vertices_to_remove.find(vertex_index_keep) != vertices_to_remove.end()
                        || vertices_to_remove.find(vertex_index_remove) != vertices_to_remove.end()) {
                    break;
                }

                // find neighbour triangle over this edge. if marked for removal, skip
                int neighbor_face_idx = triangles_neighbors[face_idx][edge_idx];
                if (neighbor_face_idx > 0 && faces_to_remove.find(neighbor_face_idx) != faces_to_remove.end()) {
                    continue;
                }

                //finaly do the collapse:
                //mark faces for removal
                faces_to_remove.insert(face_idx);
                faces_to_remove.insert(neighbor_face_idx);
                // mark one vertex for removal (with higher dot product)
                vertices_to_remove.insert(vertex_index_remove);
                // map its index to the index of the kept vertex
                vertices_index_mapping[vertex_index_remove] = vertices_index_mapping[vertex_index_keep];
                // mark the kept vertex, so that it cannot participate as a vertex for removal in any other collapse
                min_vertex_dot_product[vertex_index_keep] = -2.0f;
                // break, we are done with this triangle
                break;
            }

        }

        new_vertices.clear();
        new_vertices.reserve(mesh.vertices.size() - vertices_to_remove.size());
        for (size_t vertex_index = 0; vertex_index < mesh.vertices.size(); ++vertex_index) {
            // filter out removed vertices, add those NOT removed to the new_vertices vector, and also store their new position in the index mapping
            if (vertices_to_remove.find(vertex_index) == vertices_to_remove.end()) { //not removed
                new_vertices.push_back(mesh.vertices[vertex_index]);
                assert(vertices_index_mapping[vertex_index] == vertex_index);
                vertices_index_mapping[vertex_index] = new_vertices.size() - 1;
            }
        }

        new_indices.clear();
        new_indices.reserve(mesh.indices.size() - faces_to_remove.size());
        for (size_t face_idx = 0; face_idx < mesh.indices.size(); ++face_idx) {
            if (faces_to_remove.find(face_idx) != faces_to_remove.end()) {
                continue; //skip removed triangles
            }
            Vec3i new_triangle;
            for (int t_vertex_idx = 0; t_vertex_idx < 3; ++t_vertex_idx) {
                size_t orig_index = mesh.indices[face_idx][t_vertex_idx];
                if (vertices_to_remove.find(orig_index) == vertices_to_remove.end()) { //this vertex was not removed
                    // NOT removed vertices point to their new position, so use directly the mapping
                    new_triangle[t_vertex_idx] = vertices_index_mapping[orig_index];
                } else {
                    // Removed vertices point to the index of the kept vertex to which they collapsed.
                    // This vertex was surely not removed, so use the mapping twice (first to kept vertex, and then to its new position)
                    new_triangle[t_vertex_idx] = vertices_index_mapping[vertices_index_mapping[orig_index]];
                }
            }
            if (new_triangle[0] == new_triangle[1] || new_triangle[1] == new_triangle[2]
                    || new_triangle[2] == new_triangle[0]) {
                continue; //skip degenerate triangles
            }
            new_indices.push_back(new_triangle);
        }

        decimation_ratio = float(faces_to_remove.size()) / float(mesh.indices.size());

//        std::cout << " DECIMATION RATIO: " << decimation_ratio << std::endl;

        // finally update the mesh
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
