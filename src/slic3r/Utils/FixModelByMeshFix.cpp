#include "tmesh.h"

#include "FixModelByMeshFix.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Model.hpp"

#include <string.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <condition_variable>
#include <exception>
#include <string>
#include <thread>

#include <boost/filesystem.hpp>
#include <boost/nowide/convert.hpp>
#include <boost/nowide/cstdio.hpp>
#include <boost/thread.hpp>

#include "libslic3r/Model.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/Format/3mf.hpp"
#include "../GUI/GUI.hpp"
#include "../GUI/I18N.hpp"
#include "../GUI/MsgDialog.hpp"

#include <wx/msgdlg.h>
#include <wx/progdlg.h>

namespace Slic3r {

class RepairCanceledException: public std::exception {
public:
    const char* what() const throw () {
        return "Model repair has been cancelled";
    }
};

class RepairFailedException: public std::exception {
public:
    const char* what() const throw () {
        return "Model repair has failed";
    }
};

namespace detail {

using namespace T_MESH;

double closestPair(List *bl1, List *bl2, Vertex **closest_on_bl1, Vertex **closest_on_bl2)
        {
    Node *n, *m;
    Vertex *v, *w;
    double adist, mindist = DBL_MAX;

    FOREACHVVVERTEX(bl1, v, n)
        FOREACHVVVERTEX(bl2, w, m)
            if ((adist = w->squaredDistance(v)) < mindist)
                    {
                mindist = adist;
                *closest_on_bl1 = v;
                *closest_on_bl2 = w;
            }

    return mindist;
}

bool joinClosestComponents(Basic_TMesh *tin)
        {
    Vertex *v, *w, *gv, *gw;
    Triangle *t, *s;
    Node *n;
    List triList, boundary_loops, *one_loop;
    List **bloops_array;
    int i, j, numloops;

    // Mark triangles with connected component's unique ID
    i = 0;
    FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
    FOREACHVTTRIANGLE((&(tin->T)), t, n) if (t->info == NULL)
    {
        i++;
        triList.appendHead(t);
        t->info = (void *)(intptr_t)i;

        while (triList.numels())
        {
            t = (Triangle *)triList.popHead();
            if ((s = t->t1()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)(intptr_t)i;}
            if ((s = t->t2()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)(intptr_t)i;}
            if ((s = t->t3()) != NULL && s->info == NULL) {triList.appendHead(s); s->info = (void *)(intptr_t)i;}
        }
    }

    if (i<2)
    {
        FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;
        //   JMesh::info("Mesh is a single component. Nothing done.");
        return false;
    }

    FOREACHVTTRIANGLE((&(tin->T)), t, n)
    {
        t->v1()->info = t->v2()->info = t->v3()->info = t->info;
    }

    FOREACHVVVERTEX((&(tin->V)), v, n) if (!IS_VISITED2(v) && v->isOnBoundary())
    {
        w = v;
        one_loop = new List;
        do
        {
            one_loop->appendHead(w); MARK_VISIT2(w);
            w = w->nextOnBoundary();
        }while (w != v);
        boundary_loops.appendHead(one_loop);
    }
    FOREACHVVVERTEX((&(tin->V)), v, n) UNMARK_VISIT2(v);

    bloops_array = (List **)boundary_loops.toArray();
    numloops = boundary_loops.numels();

    double adist,
    mindist = DBL_MAX;

    gv = NULL;
    for (i = 0; i < numloops; i++)
        for (j = 0; j < numloops; j++)
            if (((Vertex*) bloops_array[i]->head()->data)->info != ((Vertex*) bloops_array[j]->head()->data)->info)
                    {
                adist = closestPair(bloops_array[i], bloops_array[j], &v, &w);
                if (adist < mindist) {
                    mindist = adist;
                    gv = v;
                    gw = w;
                }
            }

    if (gv != NULL)
        tin->joinBoundaryLoops(gv, gw, 1, 0);

    FOREACHVTTRIANGLE((&(tin->T)), t, n)
        t->info = NULL;
    FOREACHVVVERTEX((&(tin->V)), v, n)
        v->info = NULL;

    free(bloops_array);
    while ((one_loop = (List*) boundary_loops.popHead()) != NULL)
        delete one_loop;

    return (gv != NULL);
}

class Basic_TMesh_Adapter: public Basic_TMesh {
public:
    void load_indexed_triangle_set(const indexed_triangle_set &its) {

        for (const auto &vertex : its.vertices) {
            this->V.appendTail(this->newVertex(vertex.x(), vertex.y(), vertex.z()));
        }

        int nv = this->V.numels();
        Node *n = this->V.head();
        ExtVertex **tmp = (ExtVertex**) malloc(sizeof(ExtVertex*) * nv);
        for (int index = 0; index < nv; ++index) {
            tmp[index] = new ExtVertex(static_cast<Vertex*>(n->data));
            n = n->next();
        }

        for (const auto &face : its.indices) {
            if (face.x() == face.y() || face.y() == face.z() || face.z() == face.x()) {
                continue;
            }
            this->CreateIndexedTriangle(tmp, face.x(), face.y(), face.z());
        }

        closeLoadingSession(this->T.numels(), tmp, false);
    }

    indexed_triangle_set to_indexed_triangle_set() {
        indexed_triangle_set out;
        out.vertices.resize(this->V.numels());
        out.indices.resize(this->T.numels());

        Node *n = this->V.head();
        for (int vertex_idx = 0; vertex_idx < this->V.numels(); ++vertex_idx) {
            T_MESH::Vertex *v = static_cast<T_MESH::Vertex*>(n->data);
            out.vertices[vertex_idx] = Vec3f { float(v->x), float(v->y), float(v->z) };
            n = n->next();
            //store vertex index in the first coord, makes export of faces simple (inspired by saveOBJ method)
            v->x = vertex_idx;
        }

        n = this->T.head();
        for (int face_idx = 0; face_idx < this->T.numels(); ++face_idx) {
            T_MESH::Triangle *t = static_cast<T_MESH::Triangle*>(n->data);
            out.indices[face_idx] = Vec3i { int(t->v1()->x), int(t->v2()->x), int(t->v3()->x) };
            n = n->next();
        }

        return out;
    }

    bool meshclean_single_iteration(int inner_loops, const std::atomic<bool> &canceled) {
        bool ni, nd;
        Triangle *t;
        Node *m;

        nd = strongDegeneracyRemoval(inner_loops);
        if (canceled)
            throw RepairCanceledException();
        deselectTriangles();
        invertSelection();
        ni = strongIntersectionRemoval(inner_loops);
        if (canceled)
            throw RepairCanceledException();
        if (ni && nd) {
            FOREACHTRIANGLE(t, m)
                if (t->isExactlyDegenerate())
                    ni = false;
            if (ni)
                return true;
        }

        return false;
    }

private:

    void closeLoadingSession(int loaded_faces, ExtVertex **tmp, bool triangulate) {
        int i, nv = this->V.numels();
        if (tmp != NULL)
        {
            for (i = 0; i < nv; i++)
                delete (tmp[i]);
            free(tmp);
        }
        fixConnectivity();
        this->d_boundaries = this->d_handles = this->d_shells = 1;
    }
};

}

bool fix_model_by_meshfix(ModelObject &model_object, int volume_idx, wxProgressDialog &progress_dlg,
        const wxString &msg_header, std::string &fix_result) {
    std::mutex mtx;
    std::condition_variable condition;
    struct Progress {
        std::string message;
        int percent = 0;
        bool updated = false;
    } progress;
    std::atomic<bool> canceled = false;
    std::atomic<bool> finished = false;

    std::vector<ModelVolume*> volumes;
    if (volume_idx == -1)
        volumes = model_object.volumes;
    else
        volumes.emplace_back(model_object.volumes[volume_idx]);

    // Executing the calculation in a background thread, so that the COM context could be created with its own threading model.
    // (It seems like wxWidgets initialize the COM contex as single threaded and we need a multi-threaded context).
    bool success = false;
    size_t ivolume = 0;
    auto on_progress = [&mtx, &condition, &ivolume, &volumes, &progress](const char *msg, unsigned prcnt) {
        std::unique_lock<std::mutex> lock(mtx);
        progress.message = msg;
        progress.percent = (int) floor((float(prcnt) + float(ivolume) * 100.f) / float(volumes.size()));
        progress.updated = true;
        condition.notify_all();
    };
    auto worker_thread = boost::thread(
            [&model_object, &volumes, &ivolume, on_progress, &success, &canceled, &finished]() {
                try {
                    std::vector<TriangleMesh> meshes_repaired;
                    meshes_repaired.reserve(volumes.size());
                    for (ModelVolume *mv : volumes) {
                        std::vector<indexed_triangle_set> parts = its_split(mv->mesh().its);
                        unsigned percent_per_part = 95 / parts.size();
                        for (size_t part_idx = 0; part_idx < parts.size(); ++part_idx) {
                            unsigned progress_part_base = part_idx * percent_per_part;

                            detail::Basic_TMesh_Adapter tin { };
                            on_progress(L("Loading source model"), progress_part_base);
                            if (canceled)
                                throw RepairCanceledException();
                            tin.load_indexed_triangle_set(parts[part_idx]);
                            tin.boundaries();
                            on_progress(L("Join closest components"), unsigned(progress_part_base + 0.1 * percent_per_part));
                            if (canceled)
                                throw RepairCanceledException();
                            joinClosestComponents(&tin);
                            tin.deselectTriangles();
                            tin.boundaries();
                            // Keep only the largest component (i.e. with most triangles)
                            on_progress(L("Remove smallest components"), unsigned(progress_part_base + 0.2 * percent_per_part));
                            if (canceled)
                                throw RepairCanceledException();
                            tin.removeSmallestComponents();

                            // Fill holes
                            on_progress(L("Check holes"), unsigned(progress_part_base + 0.3 * percent_per_part));
                            if (canceled)
                                throw RepairCanceledException();
                            if (tin.boundaries()) {
                                on_progress(L("Patch small holes"), unsigned(progress_part_base + 0.4 * percent_per_part));
                                if (canceled)
                                    throw RepairCanceledException();
                                tin.fillSmallBoundaries(0, true);
                            }

                            on_progress(L("Geometry check"), unsigned(progress_part_base + 0.5 * percent_per_part));
                            if (canceled)
                                throw RepairCanceledException();
                            // Run geometry correction
                            if (!tin.boundaries()) {
                                int iteration = 0;
                                on_progress(L("Start iterative correction"), unsigned(progress_part_base + 0.55 * percent_per_part));
                                tin.deselectTriangles();
                                tin.invertSelection();
                                bool fixed = false;
                                while (iteration < 10 && !fixed) { //default constants taken from TMesh library
                                    fixed = tin.meshclean_single_iteration(3, canceled);
                                    on_progress(L("Fixing geometry"), progress_part_base + percent_per_part * std::min(0.9, 0.6 + iteration*0.08)); // majority of objects should finish in 4 iterations
                                    if (canceled)
                                        throw RepairCanceledException();
                                    iteration++;
                                }
                            }

                            if (tin.boundaries() || tin.T.numels() == 0) {
                                throw RepairFailedException();
                            }
                            parts[part_idx] = tin.to_indexed_triangle_set();
                        }

                        for (size_t part_idx = 1; part_idx < parts.size(); ++part_idx) {
                            its_merge(parts[0], parts[part_idx]);
                        }

                        meshes_repaired.emplace_back(std::move(parts[0]));
                    }

                    for (size_t i = 0; i < volumes.size(); ++i) {
                        volumes[i]->set_mesh(std::move(meshes_repaired[i]));
                        volumes[i]->calculate_convex_hull();
                        volumes[i]->set_new_unique_id();
                    }
                    model_object.invalidate_bounding_box();
                    --ivolume;
                    on_progress(L("Model repair finished"), 100);
                    success = true;
                    finished = true;
                } catch (RepairCanceledException& /* ex */) {
                    canceled = true;
                    finished = true;
                    on_progress(L("Model repair canceled"), 100);
                } catch (std::exception &ex) {
                    success = false;
                    finished = true;
                    on_progress(ex.what(), 100);
                }
            });

    while (!finished) {
        std::unique_lock<std::mutex> lock(mtx);
        condition.wait_for(lock, std::chrono::milliseconds(250), [&progress] {
            return progress.updated;
        });
        // decrease progress.percent value to avoid closing of the progress dialog
        if (!progress_dlg.Update(progress.percent - 1, msg_header + _(progress.message)))
            canceled = true;
        else
            progress_dlg.Fit();
        progress.updated = false;
    }

    worker_thread.join();

    if (canceled) {
        // Nothing to show.
    } else if (success) {
        fix_result = "";
    } else {
        fix_result = progress.message;
    }

    return !canceled;
}

}
