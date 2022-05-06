#include <catch2/catch.hpp>

#include <cstdlib>

#include "libslic3r/ExtrusionEntityCollection.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/ShortestPath.hpp"
#include "libslic3r/libslic3r.h"

#include "test_data.hpp"

using namespace Slic3r;

static inline Slic3r::Point random_point(float LO=-50, float HI=50) 
{
    Vec2f pt = Vec2f(LO, LO) + (Vec2d(rand(), rand()) * (HI-LO) / RAND_MAX).cast<float>();
	return pt.cast<coord_t>();
}

// build a sample extrusion entity collection with random start and end points.
static Slic3r::ExtrusionPath random_path(size_t length = 20, float LO = -50, float HI = 50)
{
    ExtrusionPath t {erPerimeter, 1.0, 1.0, 1.0};
    for (size_t j = 0; j < length; ++ j)
        t.polyline.append(random_point(LO, HI));
    return t;
}

static Slic3r::ExtrusionPaths random_paths(size_t count = 10, size_t length = 20, float LO = -50, float HI = 50)
{
    Slic3r::ExtrusionPaths p;
    for (size_t i = 0; i < count; ++ i)
        p.push_back(random_path(length, LO, HI));
    return p;
}

SCENARIO("ExtrusionEntityCollection: Polygon flattening", "[ExtrusionEntity]") {
    srand(0xDEADBEEF); // consistent seed for test reproducibility.

    // Generate one specific random path set and save it for later comparison
    Slic3r::ExtrusionPaths nosort_path_set = random_paths();

    Slic3r::ExtrusionEntityCollection sub_nosort;
    sub_nosort.append(nosort_path_set);
    sub_nosort.no_sort = true;

    Slic3r::ExtrusionEntityCollection sub_sort;
    sub_sort.no_sort = false;
    sub_sort.append(random_paths());

    GIVEN("A Extrusion Entity Collection with a child that has one child that is marked as no-sort") {
        Slic3r::ExtrusionEntityCollection sample;
        Slic3r::ExtrusionEntityCollection output;

        sample.append(sub_sort);
        sample.append(sub_nosort);
        sample.append(sub_sort);

        WHEN("The EEC is flattened with default options (preserve_order=false)") {
			output = sample.flatten();
            THEN("The output EEC contains no Extrusion Entity Collections") {
                CHECK(std::count_if(output.entities.cbegin(), output.entities.cend(), [=](const ExtrusionEntity* e) {return e->is_collection();}) == 0);
            }
        }
        WHEN("The EEC is flattened with preservation (preserve_order=true)") {
			output = sample.flatten(true);
            THEN("The output EECs contains one EEC.") {
                CHECK(std::count_if(output.entities.cbegin(), output.entities.cend(), [=](const ExtrusionEntity* e) {return e->is_collection();}) == 1);
            }
            AND_THEN("The ordered EEC contains the same order of elements than the original") {
                // find the entity in the collection
                for (auto e : output.entities)
                    if (e->is_collection()) {
                        ExtrusionEntityCollection *temp = dynamic_cast<ExtrusionEntityCollection*>(e);
                        // check each Extrusion path against nosort_path_set to see if the first and last match the same
                        CHECK(nosort_path_set.size() == temp->entities.size());
                        for (size_t i = 0; i < nosort_path_set.size(); ++ i) {
                            CHECK(temp->entities[i]->first_point() == nosort_path_set[i].first_point());
                            CHECK(temp->entities[i]->last_point() == nosort_path_set[i].last_point());
                        }
                    }
            }
        }
    }
}

TEST_CASE("ExtrusionEntityCollection: Chained path", "[ExtrusionEntity]") {
    struct Test {
        Polylines unchained;
        Polylines chained;
        Point     initial_point;
    };
    std::vector<Test> tests { 
        {
            { 
                { {0,15}, {0,18}, {0,20} },
                { {0,10}, {0,8}, {0,5} }
            },
            {
                { {0,20}, {0,18}, {0,15} },
                { {0,10}, {0,8}, {0,5} }
            },
            { 0,30 }
        },
        {
            { 
                { {4,0}, {10,0}, {15,0} },
                { {10,5}, {15,5}, {20,5} }
            },
            {
                { {20,5}, {15,5}, {10,5} },
                { {15,0}, {10,0}, {4,0} }
            },
            { 30,0 }
        },
        {
            { 
                { {15,0}, {10,0}, {4,0} },
                { {10,5}, {15,5}, {20,5} }
            },
            {
                { {20,5}, {15,5}, {10,5} },
                { {15,0}, {10,0}, {4,0} }
            },
            { 30,0 }
        },
    };
    for (const Test &test : tests) {
        Polylines chained = chain_polylines(test.unchained, &test.initial_point);
        REQUIRE(chained == test.chained);
        ExtrusionEntityCollection unchained_extrusions;
        extrusion_entities_append_paths(unchained_extrusions.entities, test.unchained,
            erInternalInfill, 0., 0.4f, 0.3f);
        ExtrusionEntityCollection chained_extrusions = unchained_extrusions.chained_path_from(test.initial_point);
        REQUIRE(chained_extrusions.entities.size() == test.chained.size());
        for (size_t i = 0; i < chained_extrusions.entities.size(); ++ i) {
            const Points &p1 = test.chained[i].points;
            const Points &p2 = dynamic_cast<const ExtrusionPath*>(chained_extrusions.entities[i])->polyline.points;
            REQUIRE(p1 == p2);
        }
    }
}

TEST_CASE("ExtrusionEntityCollection: Chained path with no explicit starting point", "[ExtrusionEntity]") {
    auto polylines = Polylines { { { 0, 15 }, {0, 18}, {0, 20} }, { { 0, 10 }, {0, 8}, {0, 5} } };
    auto target    = Polylines { { {0, 5}, {0, 8}, { 0, 10 } }, { { 0, 15 }, {0, 18}, {0, 20} } };
    auto chained   = chain_polylines(polylines);
    REQUIRE(chained == target);
}
