#ifndef slic3r_GCode_LabelObjects_hpp_
#define slic3r_GCode_LabelObjects_hpp_

#include "Print.hpp"

namespace Slic3r {

enum GCodeFlavor : unsigned char;
enum class LabelObjectsStyle;


namespace GCode {

    //std::string label_object_start(LabelObjectsStyle label_object_style, GCodeFlavor flavor, const SpanOfConstPtrs<PrintObject>& objects, int object_id, int instance_id);
    //std::string label_object_stop(LabelObjectsStyle label_object_style, GCodeFlavor flavor, int object_id, int instance_id, const std::string& name);
    //std::string label_all_objects(LabelObjectsStyle label_objects_style, GCodeFlavor flavor, const Print& print);


class LabelObjects {
public:
    enum class IncludeName {
        No,
        Yes
    };
    void init(const Print& print);
    std::string all_objects_header() const;
    std::string start_object(const PrintInstance& print_instance, IncludeName include_name) const;
    std::string stop_object(const PrintInstance& print_instance) const;

private:
    struct LabelData {
        std::string name;
        int unique_id;
        int object_id;
        int instance_id;
        bool object_has_more_instances;
    };

    LabelObjectsStyle m_label_objects_style;
    GCodeFlavor       m_flavor;
    std::unordered_map<const PrintInstance*, LabelData> m_label_data;

};


} // namespace GCode
} // namespace Slic3r

#endif // slic3r_GCode_LabelObjects_hpp_
