#ifndef slic3r_FindReplace_hpp_
#define slic3r_FindReplace_hpp_

#include "../PrintConfig.hpp"

#include <boost/regex.hpp>

namespace Slic3r {

class GCodeFindReplace {
public:
    GCodeFindReplace(const PrintConfig &print_config);

    std::string process_layer(const std::string &gcode);
    
private:
    struct Substitution {
        std::string     plain_pattern;
        boost::regex    regexp_pattern;
        std::string     format;

        bool            regexp { false };
        bool            case_insensitive { false };
        bool            whole_word { false };
    };
    std::vector<Substitution> m_substitutions;
};

}

#endif // slic3r_FindReplace_hpp_
