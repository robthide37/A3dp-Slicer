#include "FindReplace.hpp"
#include "../Utils.hpp"

namespace Slic3r {

GCodeFindReplace::GCodeFindReplace(const PrintConfig &print_config)
{
    const std::vector<std::string> &subst = print_config.gcode_substitutions.values;

    if ((subst.size() % 3) != 0)
        throw RuntimeError("Invalid length of gcode_substitutions parameter");

    m_substitutions.reserve(subst.size() / 3);
    for (size_t i = 0; i < subst.size(); i += 3) {
        boost::regex pattern;
        try {
            pattern.assign(subst[i], boost::regex::optimize); // boost::regex::icase
        } catch (const std::exception &ex) {
            throw RuntimeError(std::string("Invalid gcode_substitutions parameter, failed to compile regular expression: ") + ex.what());
        }
        m_substitutions.push_back({ std::move(pattern), subst[i + 1] });
    }
}

class ToStringIterator 
{
public:
    using iterator_category     = std::output_iterator_tag;
    using value_type            = void;
    using difference_type       = void;
    using pointer               = void;
    using reference             = void;

    ToStringIterator(std::string &data) : m_data(&data) {}

    ToStringIterator& operator=(const char val) {
        size_t needs = m_data->size() + 1;
        if (m_data->capacity() < needs)
            m_data->reserve(next_highest_power_of_2(needs));
        m_data->push_back(val);
        return *this;
    }

    ToStringIterator& operator*()     { return *this; }
    ToStringIterator& operator++()    { return *this; }
    ToStringIterator  operator++(int) { return *this; }

private:
    std::string *m_data;
};

std::string GCodeFindReplace::process_layer(const std::string &ain)
{
    std::string out;
    const std::string *in = &ain;
    std::string temp;
    temp.reserve(in->size());

    for (const Substitution &substitution : m_substitutions) {
        temp.clear();
        temp.reserve(in->size());
        boost::regex_replace(ToStringIterator(temp), in->begin(), in->end(),
            substitution.pattern, substitution.format, boost::match_default | boost::match_not_dot_newline | boost::match_not_dot_null | boost::format_all);
        std::swap(out, temp);
        in = &out;
    }

    return out;
}

}
