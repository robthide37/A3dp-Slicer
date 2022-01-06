#include "FindReplace.hpp"
#include "../Utils.hpp"

#include <cctype> // isalpha

namespace Slic3r {

GCodeFindReplace::GCodeFindReplace(const PrintConfig &print_config)
{
    const std::vector<std::string> &subst = print_config.gcode_substitutions.values;

    if ((subst.size() % 3) != 0)
        throw RuntimeError("Invalid length of gcode_substitutions parameter");

    m_substitutions.reserve(subst.size() / 3);
    for (size_t i = 0; i < subst.size(); i += 3) {
        Substitution out;
        try {
            out.plain_pattern    = subst[i];
            out.format           = subst[i + 1];
            const std::string &params = subst[i + 2];
            out.regexp           = strchr(params.c_str(), 'r') != nullptr || strchr(params.c_str(), 'R') != nullptr;
            out.case_insensitive = strchr(params.c_str(), 'i') != nullptr || strchr(params.c_str(), 'I') != nullptr;
            out.whole_word       = strchr(params.c_str(), 'w') != nullptr || strchr(params.c_str(), 'W') != nullptr;
            if (out.regexp)
                out.regexp_pattern.assign(
                    out.whole_word ? 
                        std::string("\b") + out.plain_pattern + "\b" :
                        out.plain_pattern,
                    (out.case_insensitive ? boost::regex::icase : 0) | boost::regex::optimize);
        } catch (const std::exception &ex) {
            throw RuntimeError(std::string("Invalid gcode_substitutions parameter, failed to compile regular expression: ") + ex.what());
        }
        m_substitutions.emplace_back(std::move(out));
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

template<typename FindFn>
static void find_and_replace_whole_word(std::string &inout, const std::string &match, const std::string &replace, FindFn find_fn)
{
    if (! match.empty() && inout.size() >= match.size() && match != replace) {
        std::string out;
        auto [i, j] = find_fn(inout, 0, match);
        size_t k = 0;
        for (; i != std::string::npos; std::tie(i, j) = find_fn(inout, i, match)) {
            if ((i == 0 || ! std::isalnum(inout[i - 1])) && (j == inout.size() || ! std::isalnum(inout[j]))) {
                out.reserve(inout.size());
                out.append(inout, k, i - k);
                out.append(replace);
                i = k = j;
            } else
                i += match.size();
        }
        if (k > 0) {
            out.append(inout, k, inout.size() - k);
            inout.swap(out);
        }
    }
}

std::string GCodeFindReplace::process_layer(const std::string &ain)
{
    std::string out;
    const std::string *in = &ain;
    std::string temp;
    temp.reserve(in->size());

    for (const Substitution &substitution : m_substitutions) {
        if (substitution.regexp) {
            temp.clear();
            temp.reserve(in->size());
            boost::regex_replace(ToStringIterator(temp), in->begin(), in->end(),
                substitution.regexp_pattern, substitution.format, boost::match_default | boost::match_not_dot_newline | boost::match_not_dot_null | boost::format_all);
            std::swap(out, temp);
        } else {
            if (in == &ain)
                out = ain;
            // Plain substitution
            if (substitution.case_insensitive) {
                if (substitution.whole_word)
                    find_and_replace_whole_word(out, substitution.plain_pattern, substitution.format,
                        [](const std::string &str, size_t start_pos, const std::string &match) {
                            auto begin = str.begin() + start_pos;
                            boost::iterator_range<std::string::const_iterator> r1(begin, str.end());
                            boost::iterator_range<std::string::const_iterator> r2(match.begin(), match.end());
                            auto res = boost::ifind_first(r1, r2);
                            return res ? std::make_pair(size_t(res.begin() - begin), size_t(res.end() - begin)) : std::make_pair(std::string::npos, std::string::npos);
                        });
                else
                    boost::ireplace_all(out, substitution.plain_pattern, substitution.format);
            } else {
                if (substitution.whole_word)
                    find_and_replace_whole_word(out, substitution.plain_pattern, substitution.format,
                        [](const std::string &str, size_t start_pos, const std::string &match) { 
                            size_t pos = str.find(match, start_pos);
                            return std::make_pair(pos, pos + (pos == std::string::npos ? 0 : match.size()));
                        });
                else
                    boost::replace_all(out, substitution.plain_pattern, substitution.format);
            }
        }
        in = &out;
    }

    return out;
}

}
