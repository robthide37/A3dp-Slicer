#ifndef slic3r_MapUtils_hpp_
#define slic3r_MapUtils_hpp_

#include <map>
namespace Slic3r {

/// <summary>
/// Utility for easier work with standart map 
/// </summary>
class MapUtils
{
public:
    MapUtils() = delete; // only static functions

    /// <summary>
    /// Create map with swaped key-value
    /// </summary>
    /// <param name="map">Input map</param>
    /// <returns>Map with changed key to value and vice versa</returns>
    template<typename Key, typename Value>
    static std::map<Value, Key> create_oposit(const std::map<Key, Value> &map)
    {
        std::map<Value, Key> result;
        for (const auto &it : map) result[it.second] = it.first;
        return result;
    }
};
} // namespace Slic3r
#endif slic3r_MapUtils_hpp_