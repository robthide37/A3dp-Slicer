#include "FanMover.hpp"

#include "GCodeReader.hpp"
#include "LocalesUtils.hpp"

#include <iomanip>

#include <boost/log/trivial.hpp>

/*
#include <memory.h>
#include <string.h>
#include <float.h>

#include "../libslic3r.h"
#include "../PrintConfig.hpp"
#include "../Utils.hpp"
#include "Print.hpp"
*/


namespace Slic3r {

const std::string& FanMover::process_gcode(const std::string& gcode, bool flush)
{
    m_process_output = "";

    // recompute buffer time to recover from rounding
    m_buffer_time_size = 0;
    for (auto& data : m_buffer) m_buffer_time_size += data.time;

    if(!gcode.empty())
        m_parser.parse_buffer(gcode,
            [this](GCodeReader& reader, const GCodeReader::GCodeLine& line) { /*m_process_output += line.raw() + "\n";*/ this->_process_gcode_line(reader, line); });

    if (flush) {
        while (!m_buffer.empty()) {
            write_buffer_data();
        }
    }

    return m_process_output;
}

bool is_end_of_word(char c) {
   return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == 0;
}

float get_axis_value(const std::string& line, char axis)
{
    char match[3] = " X";
    match[1] = axis;

    size_t pos = line.find(match) + 2;
    //size_t end = std::min(line.find(' ', pos + 1), line.find(';', pos + 1));
    // Try to parse the numeric value.
    const char* c = line.c_str();
    char* pend = nullptr;
    errno = 0;
    double  v = strtod(c + pos, &pend);
    if (pend != nullptr && errno == 0 && pend != c) {
        // The axis value has been parsed correctly.
        return float(v);
    }
    return NAN;
}

void change_axis_value(std::string& line, char axis, const float new_value, const int decimal_digits)
{
    char match[3] = " X";
    match[1] = axis;

    size_t pos = line.find(match) + 2;
    size_t end = std::min(line.find(' ', pos + 1), line.find(';', pos + 1));
    line = line.replace(pos, end - pos, to_string_nozero(new_value, decimal_digits));
}

int16_t get_fan_speed(const std::string &line, GCodeFlavor flavor) {
    if (line.compare(0, 4, "M106") == 0) {
        if (flavor == (gcfMach3) || flavor == (gcfMachinekit)) {
            return (int16_t)get_axis_value(line, 'P');
        } else {
            return (int16_t)get_axis_value(line, 'S');
        }
    } else if (line.compare(0, 4, "M127") == 0 || line.compare(0, 4, "M107") == 0) {
        return 0;
    } else if ((flavor == (gcfMakerWare) || flavor == (gcfSailfish)) && line.compare(0, 4, "M126") == 0) {
        return (int16_t)get_axis_value(line, 'T');
    } else {
        return -1;
    }

}

void FanMover::_put_in_middle_G1(std::list<BufferData>::iterator item_to_split, float nb_sec_since_itemtosplit_start, BufferData &&line_to_write, float max_time) {
    assert(item_to_split != m_buffer.end());
    // if the fan is at the end of the g1 and the diff is less than 10% of the delay, then don't bother
    if (nb_sec_since_itemtosplit_start > item_to_split->time * 0.9 && (item_to_split->time - nb_sec_since_itemtosplit_start) < max_time * 0.1) {
        // doesn't really need to be split, print it after
        m_buffer.insert(next(item_to_split), line_to_write);
    } else 
        // does it need to be split?
        // if it's almost at the start of the g1, and the time "lost" is less than 10%
        if (nb_sec_since_itemtosplit_start < item_to_split->time * 0.1 && nb_sec_since_itemtosplit_start < max_time * 0.1 &&
        // and the previous isn't a fan value
        (item_to_split == m_buffer.begin() || std::prev(item_to_split)->fan_speed < 0)) {
        // doesn't really need to be split, print it before
        //will also print before if line_to_split.time == 0
        m_buffer.insert(item_to_split, line_to_write);
    } else if (item_to_split->raw.size() > 2
        && item_to_split->raw[0] == 'G' && item_to_split->raw[1] == '1' && item_to_split->raw[2] == ' ') {
        float percent = nb_sec_since_itemtosplit_start / item_to_split->time;
        BufferData before = *item_to_split;
        before.time *= percent;
        item_to_split->time *= (1-percent);
        if (item_to_split->dx != 0) {
            before.dx = item_to_split->dx * percent;
            item_to_split->x += before.dx;
            item_to_split->dx = item_to_split->dx * (1-percent);
            change_axis_value(before.raw, 'X', before.x + before.dx, 3);
        }
        if (item_to_split->dy != 0) {
            before.dy = item_to_split->dy * percent;
            item_to_split->y += before.dy;
            item_to_split->dy = item_to_split->dy * (1 - percent);
            change_axis_value(before.raw, 'Y', before.y + before.dy, 3);
        }
        if (item_to_split->dz != 0) {
            before.dz = item_to_split->dz * percent;
            item_to_split->z += before.dz;
            item_to_split->dz = item_to_split->dz * (1 - percent);
            change_axis_value(before.raw, 'Z', before.z + before.dz, 3);
        }
        if (item_to_split->de != 0) {
            if (relative_e) {
                before.de = item_to_split->de * percent;
                change_axis_value(before.raw, 'E', before.de, 5);
                item_to_split->de = item_to_split->de * (1 - percent);
                change_axis_value(item_to_split->raw, 'E', item_to_split->de, 5);
            } else {
                before.de = item_to_split->de * percent;
                item_to_split->e += before.de;
                item_to_split->de = item_to_split->de * (1 - percent);
                change_axis_value(before.raw, 'E', before.e + before.de, 5);
            }
        }
        //add before then line_to_write, then there is the modified data.
        m_buffer.insert(item_to_split, before);
        m_buffer.insert(item_to_split, line_to_write);

    } else {
        //not a G1, print it before
        m_buffer.insert(item_to_split, line_to_write);
    }
}

void FanMover::_print_in_middle_G1(BufferData& line_to_split, float nb_sec_from_item_start, const std::string &line_to_write) {
    if (nb_sec_from_item_start > line_to_split.time * 0.9 && line_to_split.time < this->nb_seconds_delay / 4) {
        // doesn't really need to be split, print it after
        m_process_output += line_to_split.raw + "\n";
        m_process_output += line_to_write + (line_to_write.back() == '\n'?"":"\n");
    } else if (nb_sec_from_item_start < line_to_split.time * 0.1 && line_to_split.time < this->nb_seconds_delay / 4) {
        // doesn't really need to be split, print it before
        //will also print before if line_to_split.time == 0
        m_process_output += line_to_write + (line_to_write.back() == '\n' ? "" : "\n");
        m_process_output += line_to_split.raw + "\n";
    }else if(line_to_split.raw.size() > 2
        && line_to_split.raw[0] == 'G' && line_to_split.raw[1] == '1' && line_to_split.raw[2] == ' ') {
        float percent = nb_sec_from_item_start / line_to_split.time;
        std::string before = line_to_split.raw;
        std::string& after = line_to_split.raw;
        if (line_to_split.dx != 0) {
            change_axis_value(before, 'X', line_to_split.x + line_to_split.dx * percent, 3);
        }
        if (line_to_split.dy != 0) {
            change_axis_value(before, 'Y', line_to_split.y + line_to_split.dy * percent, 3);
        }
        if (line_to_split.dz != 0) {
            change_axis_value(before, 'Z', line_to_split.z + line_to_split.dz * percent, 3);
        }
        if (line_to_split.de != 0) {
            if (relative_e) {
                change_axis_value(before, 'E', line_to_split.de * percent, 5);
                change_axis_value(after, 'E', line_to_split.de * (1 - percent), 5);
            } else {
                change_axis_value(before, 'E', line_to_split.e + line_to_split.de * percent, 5);
            }
        }
        m_process_output += before + "\n";
        m_process_output += line_to_write + (line_to_write.back() == '\n' ? "" : "\n");
        m_process_output += line_to_split.raw + "\n";

    } else {
        //not a G1, print it before
        m_process_output += line_to_write + (line_to_write.back() == '\n' ? "" : "\n");
        m_process_output += line_to_split.raw + "\n";
    }
}

void FanMover::_remove_slow_fan(int16_t min_speed, float past_sec) {
    //erase fan in the buffer -> don't slowdown if you are in the process of step-up.
    //we began at the "recent" side , and remove as long as we don't push past_sec to 0
    auto it = m_buffer.begin();
    while (it != m_buffer.end() && past_sec > 0) {
        past_sec -= it->time;
        if (it->fan_speed >= 0 && it->fan_speed < min_speed){
            //found something that is lower than us
            it = remove_from_buffer(it);

        } else {
            ++it;
        }
    }

}

std::string FanMover::_set_fan(int16_t speed) {
    assert(m_current_extruder < 200);
    assert(speed < 256 && speed >= 0);
    std::string str  = GCodeWriter::set_fan(m_writer.gcode_config(), m_current_extruder < 200 ? m_current_extruder : 0, uint8_t(speed));
    if(!str.empty() && str.back() == '\n')
        return str.substr(0,str.size()-1);
    return str;
}


bool parse_number(const std::string_view sv, int& out)
{
    {
        // Legacy conversion, which is costly due to having to make a copy of the string before conversion.
        try {
            assert(sv.size() < 1024);
            assert(sv.data() != nullptr);
            std::string str{ sv };
            size_t read = 0;
            out = std::stoi(str, &read);
            return str.size() == read;
        }
        catch (...) {
            return false;
        }
    }
}

//FIXME: add other firmware
// or just create that damn new gcode writer arch
void FanMover::_process_T(const std::string_view command)
{
    if (command.length() > 1) {
        int eid = 0;
        if (!parse_number(command.substr(1), eid) || eid < 0 || eid > 255) {
            GCodeFlavor flavor = m_writer.gcode_config().gcode_flavor;
            // Specific to the MMU2 V2 (see https://www.help.prusa3d.com/en/article/prusa-specific-g-codes_112173):
            if ((flavor == gcfMarlinLegacy || flavor == gcfMarlinFirmware) && (command == "Tx" || command == "Tc" || command == "T?"))
                return;

            // T-1 is a valid gcode line for RepRap Firmwares (used to deselects all tools) see https://github.com/prusa3d/PrusaSlicer/issues/5677
            if ((flavor != gcfRepRap && flavor != gcfSprinter) || eid != -1)
                m_current_extruder = static_cast<uint16_t>(0);
        } else {
            m_current_extruder = static_cast<uint16_t>(eid);
        }
    }
}


void FanMover::_process_ACTIVATE_EXTRUDER(const std::string_view cmd)
{
    if (size_t cmd_end = cmd.find("ACTIVATE_EXTRUDER"); cmd_end != std::string::npos) {
        size_t extruder_pos_start = cmd.find("EXTRUDER", cmd_end + std::string_view("ACTIVATE_EXTRUDER").size()) + std::string_view("EXTRUDER").size();
        assert(cmd[extruder_pos_start - 1] == 'R');
        if (extruder_pos_start != std::string::npos) {
            //remove next char until '-' or [0-9]
            while (extruder_pos_start < cmd.size() && (cmd[extruder_pos_start] == ' ' || cmd[extruder_pos_start] == '=' || cmd[extruder_pos_start] == '\t'))
                ++extruder_pos_start;
            size_t extruder_pos_end = extruder_pos_start + 1;
            while (extruder_pos_end < cmd.size() && cmd[extruder_pos_end] != ' ' && cmd[extruder_pos_end] != '\t' && cmd[extruder_pos_end] != '\r' && cmd[extruder_pos_end] != '\n')
                ++extruder_pos_end;
            std::string_view extruder_name = cmd.substr(extruder_pos_start, extruder_pos_end-extruder_pos_start);
            // we have a "name". It may be whatever or "extruder" + X
            for (const Extruder &extruder : m_writer.extruders()) {
                if (m_writer.gcode_config().tool_name.get_at(extruder.id()) == extruder_name) {
                    m_current_extruder = static_cast<uint16_t>(extruder.id());
                    return;
                }
            }
            std::string extruder_str("extruder");
            if (extruder_str == extruder_name) {
                m_current_extruder = static_cast<uint16_t>(0);
                return;
            }
            for (const Extruder &extruder : m_writer.extruders()) {
                if (extruder_str + std::to_string(extruder.id()) == extruder_name) {
                    m_current_extruder = static_cast<uint16_t>(extruder.id());
                    return;
                }
            }
        }
        BOOST_LOG_TRIVIAL(error) << "invalid ACTIVATE_EXTRUDER gcode command: '" << cmd << "', ignored by the fam mover post-process.";
    }
}

void FanMover::_process_gcode_line(GCodeReader& reader, const GCodeReader::GCodeLine& line)
{
    // processes 'normal' gcode lines
    bool need_flush = false;
    std::string cmd(line.cmd());
    double time = 0;
    int16_t fan_speed = -1;
    if (cmd.length() > 1) {
        if (line.has_f())
            m_current_speed = line.f() / 60.0f;
        switch (::toupper(cmd[0])) {
        case 'A':
            _process_ACTIVATE_EXTRUDER(line.raw());
                break;
        case 'T':
        case 't':
            _process_T(cmd);
                break;
        case 'G':
        {
            if (::atoi(&cmd[1]) == 1 || ::atoi(&cmd[1]) == 0) {
                double distx = line.dist_X(reader);
                double disty = line.dist_Y(reader);
                double distz = line.dist_Z(reader);
                double dist = distx * distx + disty * disty + distz * distz;
                if (dist > 0) {
                    dist = std::sqrt(dist);
                    time = dist / m_current_speed;
                }
            } else if (::atoi(&cmd[1]) == 2 || ::atoi(&cmd[1]) == 3) {
                // TODO: compute real dist
                double distx = line.dist_X(reader);
                double disty = line.dist_Y(reader);
                double dist = distx * distx + disty * disty;
                if (dist > 0) {
                    dist = std::sqrt(dist);
                    time = dist / m_current_speed;
                }
            }
            break;
        }
        case 'M':
        {
            fan_speed = get_fan_speed(line.raw(), m_writer.gcode_config().gcode_flavor);
            if (fan_speed >= 0) {
                const auto fan_baseline = (m_writer.gcode_config().fan_percentage.value ? 100.0 : 255.0);
                fan_speed = 100 * fan_speed / fan_baseline;
                if (!m_is_custom_gcode) {
                    // if slow down => put in the queue. if not =>
                    if (m_current_kickstart.time > 0) {
                        assert(m_back_buffer_fan_speed == m_current_kickstart.fan_speed);
                    }
                    if (m_back_buffer_fan_speed >= fan_speed) {
                        if (m_current_kickstart.time > 0) {
                            // stop kiskstart, and slow down
                            m_current_kickstart.time = -1;
                            //this fan speed will be printed, to make and end to the kickstart
                        }
                    } else {
                        if (nb_seconds_delay > 0 && (!only_overhangs || current_role == GCodeExtrusionRole::OverhangPerimeter)) {
                            //don't put this command in the queue
                            time = -1;
                            // this M106 need to go in the past
                            //check if we have ( kickstart and not in slowdown )
                            if (kickstart > 0 && fan_speed > m_front_buffer_fan_speed) {
                                //stop current kickstart , it's not relevant anymore
                                if (m_current_kickstart.time > 0) {
                                    m_current_kickstart.time = (-1);
                                }

                                //if kickstart
                                // first erase everything lower than that value
                                _remove_slow_fan(fan_speed, m_buffer_time_size + 1);
                                // then erase everything lower that kickstart
                                _remove_slow_fan(fan_baseline, kickstart);
                                // print me
                                if (!m_buffer.empty() && (m_buffer_time_size - m_buffer.front().time * 0.1) > nb_seconds_delay) {
                                    _print_in_middle_G1(m_buffer.front(), m_buffer_time_size - nb_seconds_delay, _set_fan(100));//m_writer.set_fan(100, true)); //FIXME extruder id (or use the gcode writer, but then you have to disable the multi-thread thing
                                    remove_from_buffer(m_buffer.begin());
                                } else {
                                    m_process_output += _set_fan(100) + "\n";//m_writer.set_fan(100, true)); //FIXME extruder id (or use the gcode writer, but then you have to disable the multi-thread thing
                                }
                                //write it in the queue if possible
                                const float kickstart_duration = kickstart * float(fan_speed - m_front_buffer_fan_speed) / 100.f;
                                float time_count = kickstart_duration;
                                auto it = m_buffer.begin();
                                while (it != m_buffer.end() && time_count > 0) {
                                    time_count -= it->time;
                                    if (time_count< 0) {
                                        //found something that is lower than us
                                        _put_in_middle_G1(it, it->time + time_count, BufferData(std::string(line.raw()), 0, fan_speed, true), nb_seconds_delay);
                                        //found, stop
                                        break;
                                    }
                                    ++it;
                                }
                                if (time_count > 0) {
                                    //can't place it in the buffer, use m_current_kickstart
                                    m_current_kickstart.fan_speed = fan_speed;
                                    m_current_kickstart.time = time_count;
                                    m_current_kickstart.raw = line.raw();
                                }
                                m_front_buffer_fan_speed = fan_speed;
                            } else {
                                // first erase everything lower than that value
                                _remove_slow_fan(fan_speed, m_buffer_time_size + 1);
                                // then write the fan command
                                if (!m_buffer.empty() && (m_buffer_time_size - m_buffer.front().time * 0.1) > nb_seconds_delay) {
                                    _print_in_middle_G1(m_buffer.front(), m_buffer_time_size - nb_seconds_delay, line.raw());
                                    remove_from_buffer(m_buffer.begin());
                                } else {
                                    m_process_output += line.raw() + "\n";
                                }
                                m_front_buffer_fan_speed = fan_speed;
                            }
                        } else {
                            if (kickstart <= 0) {
                                //nothing to do
                                //we don't put time = -1; so it will printed in the buffer as other line are done
                            } else if (m_current_kickstart.time > 0) {
                                //cherry-pick this one
                                if (m_back_buffer_fan_speed >= fan_speed) {
                                    //stop kickstart
                                    m_current_kickstart.time = -1;
                                    //this will print me just after as time >=0
                                } else {
                                    // add some duration to the kickstart and use it for me.
                                    float kickstart_duration = kickstart * float(fan_speed - m_back_buffer_fan_speed) / 100.f;
                                    m_current_kickstart.fan_speed = fan_speed;
                                    m_current_kickstart.time += kickstart_duration;
                                    m_current_kickstart.raw = line.raw();
                                    //i'm printed by the m_current_kickstart
                                    time = -1;
                                }
                            } else if(m_back_buffer_fan_speed < fan_speed - 10){ //only kickstart if more than 10% change
                                //don't write this line, as it will need to be delayed
                                time = -1;
                                //get the duration of kickstart
                                float kickstart_duration = kickstart * float(fan_speed - m_back_buffer_fan_speed) / 100.f;
                                //if kickstart, write the M106 S[fan_baseline] first
                                //set the target speed and set the kickstart flag
                                put_in_buffer(BufferData(_set_fan(100)//m_writer.set_fan(100, true)); //FIXME extruder id (or use the gcode writer, but then you have to disable the multi-thread thing
                                    , 0, fan_speed, true));
                                //kickstart!
                                //m_process_output += m_writer.set_fan(100, true) + "\n";
                                //add the normal speed line for the future
                                m_current_kickstart.fan_speed = fan_speed;
                                m_current_kickstart.time = kickstart_duration;
                                m_current_kickstart.raw = line.raw();
                            }
                        }
                    }
                    //update back buffer fan speed
                    m_back_buffer_fan_speed = fan_speed;
                } else {
                    // have to flush the buffer to avoid erasing a fan command.
                    need_flush = true;
                }
            }
            break;
        }
        }
    } else {
        if(!line.raw().empty() && line.raw().front() == ';')
        {
            if (line.raw().size() > 10 && line.raw().rfind(";TYPE:", 0) == 0) {
                // get the type of the next extrusions
                std::string extrusion_string = line.raw().substr(6, line.raw().size() - 6);
                current_role                 = string_to_gcode_extrusion_role(extrusion_string);
                assert(current_role != GCodeExtrusionRole::None);
            }
            if (line.raw().size() > 16) {
                if (line.raw().rfind("; custom gcode", 0) != std::string::npos) {
                    if (line.raw().rfind("; custom gcode end", 0) != std::string::npos) {
                        m_is_custom_gcode = false;
                    } else {
                        m_is_custom_gcode = true;
                    }
                }
            }
        }
    }

    if (time >= 0) {
        BufferData& new_data = put_in_buffer(BufferData(line.raw(), time, fan_speed));
        if (line.has(Axis::X)) {
            new_data.x = reader.x();
            new_data.dx = line.dist_X(reader);
        }
        if (line.has(Axis::Y)) {
            new_data.y = reader.y();
            new_data.dy = line.dist_Y(reader);
        }
        if (line.has(Axis::Z)) {
            new_data.z = reader.z();
            new_data.dz = line.dist_Z(reader);
        }
        if (line.has(Axis::E)) {
            new_data.e = reader.e();
            if (relative_e) {
                new_data.de = line.e();
                // GCode reader doesn't know it's relative extrusion, we have to do it ourself.
                //assert(new_data.e == 0);
                new_data.e = 0;
            } else
                new_data.de = line.dist_E(reader);
        }
        assert(new_data.dx == 0 || reader.x() == new_data.x);
        assert(new_data.dx == 0 || std::abs(reader.x() + new_data.dx - line.x()) < 0.00001f);
        assert(new_data.dy == 0 || reader.y() == new_data.y);
        assert(new_data.dy == 0 || std::abs(reader.y() + new_data.dy - line.y()) < 0.00001f);
        assert(new_data.de == 0 || (relative_e?0:reader.e()) == new_data.e);
        assert(new_data.de == 0 || std::abs((relative_e?0.f:reader.e()) + new_data.de - line.e()) < 0.00001f);
        //assert(new_data.de == 0 ||(relative_e?0.f:reader.e()) + new_data.de == line.e());

        if (m_current_kickstart.time > 0 && time > 0) {
            m_current_kickstart.time -= time;
            if (m_current_kickstart.time < 0) {
                //prev is possible because we just do a emplace_back.
                _put_in_middle_G1(prev(m_buffer.end()), time + m_current_kickstart.time, BufferData{ m_current_kickstart.raw, 0, m_current_kickstart.fan_speed, true }, kickstart);
            }
        }
    }/* else {
        BufferData& new_data = put_in_buffer(BufferData("; del? "+line.raw(), 0, fan_speed));
        if (line.has(Axis::X)) {
            new_data.x = reader.x();
            new_data.dx = line.dist_X(reader);
        }
        if (line.has(Axis::Y)) {
            new_data.y = reader.y();
            new_data.dy = line.dist_Y(reader);
        }
        if (line.has(Axis::Z)) {
            new_data.z = reader.z();
            new_data.dz = line.dist_Z(reader);
        }
        if (line.has(Axis::E)) {
            new_data.e = reader.e();
            if (relative_e)
                new_data.de = line.e();
            else
                new_data.de = line.dist_E(reader);
        }
    }*/
    // puts the line back into the gcode
    //if buffer too big, flush it.
    if (time >= 0) {
        // Add EPSILON to allow to have a buffer even with 0 m_buffer_time_size, so multiple consecutive M106 can be culled.
        while (!m_buffer.empty() && (need_flush || m_buffer_time_size - m_buffer.front().time > nb_seconds_delay + EPSILON) ){
            write_buffer_data();
        }
    }
#if _DEBUG
    double sum = 0;
    for (auto& data : m_buffer) sum += data.time;
    assert( std::abs(m_buffer_time_size - sum) < 0.01);
#endif
}

void FanMover::write_buffer_data()
{
    BufferData &frontdata = m_buffer.front();
    if (frontdata.fan_speed < 0 || frontdata.fan_speed != m_front_buffer_fan_speed || frontdata.is_kickstart) {
        if (frontdata.is_kickstart && frontdata.fan_speed < m_front_buffer_fan_speed) {
            // you have to slow down! not kickstart! rewrite the fan speed.
            m_process_output += _set_fan(frontdata.fan_speed) + "\n";
            m_front_buffer_fan_speed = frontdata.fan_speed;
        } else {
            m_process_output += frontdata.raw + "\n";
            if (frontdata.fan_speed >= 0 || frontdata.is_kickstart) {
                // note that this is the only place where the fan_speed is set and we print from the buffer, as if the
                // fan_speed >= 0 => time == 0 and as this flush all time == 0 lines from the back of the queue...
                m_front_buffer_fan_speed = frontdata.is_kickstart ? 100 : frontdata.fan_speed;
            }
        }
    }
    remove_from_buffer(m_buffer.begin());
}

} // namespace Slic3r

