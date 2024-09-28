///|/ Copyright (c) Prusa Research 2022 - 2023 Tomáš Mészáros @tamasmeszaros
///|/ Copyright (c) 2023 Mimoja @Mimoja
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "SLAArchiveWriter.hpp"
#include "SLAArchiveFormatRegistry.hpp"

namespace Slic3r {

std::unique_ptr<SLAArchiveWriter>
SLAArchiveWriter::create(OutputFormat archtype, const SLAPrinterConfig &cfg)
{
    std::unique_ptr<SLAArchiveWriter> ret;
    auto factory = get_writer_factory(archtype);

    if (factory)
        ret = factory(cfg);

    return ret;
}

} // namespace Slic3r
