/*
 * Copyright (C) 2013 KLab Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

//
// Links related to EXIF, TIFF 6.0 and MPO (Multi Picture Object) format
//
// https://www.exif.org/Exif2-2.PDF
// https://www.itu.int/itudoc/itu-t/com16/tiff-fx/docs/tiff6.pdf
//
// https://en.wikipedia.org/wiki/JPEG#JPEG_Multi-Picture_Format
// http://www.cmsoft.com.br/downloads/cmsoft-stereoscopic-picture-editor-converter/3d-picture-gallery/
// https://dmitrybrant.com/2011/02/08/the-fujifilm-mpo-3d-photo-format
//

#include <windows.h>

#include "exif.h"

#ifdef __cplusplus
extern "C" {
#endif


static int Verbose = 0;
// public funtions

int exifsetverbose() {
    Verbose = 1;
    return Verbose + 1;
}

/**
    * setVerbose()
    *
    * Verbose output on/off
    *
    * parameters
    *  [in] v : 1=on  0=off
    */
void exif_Taratata(int v)
{
    Verbose = v;
}

#ifdef __cplusplus
} /* extern "C" */
#endif
