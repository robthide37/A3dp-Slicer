# What is `build_pot` for?

## Overview

`build_pot` is a GitHub Actions workflow job in the A3dp-Slicer repository that is responsible for **generating translation template files for internationalization (i18n)**.

## Purpose

The `build_pot` job extracts all translatable strings from the source code and creates POT (Portable Object Template) files that translators can use to create localized versions of the A3dp-Slicer application in different languages.

## Technical Details

### What are POT files?
- **POT** stands for "Portable Object Template"
- These are template files used in the GNU gettext internationalization system
- They contain all the translatable strings extracted from source code
- Translators use POT files as templates to create PO (Portable Object) files for specific languages

### How it works

1. **String Extraction**: The `build_pot` job runs the `gettext_make_pot` CMake target which:
   - Uses `xgettext` to scan source code for translatable strings
   - Looks for specific keywords like `L`, `_L`, `_u8L`, `L_CONTEXT`, `_L_PLURAL`
   - Extracts strings from files listed in `resources/localization/list.txt`
   - Generates `resources/localization/Slic3r.pot`

2. **Build Process**: The job includes:
   - Setting up Windows build environment
   - Installing gettext tools
   - Running `msbuild gettext_make_pot.vcxproj`

### Workflow Files

The `build_pot` job can be found in these GitHub Actions workflows:
- `.github/workflows/ccpp_win_debug.yml` (lines 31-74)
- `.github/workflows/ccpp_win.yml` (step "make .pot")  
- `.github/workflows/ccpp_win_rc.yml` (step "make .pot")

### CMake Targets

Related CMake targets in `CMakeLists.txt`:
- `gettext_make_pot`: Generates the POT file from source code strings
- `gettext_merge_community_po_with_pot`: Merges community translations with the new POT file

## Why is it separate from the main build?

The `build_pot` job runs independently because:
1. **Lightweight**: Only extracts strings, doesn't compile the full application
2. **Independent workflow**: Translation updates don't require rebuilding the entire app
3. **Translator workflow**: Allows updating translation templates without affecting the main build
4. **Automation**: Can be triggered separately when translation strings change

## Supported Languages

The repository currently supports localization for multiple languages, as seen in `resources/localization/`:
- German (de)
- Spanish (es) 
- French (fr)
- Hungarian (hu)
- Italian (it)
- Japanese (ja)
- Korean (ko)
- Dutch (nl)
- Polish (pl)
- Portuguese Brazil (pt_BR)
- Russian (ru)
- Turkish (tr)
- Ukrainian (uk)
- Chinese Simplified (zh_CN)
- Chinese Traditional (zh_TW)
- And more...

## Summary

**`build_pot` is for generating translation template files (POT files) that enable the A3dp-Slicer application to be translated into multiple languages for international users.**