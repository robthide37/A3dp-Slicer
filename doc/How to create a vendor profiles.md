# Vendors Bundles

A vendor bundle allows you to easily create many different profiles from a common base, via an inheritance graph. It's all contained in a single INI file.

In this document, you will learn to:
 * how to create a vendor bundle manually
   * write your vendor section
   * write your printer sections
   * write your print sections
   * write your filament sections
 * how to create a vendor bundle automatically (not available yet)
 * set up a repository so the slicer can get the updates, download and install them.
   * It's possible to do that simply by cloning a GitHub repository.
   * It's also possible to set up your own REST server.

# Creating a vendor bundle manually

The vendor bundle is a `.ini` file named after the vendor ID. It may link to icons, textures, and models. These are stored inside a directory with the same name as the vendor ID.

  ```
  VENDORID.ini
  VENDORID/
  ├── bed_model.stl
  ├── icon.png
  ├── texture.png
  └── ...
  ```

## ini format

The vendor bundle is a single file in 'ini' format. There are sections that begin with `[xxxx]` and each line can have a key-value pair separated by a '='. The file needs to be formatted with UTF-8 charset.
Example:
```
[section_id:section_value]
setting_key1 = setting_value1
setting_key2 = setting_value2
```
Each setting key in a section needs to be unique.

First start by copying a simple example, like in [VendorBundleExemple.ini](https://github.com/SuperSlicer-org/Profile-Template/blob/main/profiles/VENDORID.ini)

The name of your file will be the 'id' of the vendor bundle, and this id needs to be parsed by any filesystem, so only use ASCII letters (a-zA-Z), numbers, '-' and '_' (and no spaces).

## Sections
Here are the different sections of the vendor profile INI file. Each is explained in more detail afterwards:

* `[vendor]`: the first section, containing information about the vendor, to be displayed in the wizard, and to help with (auto-)updates.
* `[printer_model]`: allows grouping printer presets into a common section, with its thumbnail for the wizard. It's the entry point to install vendor presets. When installed, the slicer will remember each {vendor-printer_model-variant} key that is installed.
* `[printer]`: printer preset with its extruder(s)
* `[filament]`: FFF filament preset
* `[print]`: FFF print preset
* `[sla_material]`: SLA resin preset
* `[sla_print]`: SLA print preset

For the printer, filament, and print sections, there are two kinds of names:
1) `[print:my_print_preset_key]`: this is a preset that is usable in the slicer and can be selected in the combo box in the interface.
2) `[print:*0.2mmheight_print*]`: this is a preset that isn't visible in the slicer; it's only used as a library of settings to be applied at the end to the first kind.

The printer, filament, and print sections can inherit from another one (of the same type) via the `inherit` setting. It can take a single value or an array, separated by `;`. Example:
```
[print:*common_print*]
perimeters = 3
...

[print:*one_per_top*]
only_one_perimeter_top = 1

[print:final_print_preset]
inherit = *common_print*; *one_per_top*
fill_density = 15%
```

## vendor section

The first section is `[vendor]` and contains the general attributes of the vendor bundle.  
Here, you have to set:
 * `id` the ID of the vendor bundle; must be the same as the filename (without the `.ini`)
 * `name` the name that can be displayed to the user (short version). If `{technology}` is present inside, it will be replaced by the technology string (FFF or SLA).
 * `full_name` another name that can be displayed to the user (long version). Can contain `{technology}`.
 * `description` a description of this vendor bundle (who wrote it, for what kind of printers, etc.)
 * `config_version` the version of this bundle file. Each time you modify it, you should increase it.
 * `config_update_rest` the URL to access the REST API. If you're using GitHub, it can be simplified to `org/repository`. If left blank, this vendor profile is local and can't be upgraded automatically. You'll need to import the new version file with the GUI.
 * `technologies` a tag that can be `FFF`, `SLA`, or `FFF;SLA`. The printers inside the bundle must match one of these tags.
 * `slicer` the slicer ID this bundle was built for (almost the same as the name; you can get it in the "About" dialog in the Help menu). If the bundle is built for PrusaSlicer (or has no/empty value), it will be converted by an internal script when installed. (This script isn't perfect, and some incompatibilities may remain, often in the custom G-code sections.) If the value matches the current slicer ID, no conversion will be done—so ensure your vendor bundle is correct to avoid issues.
 * `slicer_version` the version of the slicer this bundle was built for. The format is the same as `config_version`. Only slicers with versions equal to or higher than this can upgrade to this bundle version (unless forced via the GUI).

### version
The version format is `MAJOR.MINOR.COUNTER.PATCH-BETATAG+METADATA`. (the tag for beta and metadata can be swapped)

Example of valid versions:
```
2.7
2.7.0
2.7.61.12345
2.7.1.1-alpha2
2.7.1.1-beta1
2.7.1.1+2024.01.23
2.7.1.1+2024.01.23-susi
2.7.1.1-susi+2024.01.23
```

regex: `[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+(((-[a-zA-Z0-9]+)(\+[a-zA-Z0-9]+)?)|((\+[a-zA-Z0-9]+)(-[a-zA-Z0-9]+)?))?`

But we advise sticking if possible to a simple 3-digit `1.0.0`.

For ordering (which verison is newer than the other), the tag are removed and a lexical comparison is done, with not-beta version greater than beta ones.
`1.0.1` == `1.0.1+PRUSA` > `1.0.1-rc1` > `1.0.0`

note: if the version of the configuration has a beta tag, it's possible that the slicer won't notify the user that a newer version is available.

## printer_model sections

The next sections will be `[printermodel:MODELID]`, with `MODELID` being the internal ID name for this model. In the wizard, a model is a printer with its icon and a list of possible extruders to choose from below.

The attributes are:
 * `name` Displayed name in the wizard for this printer model.
 * `variants` Iarray of strings separated by `;` to link a printer preset to a variant. These strings are displayed in the wizard in this order, with the first being the 'default' for this printer model. You need to have one `[printer]` section for each `printer_variant` defined later.
 * `technology` FFF or SLA
 * `family` name of the row in the wizard: all printer_models with the same family string are placed on the same line. These families are in the same order as they appear in this file.
 * `thumbnail` name of the PNG file in your vendor directory that contains the thumbnail to display in the wizard (e.g., `thumbnail = Voron_v4_250_thumb.png`)
 * `bed_model` 3D model file (an STL) in the icon directory, if you want to have a custom bed model (e.g., `bed_model = bed-v4-250.stl`).
 * `bed_texture` texture for the bed (can be a simple image if you don't have a 3D model) in the icon directory e.g., `bed_texture = bed-v4-250.png`).
 * `bed_with_grid` 0 or 1 if you want to have a grid displayed on top of the bed texture.
 * `default_materials` a list of material/filament IDs to be imported by default if this model is selected. You should take care with this to avoid the "no filament selected for the printers" error message.

### printer / filament / print sections

In other sections, you have two kinds of section names:
 * `[printer:*common*]` When the name has `*` on each side, it's not a real preset but a pattern, possibly missing some fields, to be applied to another preset.
 * `[printer:bestprinter]` Without any `*`, this is a real preset that can be imported and used.

## print / sla_print / printer / filament / sla_material

The section value is used as an ID to reference the printer preset for the inheritance graph.  
There are two kinds of presets: 
- The normal ones that can be displayed and used
- The hidden ones that have a `*` as the first and last character and are only used as intermediate presets in the inheritance tree.

Examples:  
`[printer:*common*]`  
`[printer:Printer project 42-id.2]`  
`[printer:my_printer_id]`

### Slicer settings

To define a setting in a print, filament, or printer preset, use its key (see tooltip in the slicer) and assign a value:

- Boolean: `0` for false, `1` for true  
- Integer: number  
- Float: number with `.` as decimal separator (`42.1234`)  
- Percentage: Float + `%` (`42.1234%`)  
- String: direct value after `=` (`start_gcode = #custom start\nprint_start EXTRUDER=200`) — special characters: `\n`, `\\`, `\r`. You can also enclose in quotes: `start_gcode = "#custom start\nprint_start EXTRUDER=200"`  
- Enum: key for the value (get from saved config), e.g., `fill_pattern = gyroid`  
- Point: FloatxFloat (`thumbnails = 64x64`)  
- Graph: vector separated by `:`. Format: `start:end:curve_type:points...`  
  Example: `overhangs_dynamic_speed = 0:3:1:0x0:50x25:100x100`

If the setting is an array, separate values with `,`:  
`thumbnails = 64x64,400x300`

To deactivate a setting, prefix with `!`:  
`first_layer_extrusion_width = !0.42`  
For arrays: `idle_temperature = !30,!30`

You must set the entire vector if you want to change one element — partial replacements are not supported.

### inheritance

The `inherits` setting allows this preset to start with values already defined. It’s an array using `;` as a separator.  
It is optional.

Examples:  
`inherits = *common*; *hot*`  
`inherits = my_printer_id; *slow*`  


The preset start with nothing defined.
The first preset in `inherits` can set the values it defines.
Same for the next ones in the `inherits` array, and it erase the current value if it's already defined.
Then the settings defined in this Preset ca erase anything already defined by the `inherits`.

Advice about inheritance:
* If you define it, put it first after the section name, to increase its visibility.
* There is two kind of presets: 
  1) The first kind defines every settings, by themselves or from an inherited preset. We call them 'full preset'
  2) The second kind only defines a subset of the settings, and so don't inherit any 'full preset'. We call them 'flavour preset'
* Any final preset (without `*` as first and last character in their id) need to be a 'full preset'.
* Hidden presets (with `*` as first and last character in their id) can be 'full preset' or 'flavour preset'.
* A final preset can inherit a final preset.
* A 'full preset' should defines all settings or inherit another 'full preset'. Any setting that isn't defined in a final Preset won't be able to be reverted to 'system value' in the slicer.
* If a 'full preset' inherit another 'full preset', it needs to be the first in the `inherits` array (anything before will be erased by it).
* It's advised to use the word 'common' in any hidden 'full preset', to be sure you only put them as the first inheritance. It's easy to have a `*common*` root that defines everything and then specialised 'commons' like `*common_voron2*` that inherit another common. And of course avoid naming any 'flavour preset' as 'common'.

## `[printer]` Section
The printer section has some special setting:
* `printer_model`: Any final Preset (without `*` as first and last character in their id) need to have this setting defined (by itself or by inheritance) to link it to a printer_model. The value of this setting need to be the same as the printer_model section value. Example: `printer_model = Voron_v4_250`
* `printer_variant`: Any final Preset (without `*` as first and last character in their id) need to have this setting defined (by itself or by inheritance) to link it to a printer_model variant. The value of this setting need to be the same as one of the `variants` values of the linked printer_model. Examples: `printer_variant = 0.4` `printer_variant = highflow 0.6`
* `default_print_profile`: the print section value that should be the default when you select this printer.
* `default_filament_profile`: the filament section value that should be the default when you select this printer.
* `default_sla_print_profile`: the print_sla section value that should be the default when you select this printer.
* `default_sla_material_profile`: the sla_material section value that should be the default when you select this printer.

## `[print]`  and `[sla_print]` Section
it's advised to set a `compatible_printers_condition` or a `compatible_printers` to let the user select this print option only when the desired printer(s) is selected. It avoids polluting the interface's print profile list with print profile from other vendors.

## `[filament]`  and `[sla_material]` Section
You may want to set `compatible_printers_condition` or a `compatible_printers` for the default material profile (like "basic PLA") to not pollute the material list with a "basic pla" for each vendor installed.
If you define a special material only defined in this vendor profile, you can let these field blank to allow any printer to use this material.

You can also use `compatible_prints` and `compatible_prints_condition` if you need to have special filament setting for a special kind of print profile (like a speedy print profile that need higher temperature & more cooling).

The `@` that can be at the left of the ID (`[filament:abs @lol]`) TODO

# Creating a vendor bundle automatically

TODO: not available yet

By using the `create_vendor_bundle.exe` tool, you can create a rough vendor bundle.
In a directory, create "project", "config", and "config_bundle" directories.
In each, put the .3mf projects, the config.ini or the config_bundle.ini you want to have in your vendor bundle.
Execute the tool, passing the root directory path as an argument. It will create a vendor bundle from all your configurations. Note that it's mostly a first step—you need to review the result and build on top of it.

You can also use the "integrate" function of the tool to include a new project/config/bundle into an existing vendor profile file.

# verifying the bundle

## principles

A real preset needs to have all fields defined, or it can lead to weird issues afterward. If it's using multiple inherits, the first one should be the one that defines all fields (like a *common*). The others afterwards should only define a subset (or they will override the previous one, making them useless). It's a good idea to have one *common* for each type that defines everything, and to inherit from it.

## template

If you don't want to set in stone all settings (because you don't really know them or want to get the new defaults when they change), then you need to create a `template.ini` file that contains your bundle. Then, using the `preset_bundle_tool.exe`, ask it to convert your file for its slicer version. It will create your bundle file and set all unset settings (in the commons) to the current default values.

## validation

You need to use the `preset_bundle_tool.exe` to validate your vendor bundle to be sure it doesn't have anything problematic (missing filaments, unused presets, broken inheritance, settings with wrong values, etc.)

# Github API

The simplest way is to fork the project `SuperSlicer-org/Profile-Template` on GitHub. And to follow the readme.

The commit messages will be available to users as the changelog, so be sure to write them clearly.  
Don't hesitate to squash your commits to have a cleaner history before creating a release (avoid having 7+ commits like "oops, revert" or "another try at this thing again" between two releases).


## create a release

* First validate your vendor bundle with the validator tool from the slicer version you have in your bundle. Fix any errors it reports.
* Check that you have increased your version number in your `.ini` file.
* Commit
* Look at your Git log to ensure the commit messages are clean and understandable.
* When it's all good, push to github
* create the release:
  * if you have the github action, go the **action** tab and check that it succeed to create the new tag.
  * go to the **Releases** page and create a new release.
  * **Automatic**: In the top "tag" field, selecting the latest tag
  * **Manual**: write your version tag: (like `1.0.0=2.7.63.0`) the version tag is composed of your `config_version` attribute and your `slicer_version` attribute, separated by a `=`. (Copy-paste them beforehand to be sure.) Agree to "create new tag on release".
  * You can write a dedicated release message or use the automatic Git-generated release message (which uses your commit messages). This message is not used (yet) in the slicer GUI.
  * Note: It's not mandatory to create a release, the new version is available as soon as the tag is created.
  
## how to setup a private server for a vendor bundle

TODO: a Python server that has the right API, to be shared on GitHub.

You can host a web server that provides the required components like the GitHub API.  
In GitHub, the base URL that is stored in `config_update_rest` is `https://api.github.com/orguser/repository`.  
If it's shortened to `orguser/repository`, it is still converted to that full version internally.  
To make the profile redirect to your server, you need to set the `config_update_rest` to your server API endpoint, like `http://SuperPrinter.com/profiles/api`

### /description

This endpoint returns a JSON object that contains the same key-values as the `[vendor]` section in your `.ini`.

Example:  
`{"vendor":{"id":"example", "name":"Example","config_version":"1.0.0", ... }}`

You can also return the contents of a `.ini` file.  
It's the entry point for the slicer when a user wants to install a bundle using only the URL.

### /tags

Note: The call should provide GET parameters (with `?`)
 * `per_page=100`
 * `page=1`  
But you can ignore them and return all tags. The pagination system is not implemented yet.

At this endpoint, you should return a JSON array of objects that contain:

 * `name`: the release tag of this version  
 * `zipball_url`: the URL to download a ZIP file that contains the repository files, with the same structure as GitHub (there is a root directory inside)  
 * `commit.sha`: the commit SHA or any unique ID for the release  
 * `commit.url`: the URL to get the release JSON (see the `/commit` endpoint)

If you don’t return all tags, you should at least return the most recent ones.

Example:
```
[
	{
		"name": "1.0.0=2.7.62.0",
		"zipball_url": "http://SuperPrinter.com/profiles/download_zip/1.0.0_2.7.62.0.zip",
		"commit": {
			"sha": "fc5caaa3d32765fb9e51c392a40107f67efcbd22",
			"url": "http://SuperPrinter.com/profiles/api/commit/fc5caaa3d32765fb9e51c392a40107f67efcbd22"
		}
	},
	{
		"name": "1.0.0-beta=2.7.62.0",
		"zipball_url": "http://SuperPrinter.com/profiles/download_zip/1.0.0-beta=2.7.62.0.zip",
		"commit": {
			"sha": "184c923a26f8e7123935b11d4f2bbf22681798ec",
			"url": "http://SuperPrinter.com/profiles/api/commit/184c923a26f8e7123935b11d4f2bbf22681798ec"
		}
	}
]
```

### /commit

Note: You can change this endpoint name, as it is not hardcoded — it's given by `commit.url` from `/tags`.

It only needs to have the field `commit.message`, which contains the changelog for this commit.  
This is only called if it's the first release of the vendor profile.

Example:
```
{
	"commit": {
		"message": "First version of Example vendor bundle."
	}
}
```

### /compare/oldtagid...newtagid

This endpoint allows the system to retrieve the changelog between two versions.

Example: `/compare/1.0.0-beta=2.7.62.0...1.0.0=2.7.62.0`

It returns a JSON object.

Example:
```
{
	"commits": [
		{
			"commit": {
				"message": "Fix the extruders retractions\nFix the extruders wipe"
			}
		},
		{
			"commit": {
				"message": "Release 1.0.0"
			}
		}
	]
}
```

The `commits` array contains a list of `commit` objects (oldest first), each containing a `message`.  
The changelog is created by concatenating all the messages in reverse order, with an extra `\n` between them.  
The message from the first tag is not included in the list—so it only includes commit messages **after the first tag and up to (and including) the last tag**.

To simplify things, you can just return a single commit object with the full changelog for the release in the message field.  
However, this approach can be misleading in the GUI, which groups releases by slicer version—so be cautious if you have a complicated release tree.
