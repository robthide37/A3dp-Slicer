//////////////////////////////////////////
// Api for SuperSlicer scripted widgets:
//
//// Functions callable ////
//
// -- to print on the console, for debugging --
// void print(string &out)
// void print_float(float)
//
// -- to get the value of real settings --
//  bool  get_bool(string &in key)
//  int   get_int(string &in key)
//    can be used by type int and enum (return the index)
//  float get_float(string &in key)
//    can be used by type float, percent and flaot_or_percent
//  bool  is_percent(string &in key)
//  void  get_string(string &in key, string &out get_val)
//    can be used by type string and enum (return the enum_value, not the label)
//
// -- to set the value of real settings --
//  void set_bool(string &in, bool new_val)
//  void set_int(string &in, int new_val)
//    if an enum, it's the index
//  void set_float(string &in, float new_val)
//    if a float_or_percent, unset the percent flag at the same time
//  void set_percent(string &in, float new_val)
//    if a float_or_percent, set the percent flag at the same time
//  void set_string(string &in, string &in new_val))
//    if an enum, it's one of the enum_value
//
//// Functions to define for each script widget //// 
//
// type bool:
//   int OPTNAME_get() 
//      will return 1 if checkd, 0 if unchecked and -1 if half-checked (not all os, will be uncehcked if not available)
//   void OPTNAME_set(bool set)
//
// type int:
//   int OPTNAME_get()
//   void OPTNAME_set(int set)
//
// type float & percent:
//   float OPTNAME_get()
//   void OPTNAME_set(float set)
//
// type float_or_percent:
//   float OPTNAME_get(bool &out is_percent)
//   void OPTNAME_set(float set, bool is_percent)
//
// type string:
//   void OPTNAME_get(string &out get)
//   void OPTNAME_set(string &in set)
//
// type enum:
//   int OPTNAME_get(string &out enum_value)
//      Only the return value is used unless it's out of bounds, then it tries to use the enum_value
//   void OPTNAME_set(string &in set_enum_value, int set_idx)
//
//

//overhangs : quick set/unset like the one in prusalicer

int s_overhangs_get()
{
	if (get_float("overhangs_width_speed") == 0) return 0;
	float width = get_float("overhangs_width");
	bool percent = is_percent("overhangs_width");
	if((percent && width > 50.f) || ((!percent) && width > 0.2f)) return 1;
	return -1;
}

void s_overhangs_set(bool set)
{
	if (set) {
		set_percent("overhangs_width_speed", 55.f);
		float width = get_float("overhangs_width");
		bool percent = is_percent("overhangs_width");
		if((percent && width < 50.f) || ((!percent) && width < 0.2f))
			set_percent("overhangs_width", 75.f);
	} else {
		set_float("overhangs_width_speed", 0.);
	}
}

void s_overhangs_reset(bool set)
{
	back_initial_value("overhangs_width_speed");
	back_initial_value("overhangs_width");
}

// "not thick bridge" like in prusaslicer

float compute_overlap()
{
	float height = get_float("layer_height");
	print("layer_height = " + height + "\n");
	float width = get_computed_float("solid_infill_extrusion_width");
	print("width = " + width + "\n");
	if(height <= 0) return 1;
	if(width <= 0) return 1;
	float solid_spacing = (width - height * 0.215);
	print("solid_spacing = " + solid_spacing + "\n");
	float solid_flow = height * solid_spacing;
	print("solid_flow = " + solid_flow + "\n");
	float bridge_spacing = sqrt(solid_flow*1.2739);
	print("bridge_spacing = " + bridge_spacing + "\n");
	print("bridge_spacing/solid_spacing = " + (bridge_spacing / solid_spacing) + "\n");
	return bridge_spacing / solid_spacing;
}

int s_not_thick_bridge_get()
{
	bool is_set = false;
	get_custom_bool(0,"not_thick_bridge", is_set);
	if(is_set){
		//set other vars
		ask_for_refresh();
		return 1;
	}
	return 0;
}

void s_not_thick_bridge_reset(bool set)
{
	set_custom_bool(0,"not_thick_bridge", false);
	back_initial_value("bridge_type");
	back_initial_value("bridge_overlap");
}

void s_not_thick_bridge_set(bool set)
{
	bool var_set = false;
	get_custom_bool(0,"not_thick_bridge", var_set);
	print("Me with value " + var_set +" has to be set to " + set+"\n");
	if (var_set != set) {
		set_custom_bool(0,"not_thick_bridge", set);
	}
	if (set) {
		if (get_int("bridge_type") != 2)
			set_int("bridge_type", 2);
		float overlap = compute_overlap();
		print("overlap = " + overlap + "\n");
		set_float("bridge_overlap", overlap);
		set_float("bridge_overlap_min", overlap);
	} else if (var_set != set) {
		back_initial_value("bridge_type");
		back_initial_value("bridge_overlap");
		back_initial_value("bridge_overlap_min");
	}
}


//test:
//	setting:script:bool:easy:depends$enforce_full_fill_volume:label$fullfill-lol:s_fullfill
//	setting:script:int:easy:depends$perimeters:label$perimeters-lol:s_perimeter
//	setting:script:float:easy:depends$top_solid_min_thickness:label$thickness-lol:s_thickness
//	setting:script:percent:easy:depends$bridge_flow_ratio:label$bridgeflow-lol:s_bridgeflow
//	setting:script:string:easy:depends$notes:label$notes-lol:s_notes
//	setting:script:enum$b$bof$m$mouaif:easy:depends$no_perimeter_unsupported_algo:label$noperi-lol:s_noperi

int s_fullfill_get()
{
	if (get_bool("enforce_full_fill_volume")) return 1;
	return 0;
}
void s_fullfill_set(bool set)
{
	set_bool("enforce_full_fill_volume", set);
}


int s_perimeter_get()
{
	return get_int("perimeters");
}
void s_perimeter_set(int set)
{
	set_int("perimeters", set);
}


float s_thickness_get()
{
	return get_float("top_solid_min_thickness");
}
void s_thickness_set(float set)
{
	set_float("top_solid_min_thickness", set);
}


float s_bridgeflow_get()
{
	return get_float("bridge_flow_ratio");
}
void s_bridgeflow_set(float set)
{
	set_percent("bridge_flow_ratio", set);
}


void s_notes_get(string &out get_val)
{
	get_string("notes", get_val);
}
void s_notes_set(string &out set_val)
{
	set_string("notes", set_val);
}


int s_noperi_get(string &out get_val)
{
	return get_int("no_perimeter_unsupported_algo") == 0 ? 0 : 1;
}
void s_noperi_set(string &out set_val, int idx)
{
	//set_int("no_perimeter_unsupported_algo", idx == 0 ? 0 : 3);
	if (idx == 0) set_int("no_perimeter_unsupported_algo",0);
	else set_string("no_perimeter_unsupported_algo", "filled");
}
