% -*- mode: SLang; mode: fold -*-

static define _get_table_names (file, env) %{{{
{
   variable list;
   
   if (file[0] != '@')
     list = file;
   else
     {
	variable fp = fopen (strtrim_beg (file, "@"), "r");
	if (fp == NULL)
	  return NULL;

	list = fgetslines (fp);
	list = array_map (String_Type, &strtrim, list);
	() = fclose (fp);
     }

   env = getenv (env);
   if (env == NULL)
     return list;

   list = array_map (String_Type, &path_concat, env, list);
   return list;
}

%}}}

public define ic_get_table_names (file)
{
   return _get_table_names (file, "IC_TABLE_DIR");
}

public define ntb_get_table_names (file)
{
   return _get_table_names (file, "NTB_TABLE_DIR");
}

private variable _IC_Table_File_Name;
private variable _NTB_Table_File_Name;

#ifexists IC_Table_File
  _IC_Table_File_Name = IC_Table_File;
#else
  _IC_Table_File_Name = "ic_table.dat";
#endif

#ifexists NTB_Table_File
  _NTB_Table_File_Name = NTB_Table_File;
#else
  _NTB_Table_File_Name = "ntb_ee_table.dat";
#endif

private variable lib_file = "libnonthermal.so";
add_compiled_function (lib_file, "sync", "syn_table.dat");
add_compiled_function (lib_file, "invc", _IC_Table_File_Name);
add_compiled_function (lib_file, "ntbrem", _NTB_Table_File_Name);

import("nonthermal");

define make_sync_table () %{{{
{
   variable file = "";
   
   switch (_NARGS)
     {
      case 0:
	% use defaults
     }   
     {
      case 1:
	file = ();
     }
     {
	% default
	_pop_n (_NARGS);
	usage ("status = make_sync_table ([filename])");
	return;
     }
   
   _sync_make_table (file);
}

%}}}

define make_invc_table () %{{{
{
   variable t = 0.0, file = "";
   
   switch (_NARGS)
     {
      case 0:
	% use defaults
     }   
     {
      case 1:
	t = ();
     }
     {
      case 2:
	(t, file) = ();
     }
     {
	% default
	_pop_n (_NARGS);
	usage ("status = make_invc_table ([temperature [, filename])");
	return;
     }
   
   if (typeof(t) != Array_Type)
     t = [t];
   
   array_map (Int_Type, &_invc_make_table, t, file);
}

%}}}

private variable NTB_Process_Type = Assoc_Type[];
NTB_Process_Type ["ee"] = NTB_ee;
NTB_Process_Type ["ep"] = NTB_ep;

static define ntb_parse_process_type (type) %{{{
{
   if (typeof(type) == Int_Type)
     return type;
   
   !if (assoc_key_exists (NTB_Process_Type, type))
     {
        vmessage ("invalid process type: %S", type);
        return NULL;
     }   
   return NTB_Process_Type[type];
}

%}}}

define make_ntbrem_table () %{{{
{
   variable file = "", process = NTB_ee;
   
   switch (_NARGS)
     {
      case 0:
	% use defaults
     }   
     {
      case 1:
	file = ();
     }
     {
      case 2:
        (file, process) = ();
     }   
     {
	% default
	_pop_n (_NARGS);
	usage ("status = make_ntbrem_table ([filename [, process]])");
	return;
     }
   
   
   if (typeof(file) != Array_Type)
     file = [file];
   
   process = array_map (Int_Type, &ntb_parse_process_type, process);
   
   array_map (Int_Type, &_ntb_make_table, file, process);
}

%}}}


define ntb_set_process_weights (pw) %{{{
{
   _ntb_set_process_weights (pw);
}

%}}}

define particle_info_struct () %{{{
{
   % the order of these struct fields must match
   % the order used in nonthermal-module.c:pop_density_info()
   variable s = struct
     {
	particle_type,
	  index, curvature, cutoff, n_GeV,
	  kT, n_th
     };

   s.particle_type = "electron";

   return s;
}

%}}}

define particle_type (t);  % declare

define struct_args (s) %{{{
{
   variable t = @s;
   t.particle_type = particle_type(s.particle_type);

   variable num, args;
   num = _push_struct_field_values (t);
   args = __pop_args (num);

   __push_args (reverse(args));
}

%}}}

private variable Particle_Types = Assoc_Type[Integer_Type];
Particle_Types ["proton"] = PROTON;
Particle_Types ["electron"] = ELECTRON;

define particle_type (t) %{{{
{
   if (typeof(t) == String_Type)
     {
	if (0 != assoc_key_exists (Particle_Types, t))
	  return Particle_Types [t];
     }

   return t;
}

%}}}

define find_gamma_min (s) %{{{
{
   return _find_gamma_min (struct_args (s));
}

%}}}

define nontherm_density (s) %{{{
{
   return s.n_GeV * _nontherm_density (struct_args(s));
}

%}}}

define nontherm_energy_density (s) %{{{
{
   return s.n_GeV * _nontherm_energy_density (struct_args(s));
}

%}}}

define force_charge_conservation (se, sp) %{{{
{
   variable pnorm;

   pnorm = conserve_charge (struct_args(se), struct_args(sp));

   if (pnorm <= 0.0)
     return -1;

   sp.n_GeV = pnorm;
   return 0;
}

%}}}

define invc_set_dilution_factors (f) %{{{
{
   _invc_set_dilution_factors (f);
}

%}}}

provide ("nonthermal");
