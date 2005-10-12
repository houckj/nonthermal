% -*- mode: SLang; mode: fold -*-

import("nonthermal");

$1 = path_concat (_nonthermal_install_prefix, "share/isis/nonthermal");
prepend_to_isis_load_path ($1);

private variable Data_Path = path_concat ($1, "data");

private define _get_table_names (file, env) %{{{
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

private define push_array_values (a) %{{{
{
   foreach (a)
     {
        ();
     }
}

%}}}

private variable IC_Keys = ["GAMMIN", "GAMMAX", "EFNMIN", "EFNMAX"];

public define ic_read_table_hook (file) %{{{
{
   variable t = struct
     {
        gamma_range, efinal_range,
          gammas, efinals, y
     };

   variable x = fits_read_table (file);
   variable k = fits_read_key_struct (file, push_array_values(IC_Keys));

   t.gamma_range = [k.gammin, k.gammax];
   t.efinal_range = [k.efnmin, k.efnmax];
   t.gammas = x.gamma;
   t.efinals = x.efinal;
   t.y = x.y;

   return t;
}

%}}}

private variable NTB_Keys = ["EPHMIN", "EPHMAX", "EKNMIN", "EKNMAX"];

public define ntb_read_table_hook (file) %{{{
{
   variable t = struct
     {
        ephoton_range, ekinetic_range,
          ephotons, ekinetics, y
     };

   variable x = fits_read_table (file);
   variable k = fits_read_key_struct (file, push_array_values(NTB_Keys));

   t.ephoton_range = [k.ephmin, k.ephmax];
   t.ekinetic_range = [k.eknmin, k.eknmax];
   t.ephotons = x.ephoton;
   t.ekinetics = x.ekinetic;
   t.y = x.y;

   return t;
}

%}}}

public define ic_get_table_names (file) %{{{
{
   return _get_table_names (file, "IC_TABLE_DIR");
}

%}}}

public define ntb_get_table_names (file) %{{{
{
   return _get_table_names (file, "NTB_TABLE_DIR");
}

%}}}

variable _sync_table_file_name;
#ifexists SYN_Table_File
_sync_table_file_name = SYN_Table_File;
#else
_sync_table_file_name = "syn_table.fits";
#endif

define _sync_table_file () %{{{
{
   return path_concat (Data_Path, _sync_table_file_name);
}

%}}}

variable _invc_table_file_name;
#ifexists IC_Table_File
_invc_table_file_name = IC_Table_File;
#else
_invc_table_file_name = "ic_table.fits";
#endif

define _invc_table_file () %{{{
{
   return path_concat (Data_Path, _invc_table_file_name);
}

%}}}

variable _ntbrem_table_file_name;
#ifexists NTB_Table_File
_ntbrem_table_file_name = NTB_Table_File;
#else
_ntbrem_table_file_name = "ntb_ee_table.fits";
#endif

define _ntbrem_table_file () %{{{
{
   return path_concat (Data_Path, _ntbrem_table_file_name);
}

%}}}

private define nonthermal_init () %{{{
{
   variable lib_file = "libnonthermal.so";
   add_compiled_function (lib_file, "sync", _sync_table_file());
   add_compiled_function (lib_file, "invc", _invc_table_file());
   add_compiled_function (lib_file, "ntbrem", _ntbrem_table_file());
   add_compiled_function (lib_file, "pizero");
}

%}}}

nonthermal_init ();

require ("sync_make_table");

define make_sync_table () %{{{
{
   variable file = _sync_table_file_name;

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

   variable t = _sync_make_table ();

   fits_write_binary_table (file, "TABLE", t);
}

%}}}

private define ic_write_table_row (i, t, fptr) %{{{
{
   variable firstrow = i;
   variable firstelem = 1;
   variable status;

   status = _fits_write_col (fptr, 1, firstrow, firstelem, t.gammas[i-1]);
   status = _fits_write_col (fptr, 2, firstrow, firstelem, t.efinals[i-1]);
   status = _fits_write_col (fptr, 3, firstrow, firstelem, t.y[i-1]);
}

%}}}

private define ic_write_table (t, file) %{{{
{
   variable fp = fits_open_file (file, "c");

   variable status, naxis2, ttype, tform, tunit, extname;
   naxis2 = length(t.gammas);
   ttype = ["gamma", "efinal", "y"];
   tform = ["1D", "1PD", "1PD"];
   tunit = NULL;
   extname = NULL;
   status = _fits_create_binary_tbl (fp, naxis2, ttype, tform, tunit, extname);

   fits_update_key (fp, IC_Keys[0], t.gamma_range[0]);
   fits_update_key (fp, IC_Keys[1], t.gamma_range[1]);
   fits_update_key (fp, IC_Keys[2], t.efinal_range[0]);
   fits_update_key (fp, IC_Keys[3], t.efinal_range[1]);

   array_map (Void_Type, &ic_write_table_row, [1:naxis2], t, fp);

   fits_close_file (fp);
}

%}}}

private define ic_make_table (t, file) %{{{
{
   ic_write_table (_invc_make_table (t), file);
}

%}}}

define make_invc_table () %{{{
{
   variable file = _invc_table_file_name;
   variable t = 0.0;

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

   array_map (Void_Type, &ic_make_table, t, file);
}

%}}}

private variable NTB_Process_Type = Assoc_Type[];
NTB_Process_Type ["ee"] = NTB_ee;
NTB_Process_Type ["ep"] = NTB_ep;

private define ntb_parse_process_type (type) %{{{
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

private define ntb_write_table_row (i, t, fptr) %{{{
{
   variable firstrow = i;
   variable firstelem = 1;
   variable status;

   status = _fits_write_col (fptr, 1, firstrow, firstelem, t.ephotons[i-1]);
   status = _fits_write_col (fptr, 2, firstrow, firstelem, t.ekinetics[i-1]);
   status = _fits_write_col (fptr, 3, firstrow, firstelem, t.y[i-1]);
}

%}}}

private define ntb_write_table (t, file) %{{{
{
   variable fp = fits_open_file (file, "c");

   variable status, naxis2, ttype, tform, tunit, extname;
   naxis2 = length(t.ephotons);
   ttype = ["ephoton", "ekinetic", "y"];
   tform = ["1D", "1PD", "1PD"];
   tunit = NULL;
   extname = NULL;
   status = _fits_create_binary_tbl (fp, naxis2, ttype, tform, tunit, extname);

   fits_update_key (fp, NTB_Keys[0], t.ephoton_range[0]);
   fits_update_key (fp, NTB_Keys[1], t.ephoton_range[1]);
   fits_update_key (fp, NTB_Keys[2], t.ekinetic_range[0]);
   fits_update_key (fp, NTB_Keys[3], t.ekinetic_range[1]);

   array_map (Void_Type, &ntb_write_table_row, [1:naxis2], t, fp);

   fits_close_file (fp);
}

%}}}

private define ntb_make_table (process, file) %{{{
{
   variable t = _ntb_make_table (ntb_parse_process_type(process));
   ntb_write_table (t, file);
}

%}}}

define make_ntbrem_table () %{{{
{
   variable file = _ntbrem_table_file_name;
   variable process = NTB_ee;

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

   array_map (Void_Type, &ntb_make_table, process, file);
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

private define struct_args (s) %{{{
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

define find_momentum_min (s) %{{{
{
   return _find_momentum_min (struct_args (s));
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

define force_charge_conservation () %{{{
{
   variable pnorm, se, sp;
   variable method = 0;
   
   switch (_NARGS)
     {
      case 2:
        (se, sp) = ();
     }
     {        
      case 3:
        (se, sp, method) = ();
     }
     {
        % default:
        vmessage ("Usage:  force_charge_conservation (e, p [, method]);");
        vmessage ("        method = 0 means equal injection densities (ne_inj = np_inj)");
        vmessage ("        method = 1 means equal integrated nonthermal densities (ne_tot=np_tot)");
        return;
     }   

   pnorm = conserve_charge (struct_args(se), struct_args(sp), method);
   if (pnorm < 0)
     return pnorm;

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
