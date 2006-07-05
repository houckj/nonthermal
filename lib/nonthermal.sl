% -*- mode: SLang; mode: fold -*-

import("nonthermal");

$1 = path_concat (path_dirname (__FILE__), "nonthermal");
prepend_to_isis_load_path ($1);
private variable Data_Path = path_concat ($1, "data");

private define push_array_values (a) %{{{
{
   foreach (a)
     {
        ();
     }
}

%}}}

public define invc_table_init_hook (file) %{{{
{
   variable tf = fits_read_table (sprintf ("%s[TABLE]", file));

   variable keys_hdu = sprintf ("%s[TABLE]", file);
   variable key_names = ["gammin", "gammax", "efnmin", "efnmax",
                         "sigma0", "omega0", "xepsilon"];
   variable s = fits_read_key_struct (keys_hdu, push_array_values(key_names));
   variable keys = Double_Type[length(key_names)];
   _for (0, length(key_names)-1, 1)
     {
        variable k = ();
        keys[k] = get_struct_field (s, key_names[k]);
     }

   variable tx = fits_read_table (sprintf ("%s[XGRID]", file));
   variable ty = fits_read_table (sprintf ("%s[YGRID]", file));
   variable t = struct {xgrid,ygrid,f};
   t.xgrid = tx.xgrid;
   t.ygrid = ty.ygrid;
   t.f = tf.f;

   variable b = struct
     {
        x, y, f
     };
   b.x = t.xgrid;
   b.y = t.ygrid;
   b.f = t.f;
   variable p = bspline_open_intrin (NULL, b, 6, 6);

   return (p, keys);
}

%}}}

public define ntbrem_table_init_hook (file) %{{{
{
   variable t = fits_read_table (file);
   variable bdry_hdu = sprintf ("%s[BOUNDARY]", file);
   variable bdry = fits_read_table (bdry_hdu);
   variable key_names = ["ekinmin", "ekinmax", "ephmin", "ephmax", "sigma0", "xepsilon"];
   variable s = fits_read_key_struct (bdry_hdu, push_array_values(key_names));
   variable keys = Double_Type[length(key_names)];
   _for (0, length(key_names)-1, 1)
     {
        variable k = ();
        keys[k] = get_struct_field (s, key_names[k]);
     }
   variable b = struct
     {
        x, y, f
     };
   b.x = t.xgrid;
   b.y = t.ygrid;
   b.f = t.f;
   variable p = bspline_open_intrin (NULL, b, 3, 3);
   return (p, bdry.xleft, bdry.yleft, keys);
}

%}}}

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

private define pdf_default (l,h,p) {return ("default", p);}
private define pdf_default_contin (x,p) {return ("default", p);}

private define pdf_etot (l,h,p) {return ("etot", p);}
private define pdf_etot_contin (x,p) {return ("etot", p);}

private define pdf_mori (l,h,p) {return ("mori", p);}
private define pdf_mori_contin (x,p) {return ("mori", p);}

private define pdf_dermer (l,h,p) {return ("dermer", p);}
private define pdf_dermer_contin (x,p) {return ("dermer", p);}

private define pdf_ke_cutoff (l,h,p) {return ("ke_cutoff", p);}
private define pdf_ke_cutoff_contin (x,p) {return ("ke_cutoff", p);}

private define pdf_cbreak (l,h,p) {return ("cbreak", p);}
private define pdf_cbreak_contin (x,p) {return ("cbreak", p);}

private define pdf_param_defaults (i, val, freeze, min, max) %{{{
{
   return (val[i], freeze[i], min[i], max[i]);
}

%}}}

private define init_pdf_options () %{{{
{
   add_slang_function ("pdf_default", [&pdf_default, &pdf_default_contin],
                       ["index", "curvature", "cutoff [TeV]"]);
   set_param_default_hook ("pdf_default", &pdf_param_defaults,
                           [2.0, 0.0, 10.0], [0, 1, 0],
                           [1.0, -1.0, 1.0],
                           [3.0, 1.0, 1.e2]);

   add_slang_function ("pdf_etot", [&pdf_etot, &pdf_etot_contin], ["index"]);
   set_param_default_hook ("pdf_etot", &pdf_param_defaults,
                           [2.0], [0], [1.0], [3.0]);

   add_slang_function ("pdf_mori", [&pdf_mori, &pdf_mori_contin]);

   add_slang_function ("pdf_dermer", [&pdf_dermer, &pdf_dermer_contin],
                       ["index", "curvature"]);
   set_param_default_hook ("pdf_dermer", &pdf_param_defaults,
                           [2.0, 0.0], [0, 1],
                           [1.0, -1.0],
                           [3.0, 1.0]);

   add_slang_function ("pdf_ke_cutoff", [&pdf_ke_cutoff, &pdf_ke_cutoff_contin],
                       ["index", "curvature", "cutoff [TeV]"]);
   set_param_default_hook ("pdf_ke_cutoff", &pdf_param_defaults,
                           [2.0, 0.0, 10.0], [0, 1, 0],
                           [1.0, -1.0, 1.0],
                           [3.0, 1.0, 1.e2]);

   add_slang_function ("pdf_cbreak", [&pdf_cbreak, &pdf_cbreak_contin],
                       ["index", "curvature", "cutoff [TeV]", "break [TeV]"]);
   set_param_default_hook ("pdf_cbreak", &pdf_param_defaults,
                           [2.0, 0.0, 10.0, 1.0], [0, 1, 0, 1],
                           [1.0, -1.0, 1.0, 0.0],
                           [3.0, 1.0, 1.e2, 1.e2]);
}

%}}}

define add_pdf () %{{{
{
   variable msg = "add_pdf (file, name, param_names, value, freeze, min, max)";
   variable file, name, param_names, value, freeze, min, max;
   
   if (_NARGS != 7)
     {
        _pop_n (_NARGS);
        usage(msg);
        return;
     }

   (file, name, param_names, value, freeze, min, max) = ();

   add_user_pdf_intrin (file, "$name"$, "");

   variable t
     = [
        "private define pdf_$name (l,h,p) {return (\"$name\", p);}"$,
        "private define pdf_${name}_contin (x,p) {return (\"$name\", p);}"$,
        "add_slang_function (\"pdf_$name\", [&pdf_$name, &pdf_${name}_contin], [%s]);"$
        ];

   t = strjoin (t, " ");
   param_names = array_map (String_Type, &make_printable_string, param_names);
   t = sprintf (t, strjoin (param_names, ","));
   eval(t, "Global");

   set_param_default_hook ("pdf_$name"$, &pdf_param_defaults,
                           value, freeze, min, max);
}

%}}}

private define add_function () %{{{
{
   variable lib_file = "libnonthermal.so";
   variable name, arg;

   try
     {
        if (_NARGS == 1)
          {
             name = ();
             add_compiled_function (lib_file, name);
          }
        else
          {
             (name, arg) = ();
             add_compiled_function (lib_file, name, arg);
          }
     }
   catch AnyError:
     {
        vmessage ("*** %s model is not loaded", name);
        del_function (name);
     }
}

%}}}

private define nonthermal_init () %{{{
{
   add_function ("sync", _sync_table_file());
   add_function ("invc", _invc_table_file());
   add_function ("ntbrem", _ntbrem_table_file());
   add_function ("pizero");

   init_pdf_options();
}

%}}}

nonthermal_init ();

autoload ("_invc_make_table", "invc_make_table");
autoload ("_ntbrem_make_table", "ntbrem_make_table");
autoload ("_sync_make_table", "sync_make_table");

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
	usage ("make_sync_table ([filename])");
	return;
     }

   variable t = _sync_make_table ();

   fits_write_binary_table (file, "TABLE", t);
}

%}}}

define make_invc_table () %{{{
{
   variable file = _invc_table_file_name;

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
	usage ("make_invc_table ([filename])");
	return;
     }

   _invc_make_table (file);
}

%}}}

define make_ntbrem_table () %{{{
{
   variable file = _ntbrem_table_file_name;

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
	usage ("make_ntbrem_table ([filename])");
	return;
     }

   _ntbrem_make_table (file);
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
        pdf_name, params, particle_type, n_GeV, kT, n_th
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
