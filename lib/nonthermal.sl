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

private define pdf_param_defaults (i, val, freeze, min, max) %{{{
{
   return (val[i], freeze[i], min[i], max[i]);
}

%}}}

private define add_pdf_fitfun () %{{{
{
   variable name, param_names=NULL, value, freeze, min, max;

   switch (_NARGS)
     {
      case 1:
        name = ();
     }
     {
      case 6:
        (name, param_names, value, freeze, min, max) = ();
     }
     {
        %default:
        _pop_n (_NARGS);
        usage ("add_pdf_fitfun (name [, param_names, value, freeze, min, max])");
        return;
     }

   variable t
     = ["private define $name (l,h,p) {return (\"$name\", p);}"$,
        "private define ${name}_contin (x,p) {return (\"$name\", p);}"$,
        "add_slang_function (\"$name\", [&$name, &${name}_contin] %s);"$
        ];

   t = strjoin (t, " ");

   if (param_names != NULL)
     {
        param_names = array_map (String_Type, &make_printable_string, param_names);
        param_names = strjoin (param_names, ",");
        t = sprintf (t, ", [$param_names]"$);
     }
   else t = sprintf (t, "");

   eval(t, "Global");

   if (param_names != NULL)
     {
        set_param_default_hook ("$name"$, &pdf_param_defaults,
                                value, freeze, min, max);
     }
}

%}}}

private define init_pdfs () %{{{
{
   add_pdf_fitfun ("default", ["index", "curvature", "cutoff [TeV]"],
                   [2.0, 0.0, 10.0], [0, 1, 0],
                   [1.0, -1.0, 1.0],
                   [3.0, 1.0, 1.e2]);

   add_pdf_fitfun ("etot", ["index"],
                   [2.0], [0], [1.0], [3.0]);

   add_pdf_fitfun ("mori");

   add_pdf_fitfun ("dermer", ["index", "curvature"],
                   [2.0, 0.0], [0, 1],
                   [1.0, -1.0],
                   [3.0, 1.0]);

   add_pdf_fitfun ("ke_cutoff", ["index", "curvature", "cutoff [TeV]"],
                   [2.0, 0.0, 10.0], [0, 1, 0],
                   [1.0, -1.0, 1.0],
                   [3.0, 1.0, 1.e2]);

   add_pdf_fitfun ("cbreak", ["index", "curvature", "cutoff [TeV]", "break [TeV]"],
                   [2.0, 0.0, 10.0, 1.0], [0, 1, 0, 1],
                   [1.0, -1.0, 1.0, 0.0],
                   [3.0, 1.0, 1.e2, 1.e2]);
}

%}}}

%!%+
%\function{add_pdf}
%\synopsis{Add a particle distribution function}
%\usage{add_pdf (String_Type libname, String_Type pdfname [, param_names[], param_values[], freeze[], min[], max[]])}
%\description
%  The \var{add_pdf} function can be used to load a particle
%  distribution function (PDF) from a shared library.  An example
%  implementation of such an object is provided in the
%  \var{examples} directory of the source code distribution.
%  
%  The first two arguments give the name of the shared library
%  and the name of the PDF.  If the PDF has any fit parameters,
%  the remaining array arguments should provide the parameter names,
%  default values, default freeze/thaw state and default allowed
%  minimum and maximum values.
%  
%\seealso{}
%!%-
define add_pdf () %{{{
{
   variable msg = "add_pdf (file, name [, param_names, value, freeze, min, max])";
   variable file, name, param_names=NULL, value, freeze, min, max;

   switch (_NARGS)
     {
      case 2:
        (file, name) = ();
     }
     {
      case 7:
        (file, name, param_names, value, freeze, min, max) = ();
     }
     {
        % default:
        _pop_n (_NARGS);
        usage(msg);
        return;
     }

   add_user_pdf_intrin (file, "$name"$, "");

   if (param_names != NULL)
     add_pdf_fitfun (name, param_names, value, freeze, min, max);
   else
     add_pdf_fitfun (name);
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

   init_pdfs();
}

%}}}

nonthermal_init ();

autoload ("_invc_make_table", "invc_make_table");
autoload ("_ntbrem_make_table", "ntbrem_make_table");
autoload ("_sync_make_table", "sync_make_table");

%!%+
%\function{make_sync_table}
%\synopsis{Generate a synchrotron lookup table}
%\usage{make_sync_table ([String_Type file])}
%\description
%  The \var{make_sync_table} function can be used to generate
%  a new lookup table for the angular integration that arises
%  in the computation of the \var{sync} model.
%  
%  The file name may be provided as a parameter.  If no
%  file name is provided, the default file name is used.
%  The output file is a FITS bintable.
%  
%  To customize details of the table computation, modify
%  the script \var{lib/sync_make_table.sl} which is
%  installed in \var{${prefix}/share/isis/nonthermal}.
%  
%\seealso{_sync_angular_integral}
%!%-
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

%!%+
%\function{make_invc_table}
%\synopsis{Generate an inverse Compton lookup table}
%\usage{make_invc_table ([String_Type file])}
%\description
%  The \var{make_invc_table} function can be used to generate
%  a new lookup table for the integral over incident photons that 
%  arises in the computation of the \var{invc} model.
%  
%  The file name may be provided as a parameter.  If no
%  file name is provided, the default file name is used.
%  The output file is a FITS bintable.
%  
%  To customize details of the table computation, modify
%  the script \var{lib/invc_make_table.sl} which is
%  installed in \var{${prefix}/share/isis/nonthermal}.
%  
%  A blackbody radiation field is used by default.  To use
%  a different radiation field, it is necessary to modify
%  the source code [replacing src/bbody.c and making any
%  other necessary changes]
%  
%\seealso{_invc_photon_integral}
%!%-
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

%!%+
%\function{make_ntbrem_table}
%\synopsis{Generate an ee Bremsstrahlung lookup table}
%\usage{make_ntbrem_table ([String_Type file])}
%\description
%  The \var{make_ntbrem_table} function can be used to generate
%  a new lookup table for the electron-electron bremsstrahlung 
%  differential cross-sections that arise in the computation
%  of the \var{ntbrem} model.
%  
%  The file name may be provided as a parameter.  If no
%  file name is provided, the default file name is used.
%  The output file is a FITS bintable.
%  
%  To customize details of the table computation, modify
%  the script \var{lib/ntbrem_make_table.sl} which is
%  installed in \var{${prefix}/share/isis/nonthermal}.
%  
%\seealso{_ee_haug1, _ee_haug1_lab}
%!%-
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

%!%+
%\function{ntb_set_process_weights}
%\synopsis{Set the relative weights of ee and ep bremsstrahlung}
%\usage{ntb_set_process_weights (wt_ee, wt_ep)}
%\description
%  The \var{ntb_set_process_weights} function can be used to 
%  set the weights used for electron-electron and electron-proton
%  bremsstrahlung.  The total contribution from both processes is
%#v+
%    S(E) = wt_ee * S_ee(E) + wt_ep * S_ep(E)
%#v-
%  The fluxes are scaled relative to the ambient number density 
%  of nuclei, \var{n}. The default weights are appropriate for a 
%  fully ionized gas with cosmic abundances so that most 
%  electrons come from hydrogen and helium.  In that case, 
%  the number fractions of hydrogen and helium are roughly
%#v+
%     X_H  = 1.0/1.1 = 0.9091,
%     X_HE = 0.1/1.1 = 0.0909.
%#v-
%  The electron-electron bremsstrahlung contribution,
%  proportional to the ambient electron density has weight
%#v+
%     wt_ee = X_H + 2*X_HE = 1.09091
%#v-
%  and the total electron-proton bremsstrahlung contribution,
%  proportional to \var{Z^2n}, has weight
%#v+
%     wt_ep = X_H + 4*X_HE = 1.27273
%#v-
%
%\seealso{}
%!%-
define ntb_set_process_weights (ee, ep) %{{{
{
   _ntb_set_process_weights (ee, ep);
}

%}}}

%!%+
%\function{particle_info_struct}
%\synopsis{Get a template structure for specifying the PDF}
%\usage{Struct_Type = particle_info_struct ()}
%\description
%  The \var{particle_info_struct} returns a template 
%  structure with the following fields:
%#v+
%   __Name__       __Definition__
%   pdf_name        name of the particle distribution function
%   params          parameter values
%   particle_type   electron=0, proton=1
%   n_GeV           density at 1 GeV [cm^-3 GeV^-1]
%   kT              temperature of thermal particles
%   n_th            density of thermal particles
%#v-
%
%\seealso{}
%!%-
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

%!%+
%\function{find_momentum_min}
%\synopsis{Find the momentum at which the thermal and nonthermal PDFs intersect}
%\usage{pc_min = find_momentum_min (Struct_Type p)}
%\description
%  The \var{find_momentum_min} function takes a structure
%  defining the thermal and nonthermal components of the
%  particle distribution function and computes the momentum
%  at which those two distributions contribute equally.
%  
%  If the thermal PDF is everywhere below the nonthermal PDF,
%  the momentum at the thermal peak is returned.
%
%\seealso{particle_info_struct}
%!%-
define find_momentum_min (s) %{{{
{
   return _find_momentum_min (struct_args (s));
}

%}}}

%!%+
%\function{nontherm_density}
%\synopsis{Compute the density of nonthermal particles}
%\usage{n = nontherm_density (Struct_Type p)}
%\description
%  The \var{nontherm_density} function takes a structure
%  defining the thermal and nonthermal components of the
%  particle distribution function (PDF) and computes the integral
%  over the nonthermal PDF.  The lower limit of this
%  integral is the momentum at which the thermal and nonthermal
%  PDFs intersect. If the thermal PDF is everywhere below the 
%  nonthermal PDF, the momentum at the thermal peak is used.
%
%\seealso{find_momentum_min}
%!%-
define nontherm_density (s) %{{{
{
   return s.n_GeV * _nontherm_density (struct_args(s));
}

%}}}

%!%+
%\function{nontherm_energy_density}
%\synopsis{Compute the energy density of nonthermal particles}
%\usage{n = nontherm_energy_density (Struct_Type p)}
%\description
%  The \var{nontherm_energy_density} function takes a structure
%  defining the thermal and nonthermal components of the
%  particle distribution function (PDF) and computes the 
%  energy density in nonthermal particles by integrating over 
%  the nonthermal PDF.  The lower limit of this
%  integral is the momentum at which the thermal and nonthermal
%  PDFs intersect. If the thermal PDF is everywhere below the 
%  nonthermal PDF, the momentum at the thermal peak is used.
%
%\seealso{find_momentum_min}
%!%-
define nontherm_energy_density (s) %{{{
{
   return s.n_GeV * _nontherm_energy_density (struct_args(s));
}

%}}}

%!%+
%\function{force_charge_conservation}
%\synopsis{Enforce charge conservation by adjusting the proton norm}
%\usage{force_charge_conservation (electrons, protons [, method])}
%\description
%  The \var{force_charge_conservation} function take structures
%  defining the electron and proton particle distribution
%  functions and adjusts the normalization of the nonthermal
%  proton distribution (protons.n_GeV) to enforce some definition 
%  of charge conservation.  The supported definitions are
%#v+
%   method=0 equal density of electrons and protons at some chosen 
%            injection kinetic energy
%   method=1 equal integrated density of nonthermal electrons and protons
%#v-
%  The default is \var{method=0}.
%
%\seealso{particle_info_struct}
%!%-
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
