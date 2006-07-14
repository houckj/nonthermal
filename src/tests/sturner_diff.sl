
variable IC_Table_File = "/nfs/cxc/h1/houck/src/modules/nonthermal/src/tests/ic_table_2.7K.fits";
vmessage ("Using IC table = %s", IC_Table_File);

%prepend_to_isis_load_path ("/tmp/i686/share/isis");
%prepend_to_isis_module_path ("/tmp/i686/lib/isis/modules");
require ("nonthermal");

prepend_to_isis_load_path ("/vex/d0/i686/share/slsh/local-packages");
require ("xfig");

define fixup_units (e_MeV, index)
{
   %  Sturner's spectra have a residual dependence on e- power-law index
   %  and were computed using MeV scaling instead of GeV scaling:
   %  His N(E) = (pc/ 1 MeV)^(-Gamma) exp (-E/Ecut) cm^-3 MeV^-1
   %  My  N(E) = (pc/ 1 GeV)^(-Gamma) exp ((1 GeV - E)/Ecut) cm^-3 GeV^-1
   %  so:
   variable per_gev = 10.0^(3.0*index);
   variable MeV = 1.0e6;

   %tic;
   variable y = e_MeV^2 * get_cfun (e_MeV * MeV) / per_gev;
   %vmessage ("time = %S", toc);

   return (e_MeV, y);
}

define dlogy (y)
{
   return 1.0 - shift(y,1)/y;
}

define sync_comparison (w, index)
{
   vmessage ("index = %g", index);

   variable nn = int (index*10);

   variable dir;
   %dir = "/nfs/triassic/d0/gea/science/analysis/sn1006/chandra/00732/anal01/models/";
   dir = "sturner/";
   variable radsor_file, e_s, s_t;
   radsor_file = sprintf ("%s/radsor%d.dat", dir, nn);
   (e_s,s_t) = readcol (radsor_file, 1,2);

   vmessage ("sync");
   fit_fun ("sync(1,ke_cutoff(1))");
   set_par ("sync(1).B_tot", 1.0);
   set_par ("ke_cutoff(1).index",index);
   set_par ("ke_cutoff(1).curvature", 0);
   set_par ("ke_cutoff(1).cutoff", 10.0);

   variable i;
   variable s_h;

   (, s_h) = fixup_units (e_s, index);
   i = where (s_t > 0);
   xfig_plot_lines (w,
                    log10(1.e6*e_s[i]),
                    s_h[i] / s_t[i]);
}

define gamma_comparison (w, index)
{
   vmessage ("index = %g", index);

   variable nn = int (index*10);

   variable dir;
   %dir = "/nfs/triassic/d0/gea/science/analysis/sn1006/chandra/00732/anal01/models/";
   dir = "sturner/";

   variable gamsor_file, e_ic, ic_t, ntbrems_t;
   gamsor_file = sprintf ("%s/gamsor%d_cmbonly.dat", dir, nn);
   (e_ic,ic_t, ntbrems_t) = readcol (gamsor_file, 1,2, 3);

   variable i;

   vmessage ("invc");
   fit_fun ("invc(1,ke_cutoff(1))");
   set_par ("invc(1).T_photon", 2.725);   % Sturner used 2.7 not 2.725
   set_par ("ke_cutoff(1).index",index);
   set_par ("ke_cutoff(1).curvature", 0);
   set_par ("ke_cutoff(1).cutoff", 10.0);

   variable ic_h;
   ( ,ic_h) = fixup_units (e_ic, index);
   i = where (ic_t > 0);
   xfig_plot_set_line_color (w, "black"); %"red");
   xfig_plot_lines (w,
                    log10(1.e6*e_ic[i]),
                    ic_h[i]/ic_t[i]);

   vmessage ("ntbrem");
   fit_fun ("ntbrem(1,ke_cutoff(1))");
   set_par ("ntbrem(1).norm", 0.11);
   set_par ("ke_cutoff(1).index",index);
   set_par ("ke_cutoff(1).curvature", 0);
   set_par ("ke_cutoff(1).cutoff", 10.0);

   variable ntb_h;
   %ntb_set_process_weights (1.2, 1.4);
   ntb_set_process_weights (1.090, 1.182);
   (, ntb_h) = fixup_units (e_ic, index);
   i = where (ntbrems_t > 0);

   xfig_plot_set_line_color (w, "black"); %"blue");
   xfig_plot_lines (w,
                    log10(1.e6*e_ic[i]),
                    ntb_h[i]/ntbrems_t[i]);

   xfig_plot_set_line_color (w, "black");
}

define make_y1_label (n)
{
   return sprintf ("{\Large $\bm{%5.3g}$}"R, n);
}

% Compare with Steve Sturner's numbers

define main ()
{
   variable indices = [1.8, 2.0, 2.3];
   variable line_styles = [1, 0, 2];

   variable Xsize = 6*2.54;
   variable Ysize = 4*2.54;

   variable xmin, xmax, ymin, ymax;
   xmin = -9.0;
   xmax =  5.0;
   ymin = 0.9;
   ymax = 1.03;

   variable w;
   w = xfig_plot_new (Xsize, Ysize);
   xfig_plot_set_line_thickness (w, 1);
   xfig_plot_define_world (w, xmin, xmax, ymin, ymax);
   
#ifexists xfig_plot_set_x_tic_labels_font
   xfig_plot_set_x_tic_labels_font (w, "Large", "\Large"R, NULL);
   xfig_plot_set_y_tic_labels_font (w, "Large", "\Large"R, NULL);
#endif
   
   xfig_plot_add_x_axis (w, 0, "\Large $\log_{10} \omega $ [eV]"R);
   xfig_plot_add_y_axis (w, 0, "\Large $S_{\rm sync} / S_{\rm Sturner}$"R);

#iftrue
   variable y1_major_tics, y1_minor_tics, y1_tic_labels;
   y1_major_tics = [ymin:ymax+0.02:0.02];
   y1_minor_tics = [ymin:ymax:0.01];
   y1_tic_labels = array_map (String_Type, &make_y1_label, y1_major_tics);
   xfig_plot_set_y1_tics (w, y1_major_tics, y1_tic_labels, y1_minor_tics);
#endif
   
   xfig_plot_set_line_thickness (w, 3);
   xfig_plot_set_line_color (w, "black");
   
   variable i;

   _for (0, 2, 1)
     {
        i = ();
        xfig_plot_set_line_style (w, line_styles[i]);
        sync_comparison (w, indices[i]);
     }
   
   xfig_plot_set_line_thickness (w, 1);
   xfig_plot_set_line_color (w, "black");
   xfig_plot_set_line_style (w, 0);
   xfig_plot_lines (w, [xmin,xmax], [1, 1]);
   
   variable legend, labels, colors, linestyles, thicknesses, width;
   labels = ["{\\large $\\Gamma$ = 2.3}", 
             "{\\large $\\Gamma$ = 2.0}",
             "{\\large $\\Gamma$ = 1.8}"];
   colors = "black";
   linestyles = [1, 0, 2];
   thicknesses = [3, 3, 3];
   width = 1;
   legend = xfig_new_legend (labels, colors, linestyles, thicknesses, width);
   xfig_plot_add_object (w, legend, -6.75, 0.92);

   xfig_render_object (w, "sturner_diff_sync.eps");

   xmin =  3.0;
   xmax = 14.5;
   ymin = 0.9;
   ymax = 1.1;

   w = xfig_plot_new (Xsize, Ysize);
   xfig_plot_set_line_thickness (w, 1);
   xfig_plot_define_world (w, xmin, xmax, ymin, ymax);
   
#ifexists xfig_plot_set_x_tic_labels_font
   xfig_plot_set_x_tic_labels_font (w, "Large", "\Large"R, NULL);
   xfig_plot_set_y_tic_labels_font (w, "Large", "\Large"R, NULL);
#endif
   
   xfig_plot_add_x_axis (w, 0, "\Large $\log_{10} \omega $ [eV]"R);
   xfig_plot_add_y_axis (w, 0, "\Large $S_{\rm invc,~ntbrem}/ S_{\rm Sturner}$"R);

#iftrue
   y1_major_tics = [ymin:ymax+0.05:0.05];
   y1_minor_tics = [ymin:ymax:0.01];
   y1_tic_labels = array_map (String_Type, &make_y1_label, y1_major_tics);
   xfig_plot_set_y1_tics (w, y1_major_tics, y1_tic_labels, y1_minor_tics);
#endif
   
   xfig_plot_set_line_thickness (w, 3);
   xfig_plot_set_line_color (w, "black");

   _for (0, 2, 1)
     {
        i = ();
        xfig_plot_set_line_style (w, line_styles[i]);
        gamma_comparison (w, indices[i]);
     }

   xfig_plot_set_line_thickness (w, 1);
   xfig_plot_set_line_color (w, "black");
   xfig_plot_set_line_style (w, 0);
   xfig_plot_lines (w, [xmin,xmax], [1, 1]);

   variable font = xfig_make_font ("\\bf", "\\Huge", 0x000000);
   xfig_plot_add_object (w, xfig_new_text ("IC", font), 7.5, 0.99, -0.5, 0);
   xfig_plot_add_object (w, xfig_new_text ("NB", font), 8.38, 1.026, -0.5, 0);

   labels = ["{\\large $\\Gamma$ = 2.3}", 
             "{\\large $\\Gamma$ = 2.0}",
             "{\\large $\\Gamma$ = 1.8}"];
   colors = "black";
   linestyles = [1, 0, 2];
   thicknesses = [3, 3, 3];
   width = 1;
   legend = xfig_new_legend (labels, colors, linestyles, thicknesses, width);
   xfig_plot_add_object (w, legend, 11.5, 1.07); %4.75, 1.07);
   
   xfig_render_object (w, "sturner_diff_gamma.eps");
}

main ();

