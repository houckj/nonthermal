%prepend_to_isis_load_path ("/nfs/vex/d0/i686/share/isis");
%prepend_to_isis_module_path ("/nfs/vex/d0/i686/lib/isis/modules");
require ("nonthermal");

prepend_to_isis_load_path ("/vex/d0/i686/share/slsh/local-packages");
require ("xfig");

nontherm_pdf ("mori");

fit_fun("pizero(1)");
set_par ("pizero(1).norm", 1.0);
set_par ("pizero(1).index", 2.75);  % ignored
set_par ("pizero(1).cutoff", 1.e10, 0, 0, 0);

%Pizero_Interpolate=1;

variable e_mori_gev =
  [0.01, 0.01585, 0.02512, 0.03981, 0.06310,         
    0.1, 0.1585,  0.2512,  0.3981,  0.6310,       
   1.e0, 1.585,   2.512,   3.981,   6.310,      
   1.e1, 1.585e1, 2.512e1, 3.981e1, 6.310e1,
   1.e2, 1.585e2, 2.512e2, 3.981e2, 6.310e2,
   1.e3, 1.585e3, 2.512e3, 3.981e3, 6.310e3,
   1.e4, 1.585e4, 2.512e4, 3.981e4, 6.310e4,
   1.e5];
variable y_mori_median =
  [9.960e-26, 2.019e-25, 3.540e-25, 4.726e-25, 5.223e-25,
   4.937e-25, 3.937e-25, 2.413e-25, 1.175e-25, 5.116e-26, 
   1.961e-26, 6.861e-27, 2.125e-27, 6.348e-28, 1.910e-28, 
   5.698e-29, 1.673e-29, 4.854e-30, 1.398e-30, 4.015e-31,
   1.153e-31, 3.324e-32, 9.626e-33, 2.805e-33, 8.224e-34,
   2.423e-34, 7.175e-35, 2.133e-35, 6.365e-36, 1.906e-36,
   5.727e-37, 1.726e-37, 5.221e-38, 1.585e-38, 4.816e-39,
   1.440e-39];

#iffalse
variable y_mori_min =
  [5.872e-26, 2.581e-25, 2.088e-25, 1.333e-25,
   6.908e-26, 3.206e-26, 1.301e-26, 4.765e-27,
   1.546e-27, 4.812e-28, 1.497e-28, 4.579e-29,
   1.361e-29, 3.952e-30, 1.133e-30, 3.238e-31,
   9.277e-32, 1.941e-34, 4.587e-37, 1.153e-39];
variable y_mori_max =
  [1.380e-25, 7.652e-25, 6.006e-25, 3.549e-25,
   1.630e-25, 6.684e-26, 2.447e-26, 8.314e-27,
   2.524e-27, 7.494e-28, 2.266e-28, 6.822e-29,
   2.016e-29, 5.852e-30, 1.678e-30, 4.799e-31,
   1.375e-31, 2.878e-34, 6.798e-37, 1.709e-39];
#endif
variable e = 10.0^(9 + [-2:6:0.1]);
%variable y = get_cfun(e);
variable yy = get_cfun (e_mori_gev * 1.e9);

variable Xsize = 6*2.54;
variable Ysize = 4*2.54;

variable xmin, xmax, ymin, ymax;
xmin = -2.5 + 9;
xmax =  5.5 + 9;
ymin = 0.8;
ymax = 2.0;

variable w = xfig_plot_new (Xsize, Ysize);
xfig_plot_set_line_thickness (w, 1);
xfig_plot_define_world (w, xmin, xmax, ymin, ymax);

#ifexists xfig_plot_set_x_tic_labels_font
xfig_plot_set_x_tic_labels_font (w, "Large", "\Large"R, NULL);
xfig_plot_set_y_tic_labels_font (w, "Large", "\Large"R, NULL);
#endif

xfig_plot_add_x_axis (w, 0, "\Large $\log_{10} \omega $ [eV]"R);
xfig_plot_add_y_axis (w, 0, "\Large $S_\pi / S_{\pi,{\rm Mori}}$"R);

xfig_plot_set_line_thickness (w, 3);
xfig_plot_set_line_color (w, "black");

%pointstyle(-5);
%set_line_width(3);
xfig_plot_set_point_size (w, 5);
xfig_plot_lines (w, 
                 9 + log10(e_mori_gev),
                 yy / y_mori_median);
#iffalse
xfig_plot_points (w,
                 9 + log10(e_mori_gev),
                 yy / y_mori_median);
#endif

xfig_plot_set_line_thickness (w, 1);
%pointstyle(-1);

xfig_plot_lines (w, 
                 9 + log10([1.e-4, 1.e6]), 
                 ones(2));

xfig_render_object (w, "mori97_table1.eps");

variable err =  1 - yy/y_mori_median;
writecol ("mori_compare.txt", e_mori_gev, yy, y_mori_median, err);
vmessage ("fractional errors;  mean/median/max = %g/%g/%g",
          mean(abs(err)), median(abs(err)), max(abs(err)));

