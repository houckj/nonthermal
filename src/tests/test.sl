require ("nonthermal");

define fixup_units (e_MeV, index)
{
   %  Sturner's spectra have a residual dependence on e- power-law index
   %  and were computed using MeV scaling instead of GeV scaling:
   %  His N(E) = (pc/ 1 MeV)^(-Gamma) exp (-E/Ecut) cm^-3 MeV^-1
   %  My  N(E) = (pc/ 1 GeV)^(-Gamma) exp ((1 GeV - E)/Ecut) cm^-3 GeV^-1
   %  so:
   variable per_gev = 10.0^(3.0*(index-1.0));
   variable MeV = 1.0e6;

   variable y = e_MeV^2 * (1.e-3 * get_cfun (e_MeV * MeV)) / per_gev;

   return (e_MeV, y);
}

define dlogy (y)
{
   return 1.0 - shift(y,1)/y;
}


define plot_comparison (index)
{
   vmessage ("index = %g", index);

   variable nn = int (index*10);

   variable dir;
   %dir = "/nfs/triassic/d0/gea/science/analysis/sn1006/chandra/00732/anal01/models/";
   dir = "sturner/";
   variable radsor_file, e_s, s_t;
   radsor_file = sprintf ("%s/radsor%d.dat", dir, nn);
   (e_s,s_t) = readcol (radsor_file, 1,2);

   variable gamsor_file, e_ic, ic_t, ntbrems_t;
   gamsor_file = sprintf ("%s/gamsor%d_cmbonly.dat", dir, nn);
   (e_ic,ic_t, ntbrems_t) = readcol (gamsor_file, 1,2, 3);

   xlog;
   ylog;

   xlabel ("E [MeV]");
   ylabel (latex2pg("E^2 Q [MeV /s /cm^2]"));
   title (sprintf ("index = %g", index));

   xrange (min([e_s, e_ic]), max([e_s, e_ic]));
   yrange (1.e-30, 1.e-10);

   fit_fun ("sync(1,ke_cutoff(1))");
   set_par ("sync(1).B_tot", 1.0);
   set_par ("ke_cutoff(1).index",index);
   set_par ("ke_cutoff(1).curvature", 0);
   set_par ("ke_cutoff(1).cutoff", 10.0);
   plot (fixup_units (e_s, index), 1);
   oplot(e_s, abs(s_t), red);

   fit_fun ("invc(1,ke_cutoff(1))");
   set_par ("invc(1).T_photon", 2.725);
   set_par ("ke_cutoff(1).index",index);
   set_par ("ke_cutoff(1).curvature", 0);
   set_par ("ke_cutoff(1).cutoff", 10.0);
   oplot (fixup_units (e_ic, index), 1);
   oplot(e_ic, abs(ic_t), red);

   fit_fun ("ntbrem(1,ke_cutoff(1))");
   set_par ("ntbrem(1).norm", 0.1);
   set_par ("ke_cutoff(1).index",index);
   set_par ("ke_cutoff(1).curvature", 0);
   set_par ("ke_cutoff(1).cutoff", 10.0);
   % Sturner used target proton density = 0.1 cm^-3
   % process weight is Z^2 * number_density_fraction
   variable y, y_ee, y_ep;
#iftrue
#iffalse   
   variable X_He = 0.1;
   variable X_H = 1.0 - X_He;
   variable x = [(X_H + 2*X_He), X_H, X_He];  % e, p, He
   variable Z = [1.0, 1.0, 2.0];              % e, p, He
   variable weights = x*Z^2;
   ntb_set_process_weights (weights[0], weights[1]+weights[2]);
#else
   ntb_set_process_weights (1.2, 0.0);
   (, y_ee) = fixup_units (e_ic, index);
   oplot (e_ic, y_ee, green);
   ntb_set_process_weights (0.0, 1.4);
   (, y_ep) = fixup_units (e_ic, index);
   oplot (e_ic, y_ep, blue);
   ntb_set_process_weights (1.2, 1.4);
   (, y) = fixup_units (e_ic, index);
   oplot (e_ic, y, 1);
   oplot(e_ic, abs(ntbrems_t), red);
   
#endif  
#else
   ntb_set_process_weights (1.0, 0.0);
   (, y) = fixup_units (e_ic, index);
   oplot (e_ic, y, 1);
   oplot(e_ic, abs(ntbrems_t), red);
#endif   

   plot_pause;
   
#iftrue
   variable i = where (abs(ntbrems_t) > 0);
   ylin;yrange (0.1,2);
   ylabel (latex2pg("Q_{jch}/Q_{ss}"));
   plot (e_ic[i], y[i]/abs(ntbrems_t[i]));
   oplot (e_ic[i], y_ee[i]/abs(ntbrems_t[i]), green);
   oplot (e_ic[i], y_ep[i]/abs(ntbrems_t[i]), blue);
   variable p = get_plot_info();
   oplot (10^[p.xmin, p.xmax], [1, 1], red);

   yrange (-1,1);
   ylabel (latex2pg ("\\Delta Q/Q"));
   plot (e_ic[i], dlogy(ntbrems_t[i]));
   oplot (e_ic[i], dlogy(y_ee[i]), green);
   oplot (e_ic[i], dlogy(y_ep[i]), blue);
   oplot (e_ic[i], dlogy(y_ep[i]+y_ee[i]), purple);   
   oplot ([0.511, 0.511], [-1,1], red);

   % What's going on near photon energy = electron rest energy?
   ntb_set_process_weights (1.2, 0.0);
   xrange(0.01,100);
   yrange; ylog;
   variable yy = e_ic^2 * get_cfun (e_ic * 1.e6);
   i = where (yy > 0); ylabel ("Q");
   plot(e_ic[i], yy[i]);
   yrange; ylin; ylabel (latex2pg ("\\Delta Q/Q"));
   plot(e_ic[i], dlogy(yy[i]));
   
   plot_pause;
   
#endif
}

% Compare with Steve Sturner's numbers

define main ()
{
   variable device, indices;
#iftrue
   device = "t.ps/cps";
   indices = [1.8, 2.0, 2.3];
#else
   device = "x.ps/cps";
   indices = [1.8, 2.0, 2.3];
   Ntb_Interpolate=0;   
#endif   
   
   variable id = plot_open (device);
   foreach (indices)
     {
        variable index = ();
        plot_comparison (index);
     }

   plot_close(id);
}

main ();

