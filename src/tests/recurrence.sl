require ("nonthermal");
require ("histogram");
require ("xfig");

variable Energies, Model_Name, PDF_Name = "default";

define fcn (gamma, alpha)
{
   fit_fun (sprintf ("%s(1, %s(1))", Model_Name, PDF_Name));

   variable cutoff_tev = 1.0 / alpha;
   set_par (PDF + "(1).index", gamma, 0, 0, 0);
   set_par (PDF + "(1).cutoff", cutoff_tev, 0, 0, 0);

   return get_cfun (Energies);
}

define recurrence (gam, alpha, eps)
{
#iftrue
   % FIXME - short-circuit
   if (Model_Name != "pizero")
     return Energies * 0.0;
#endif

   variable sm, s0, sp, ss, dsda;

   % E0 is 1 GeV
   variable e0_tev = 1.e-3;

   variable am = alpha * (1.0 - 0.5*eps);
   variable a0 = alpha;
   variable ap = alpha * (1.0 + 0.5*eps);
   variable da = alpha * eps;
   
   sm = fcn (gam, am);
   s0 = fcn (gam, a0);
   sp = fcn (gam, ap);
   ss = fcn (gam-1.0, a0);
   dsda = (sp - sm) / da;

   variable err, i;
   err = Double_Type[length(s0)];
   i = where (s0 != 0);

#iftrue
   err[i] = ((1.0/e0_tev)/s0[i]) * dsda[i] + ss[i]/s0[i] - 1.0;
#iftrue
   % FIXME -- compute something more like a fractional error
   err[i] /=sqrt((((1.0/e0_tev)/s0[i]) * dsda[i])^2 + (ss[i]/s0[i])^2 + 1.0);
#endif

#else
   variable ssm, ssp;
   ssm = fcn (gam-1.0, am);
   ssp = fcn (gam-1.0, ap);
   err[i] = ((1.0/e0_tev)/s0[i]) * dsda[i]
     - 0.5 * ( sp[i] +  sm[i]) / s0[i]
     + 0.5 * (ssp[i] + ssm[i]) / s0[i];
#endif
   
#iffalse
   variable ii = where (abs(err[i]) > 1.e-3);
   if (length(ii) > 0)
     {
        vmessage ("gam = %15.9e  a0 = %15.9e", gam, a0);
        foreach (ii)
          {
             variable k = ();
             variable ik = i[k];

             () = fprintf (stdout, "%15.9e %17.12e %17.12e %17.12e %17.12e %15.9e\n",
                           Energies[ik],
                           sm[ik], s0[ik], sp[ik], ss[ik],
                           err[ik]);
%             () = fprintf (stdout, "%17.12e %17.12e %17.12e %17.12e %15.9e\n",
%                           Energies[ik], sp[ik]-sm[ik], sp[ik], sm[ik], err[ik]);
          }
        %plot_pause;
     }
#endif

#iffalse
   vmessage ("gam = %15.9e  am = %15.9e", gam, am);
   xlog;
   ylog;
   plot(Energies[i], s0[i]);
   plot_pause;
   oplot(Energies[i], sm[i]);
   plot_pause;
   oplot(Energies[i], sp[i]);
   plot_pause;
   oplot(Energies[i], ss[i]);
   plot_pause;
   set_float_format ("%16.9e");
   writecol (stdout, Energies[i], sm[i], s0[i], sp[i], ss[i]);
   plot_pause;
#endif

#iffalse
   xlog;
   ylog;
   %xrange (1.e11, 1.e14);
   plot (Energies[i], abs(err[i]));
   plot_pause;

   variable j = where (0 < Energies[i]);% and Energies[i] < 1.e14);
   variable t1 = ((1.0/e0_tev)/s0[i]) * dsda[i];
   variable t2 = ss[i]/s0[i];
   set_float_format ("%20.15e");
   writecol (stdout, Energies[i[j]], t1[j], t2[j], err[i[j]]);

#iffalse
   ylin;
   plot (Energies[i], 1.0 - abs(t2)/abs(t1));
   plot_pause;
   variable ee = Energies[i];
   variable d, s;
   d = t2 - shift(t2,1);
   s = t2 + shift(t2,1);
   plot (ee[[:-2]], (d/s)[[:-2]]);
   plot_pause;
   d = (dsda - shift(dsda,1))[i];
   s = (dsda + shift(dsda,1))[i];
   plot (ee[[:-2]], (d/s)[[:-2]]);
   plot_pause;
#endif
#endif

   return err;
}

define quick_main ()
{
   variable gam, alpha, eps, d;
   gam = 2.2;
   alpha = 0.001;    % 1/cutoff_tev
   eps = 2.5e-6;

   variable e_break_eV = 0.066666 * (1.0/alpha)^2;  % B_tot = 1 microgauss
   variable sync_energies = 10.0^[-5: log10(1.e5*e_break_eV) :0.01];

   variable gamma_energies;
#iftrue
   gamma_energies = 10.0^[5.0:log10(3.0e12/alpha):0.01];
#else
   gamma_energies = [
                     7.07945e+11,
                     7.0794578438e+11
                     %7.079457843841373e+11
                     ];
#endif

   variable d_sync, d_invc, d_ntbrem_ee, d_ntbrem_ep, d_pizero;

   Energies = sync_energies;
   Model_Name = "sync";
   d_sync = recurrence (gam, alpha, eps);

   Energies = gamma_energies;
   Model_Name = "invc";
   d_invc = recurrence (gam, alpha, eps);

   Model_Name = "ntbrem";
   ntb_set_process_weights (1.0, 0.0);
   d_ntbrem_ee = recurrence (gam, alpha, eps);

   ntb_set_process_weights (0.0, 1.0);
   d_ntbrem_ep = recurrence (gam, alpha, eps);

   Energies = gamma_energies;
   % pizero has some ugly glitches below ~200 MeV
   Model_Name = "pizero";
   d_pizero = recurrence (gam, alpha, eps);
}

define plot_main()
{
   variable gam, alpha, eps, d;
   gam = 2.2;
   alpha = 1.0/20;    % 1/cutoff_tev
   eps = 1.25e-6;

   variable sync_energies = 10.0^[-5:5:0.01];
   variable gamma_energies = 10.0^[5:14:0.01]; %[8.3:14:0.01];

   variable d_sync, d_invc, d_ntbrem_ee, d_ntbrem_ep, d_pizero;

   Energies = sync_energies;
   Model_Name = "sync";
   d_sync = recurrence (gam, alpha, eps);

   Energies = gamma_energies;
   Model_Name = "invc";
   d_invc = recurrence (gam, alpha, eps);

   Model_Name = "ntbrem";
   ntb_set_process_weights (1.0, 0.0);
   d_ntbrem_ee = recurrence (gam, alpha, eps);

   ntb_set_process_weights (0.0, 1.0);
   d_ntbrem_ep = recurrence (gam, alpha, eps);

   % pizero has some ugly glitches below ~200 MeV
   Model_Name = "pizero";
   d_pizero = recurrence (gam, alpha, eps);

   variable Xsize = 6*2.54;
   variable Ysize = 4*2.54;

   variable xmin, xmax, ymin, ymax;
   variable w;

   variable Scale = 1.e4;

   xmax = log10 (sync_energies[-1]);
   xmin = log10 (sync_energies[0]);
   ymax = max(d_sync*Scale);
   ymin = min(d_sync*Scale);

   w = xfig_plot_new (Xsize, Ysize);
   xfig_plot_set_line_thickness (w, 1);
   xfig_plot_define_world (w, xmin, xmax, ymin, ymax);
   xfig_plot_add_x_axis (w, 0, "$\log_{10} E_\gamma $ [eV]"R);
   xfig_plot_add_y_axis (w, 0, "$\epsilon \times 10^4$"R);

   xfig_plot_set_line_thickness (w, 2);

   xfig_plot_set_line_color (w, "black");
   xfig_plot_lines (w,
                    log10(sync_energies),
                    d_sync * Scale);

   xfig_render_object (w, "recur_sync.eps");

   xmax = log10 (gamma_energies[-1]);
   xmin = log10 (gamma_energies[0]);
   ymax = max([d_invc, d_ntbrem_ep, d_ntbrem_ee]*Scale);
   ymin = min([d_invc, d_ntbrem_ep, d_ntbrem_ee]*Scale);

   w = xfig_plot_new (Xsize, Ysize);
   xfig_plot_set_line_thickness (w, 1);
   xfig_plot_define_world (w, xmin, xmax, ymin, ymax);
   xfig_plot_add_x_axis (w, 0, "$\log_{10} E_\gamma $ [eV]"R);
   xfig_plot_add_y_axis (w, 0, "$\epsilon \times 10^4$"R);

   xfig_plot_set_line_thickness (w, 2);

   xfig_plot_inc_line_depth (w, 20);
   xfig_plot_set_line_color (w, "blue");
   xfig_plot_lines (w,
                    log10(gamma_energies),
                    d_invc*Scale);

   xfig_plot_inc_line_depth (w, -10);
   xfig_plot_set_line_color (w, "black");
   xfig_plot_lines (w,
                    log10(gamma_energies),
                    d_ntbrem_ep*Scale);

   xfig_plot_inc_line_depth (w, -10);
   xfig_plot_set_line_color (w, "red");
   xfig_plot_lines (w,
                    log10(gamma_energies),
                    d_ntbrem_ee*Scale);

   xfig_plot_set_line_color (w, "green");
   xfig_plot_lines (w,
                    log10(gamma_energies),
                    d_pizero*Scale);

   xfig_render_object (w, "recur_gamma.eps");

}

define converge_main()
{
   variable fp = fopen ("recurrence_converge.dat", "w");

   variable gam, alpha, eps, d;
   gam = 2.5;
   alpha = 1.0/20;    % 1/cutoff_tev

   variable sync_energies = 10.0^[-5:5:0.01];
   variable gamma_energies = 10.0^[5:14:0.01]; %[8.3:14:0.01];

   variable d_sync, d_invc, d_ntbrem_ee, d_ntbrem_ep, d_pizero;
   variable epsilons = reverse (10.0^[-7:-1:0.2]);

   foreach (epsilons)
     {
        eps = ();

        vmessage ("eps = %g", eps);

        Energies = sync_energies;
        Model_Name = "sync";
        d_sync = recurrence (gam, alpha, eps);

        Energies = gamma_energies;
        Model_Name = "invc";
        d_invc = recurrence (gam, alpha, eps);

        Model_Name = "ntbrem";
        ntb_set_process_weights (1.0, 0.0);
        d_ntbrem_ee = recurrence (gam, alpha, eps);
        ntb_set_process_weights (0.0, 1.0);
        d_ntbrem_ep = recurrence (gam, alpha, eps);

        Model_Name = "pizero";
        d_pizero = recurrence (gam, alpha, eps);

        () = fprintf (fp, "%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                      eps,
                      sqrt(mean(d_sync^2)),
                      sqrt(mean(d_invc^2)),
                      sqrt(mean(d_ntbrem_ee^2)),
                      sqrt(mean(d_ntbrem_ep^2)),
                      sqrt(mean(d_pizero^2))
                      );
        () = fflush(fp);
     }

   () = fclose (fp);
}

require ("maplib");

define map_main()
{
   variable map_file = "recur_map_3e_1000_frac_pizero.dat";
   variable hist_file = "recur_map_hist_3e_1000_frac_pizero.dat";

   variable fp = fopen (map_file, "w");

   variable eps = 2.5e-6;  % invc best:  1.25e-6:1.25e-5   2.5e-6 is good

   variable d_sync, d_invc, d_ntbrem_ee, d_ntbrem_ep, d_pizero;
   variable sync_energies, gamma_energies;

   variable gammas = [1.8:4.0:0.125];
   variable cutoffs = 10.0^[1.0: 3.0 :0.125];  % also used 2.75
   variable alphas = reverse (1.0/cutoffs);

   variable gg, aa;
   (gg, aa) = maplib_meshgrid (gammas, alphas);

   variable n = length (gg);
   vmessage ("n = %d", n);

   reshape(gg,n);
   reshape(aa,n);

   variable rr_grid = 10.0^[-5.0:5.0:0.0125];
   rr_grid = [0.0, rr_grid];
#iftrue
   % adjust for "frac error" recurrence  [FIXME]
   rr_grid /= 1.e5;
#endif
   variable Nrr = Array_Type[5];
   Nrr[*] = Double_Type[length(rr_grid)];

   _for (0, n-1, 1)
     {
        variable i = ();

        variable gamma, alpha;
        gamma = gg[i];
        alpha = aa[i];

        % only consider gamma-ray spectra below N*E_cutoff.
        gamma_energies = 10.0^[5: log10(3*1.e12/alpha) :0.0125];

        % only consider sync spectra below N * E_break
        % in this case, N is really huge because the sync
        % spectrum is so broad -- maybe N = 1.e5 or more;
        variable e_break_eV = 0.066666 * (1.0/alpha)^2;  % B_tot = 1 microgauss
        sync_energies = 10.0^[-5: log10(1.e5*e_break_eV) :0.0125];

        Energies = sync_energies;
        Model_Name = "sync";
        d_sync = recurrence (gamma, alpha, eps);

        Energies = gamma_energies;
        Model_Name = "invc";
        d_invc = recurrence (gamma, alpha, eps);

        Model_Name = "ntbrem";
        ntb_set_process_weights (1.0, 0.0);
        d_ntbrem_ee = recurrence (gamma, alpha, eps);
        ntb_set_process_weights (0.0, 1.0);
        d_ntbrem_ep = recurrence (gamma, alpha, eps);

        Model_Name = "pizero";
        d_pizero = recurrence (gamma, alpha, eps);

        Nrr[0] += hist1d (abs(d_sync), rr_grid);
        Nrr[1] += hist1d (abs(d_invc), rr_grid);
        Nrr[2] += hist1d (abs(d_ntbrem_ee), rr_grid);
        Nrr[3] += hist1d (abs(d_ntbrem_ep), rr_grid);
        Nrr[4] += hist1d (abs(d_pizero), rr_grid);

        variable rms_sync, rms_invc, rms_ee, rms_ep, rms_pi;
        rms_sync = sqrt(mean(d_sync^2));
        rms_invc = sqrt(mean(d_invc^2));
        rms_ee = sqrt(mean(d_ntbrem_ee^2));
        rms_ep = sqrt(mean(d_ntbrem_ep^2));
        rms_pi = sqrt(mean(d_pizero^2));

        () = fprintf (fp, "%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                      gamma, alpha,
                      rms_sync, rms_invc, rms_ee, rms_ep, rms_pi);
        () = fflush(fp);

        % update in each loop so I can monitor progress while
        % it runs
        set_float_format ("%12.5e");
        writecol (hist_file,
                  rr_grid, Nrr[0], Nrr[1], Nrr[2], Nrr[3], Nrr[4]);
     }

   () = fclose (fp);

}

%quick_main ();
%converge_main();
%plot_main();
map_main ();
