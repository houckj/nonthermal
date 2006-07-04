require ("nonthermal");
require ("gsl");

Syn_Interpolate = 1;

define sync_numeric (nu, s, B)
{
   fit_fun ("sync(1,pdf_etot(1))");
   set_par ("sync(1).norm", 1);
   set_par ("sync(1).B_tot", B, 0, 0, 0);
   set_par ("pdf_etot(1).index", s, 0, 0, 0);

   variable energy_erg = CONST_CGSM_PLANCKS_CONSTANT_H * nu;

   % ph /s /cm^2 /GeV
   variable f = get_cfun (energy_erg / CONST_CGSM_ELECTRON_VOLT);

   % erg /s /cm^2 /GeV
   return f * energy_erg;
}

% Analytic result from Blumenthal and Gould (1970)

define a_bg (s)
{
   variable num, den;
   num = 2^((s-1.0)/2) * sqrt(3.0)
     * gamma ((3*s-1.0)/12) * gamma ((3*s+19.0)/12) * gamma ((s+5.0)/4);
   den = 8 * sqrt(PI) * (s+1.0) * gamma ((s+7.0)/4);
   return num/den;
}

define sync_analytic_bg (nu, s, B)
{
   % cgs units

   variable m_e = CONST_CGSM_MASS_ELECTRON;
   variable c   = CONST_CGSM_SPEED_OF_LIGHT;
   variable h   = CONST_CGSM_PLANCKS_CONSTANT_H;
   variable eV  = CONST_CGSM_ELECTRON_VOLT;

   variable electron_rest_energy = m_e * c^2;
   variable electron_charge = 4.803250e-10; % esu

   variable GeV = 1.e9 * eV;
   variable B_tot = B * 1.e-6;

   variable k = (GeV / electron_rest_energy)^s;
   variable c1 = 4*PI * k * electron_charge^3 / electron_rest_energy;
   variable c2 = 3*electron_charge / (4*PI * m_e * c);

   % erg/ s /cm^2 /Hz
   variable f = c1 * a_bg(s) * B_tot / (nu/(c2*B_tot))^((s-1.0)/2);

   % detailed code using a different value of the fine-structure
   % constant -- correct for that:
   variable alpha_gsl = 7.297352533e-3;
   variable alpha_def = electron_charge^2 / (h/(2*PI)) / c;
   f *= alpha_gsl / alpha_def;

   variable fixup_units = electron_rest_energy / h;
   return f * fixup_units;
}

define plot_fixed_B (B)
{
   % Hz
   variable nu = 10.0^[7.0:19.5:0.125];  %[9.0]; %

   xlabel (latex2pg ("\\nu [Hz]"));
   ylabel (latex2pg ("1 - sync / S_{analytic}(\\nu)"));
   title (latex2pg (sprintf ("\Gamma = [1 - 5], B = %g \mu G"R, B)));

   xrange (nu[0], nu[-1]);
   xlog;
   ylog; yrange (1.e-16, 1.e-10);
   %yrange (-4.e-5, 4.e-5);

   plot_box;
   color(1);

   variable fa, fn, s;

#iffalse
   foreach ([4.5])
#else     
   foreach ([1.0:4.0:0.5])
#endif       
     {
        s = ();
        fa = sync_analytic_bg (nu, s, B);
        fn = sync_numeric (nu, s, B);
        vmessage ("err_rms= %12.4g  err_max=%12.4g Gamma=%3.2f  B=%12.3g",
                  sqrt(mean((1-fn/fa)^2)),
                  max(abs(fn - fa)/fa),
                  s,
                  B);
        oplot (nu, abs((fa - fn)/fa));        
     }

}

define main()
{
   variable id = plot_open ("analytic_sync.ps/cps");

#iffalse
   foreach ([4.5])
#else
   foreach ([0.0:4.5:0.5])
#endif       
     {
        variable x = ();
        plot_fixed_B (10.0^x);
     }

   plot_close (id);
}

main();
