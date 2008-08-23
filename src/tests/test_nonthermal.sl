%
if (NULL != stat_file ("nonthermal-module.so"))
{
   $1 = getcwd();
   prepend_to_isis_module_path ($1);
   prepend_to_isis_load_path (path_concat ($1, "../lib"));
}
%

require ("nonthermal");
require ("readascii");

define log_egrid (e_min_eV, e_max_eV, n)
{
   variable elo, ehi;
   (elo, ehi) = linear_grid (log10(e_min_eV), log10(e_max_eV), n);
   elo = 10.0^elo;
   ehi = 10.0^ehi;
   return elo, ehi;
}

define new_test (spec, spec_pars, pdf, pdf_pars, elo_eV, ehi_eV)
{
   variable t = struct
     {
        spec, spec_pars, pdf, pdf_pars,
        elo_eV, ehi_eV, f
     };

   t.spec = spec;
   t.spec_pars = spec_pars;
   t.pdf = pdf;
   t.pdf_pars = pdf_pars;
   t.elo_eV = elo_eV;
   t.ehi_eV = ehi_eV;

   return t;
}

define eval_fun_ev (elo_ev, ehi_ev)
{
   return reverse (eval_fun (_A(elo_ev/1.e3, ehi_ev/1.e3)));
}

define compute_spectrum (t)
{
   variable spec = t.spec, pdf = t.pdf;
   fit_fun ("$spec((1), $pdf(1))"$);
   if (t.spec_pars != NULL)
     {
        set_par ("$spec(1)"$, t.spec_pars);
     }
   if (t.pdf_pars != NULL)
     {
        set_par ("$pdf(1)"$, t.pdf_pars);
     }
   t.f = eval_fun_ev (t.elo_eV, t.ehi_eV);

#iffalse
   xlog;
   ylog;
   hplot (t.elo_eV, t.ehi_eV, t.f);
   plot_pause;
#endif
}

define save_test (t, dir, id)
{
   variable prefix = sprintf ("%s_%s_%d", t.spec, t.pdf, id);

   vmessage(" ===> ${dir}/${prefix}.p"$);
   save_par ("${dir}/${prefix}.p"$);

   variable fp = fopen ("$dir/$prefix.dat"$, "w");
   if (fp == NULL)
     throw ApplicationError, "*** Failed opening $prefix.dat"$;

   variable i, n = length(t.elo_eV);
   _for i (0, n-1, 1)
     {
        () = fprintf (fp, "%0.17e %0.17e %0.17e\n",
                      t.elo_eV[i], t.ehi_eV[i], t.f[i]);
     }

   if (-1 == fclose (fp))
     throw ApplicationError, "*** Failed closing $prefix.dat"$;
}

private variable Test_Counter = 0;

define generate_tests (dir, tests)
{
   if ((dir != NULL) && (NULL == stat_file(dir)))
     {
        if (0 != mkdir (dir, 0777))
          {
             if (errno != EEXIST)
               throw ApplicationError, "Failed creating directory '$dir'"$;
          }
     }

   variable t;
   foreach t (tests)
     {
        tic;
        compute_spectrum (t);
        vmessage ("%S sec:   %s", toc, get_fit_fun());
        if (dir != NULL)
          {
             Test_Counter += 1;
             save_test (t, dir, Test_Counter);
          }
     }
}

define perform_tests (dir)
{
   variable pf, pfiles = glob ("$dir/*.p"$);

   vmessage ("\nRunning %d tests - this may take a while...",
             length(pfiles));

   variable dt, n = 0;

   foreach pf (pfiles)
     {
        variable test_file, i, err, elo, ehi, tf, f;

        (test_file, ) = strreplace (pf, ".p", ".dat", 1);

        if (readascii (test_file, &elo, &ehi, &tf; format="%le %le %le") < 0)
          throw ApplicationError, "Failed reading $test_file"$;

        load_par (pf);
        tic;
        f = eval_fun_ev (elo, ehi);
        dt = toc();

        i = where (tf != 0);
        err = 1.0 - f[i]/tf[i];

        n += 1;
        vmessage ("%2d. err: mean=%9.3e max=%9.3e [%7.3f sec] %s",
                  n, mean(abs(err)), max(abs(err)), dt,
                  path_basename (test_file));
     }
}

define test_list (pdf, pdf_pars)
{
   variable tl = list_new();
   variable elo, ehi;

   (elo, ehi) = log_egrid (1.e-3, 1.e5, 1000);
   list_append (tl, new_test ("sync", [1.0, 1.0], pdf, pdf_pars, elo, ehi));
   list_append (tl, new_test ("sync", [1.0, 10.0], pdf, pdf_pars, elo, ehi));

   (elo, ehi) = log_egrid (1.e2, 1.e13, 1200);
   list_append (tl, new_test ("ntbrem", NULL, pdf, pdf_pars, elo, ehi));

   (elo, ehi) = log_egrid (1.e6, 1.e13, 1000);
   list_append (tl, new_test ("invc", NULL, pdf, pdf_pars, elo, ehi));
   list_append (tl, new_test ("pizero", NULL, pdf, pdf_pars, elo, ehi));

   return tl;
}

define isis_main ()
{
   variable dir = path_concat (path_dirname(__FILE__), "../../data/tests");

   if (__argc > 3)
     {
        usage ("%s DIR [init]", __argv[0]);
        exit(0);
     }

   if (__argc > 1)
     dir = __argv[1];

   if (__argc < 3 || __argv[2] != "init")
     {
        perform_tests (dir);
        exit(0);
     }

   message ("Initializing test cases in '$dir':"$);
   variable pdf = "full1";
   generate_tests (dir, test_list (pdf, [ 1.0, 1.0, 2.0, 0.0 , 10.0]));
   generate_tests (dir, test_list (pdf, [50.0, 1.0, 2.0, 0.0 , 10.0]));
   generate_tests (dir, test_list (pdf, [ 1.0, 0.1, 2.0, 0.0 , 10.0]));
   generate_tests (dir, test_list (pdf, [ 1.0, 1.0, 2.2, 0.0 , 10.0]));
   generate_tests (dir, test_list (pdf, [ 1.0, 1.0, 2.0, 0.05, 10.0]));
   generate_tests (dir, test_list (pdf, [ 1.0, 1.0, 2.0, 0.0 , 30.0]));
}
