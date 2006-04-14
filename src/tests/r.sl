require ("gsl");
require ("nonthermal");

define W (x,a,b)
{
   return exp(-x/2.0) * x^(b + 0.5) * hyperg_U (0.5+b-a, 1+2*b, x);
}

define R(x)
{
   return 0.5*PI*x * (W (x, 0, 4.0/3) * W (x, 0, 1.0/3)
                      - W (x, 0.5, 5.0/6) * W (x, -0.5, 5.0/6));
}

define uu(kappa,mu,x)
{
   return hyperg_U (0.5+mu-kappa, 1.0+2*mu, x);
}

define R2(x)
{
   return 0.5*PI*exp(-x)*x^(11.0/3) *
     (uu(0,4.0/3,x) * uu(0, 1.0/3, x) - uu(0.5, 5.0/6, x)*uu(-0.5, 5.0/6, x));
}

define R_asymptotic(x)
{
   variable i, j, k;
   variable r = Double_Type[length(x)];
   
   variable xlo, xhi;
   xlo = 1.e-21;
   xhi = 30.0;
   
   % These asymptotic forms don't seem very accurate...

   i = where (x < xlo);
   if (length(i) > 0)
     {
        r[i] = 1.808418021102803e+00 * x[i]^(1.0/3);
     }   

   j = where (xlo <= x and x < xhi);
   if (length(j) > 0)
     r[j] = R2(x[j]);

   k = where (xhi <= x);
   if (length(k) > 0)
     r[k] = 0.5 * PI * exp(-x[k]) * (1.0 - 99.0/162.0/x[k]);

   return r;
}

define main()
{
   variable x = 10.0^[-38.0:3.0:0.1];
   variable r = R(x);

   xlog;
   ylog;
   plot (x, abs(1.0 - R_asymptotic(x)/R2(x)));

#iffalse
   variable y = array_map (Double_Type, &_sync_angular_integral, x, 0);
   variable yi = array_map (Double_Type, &_sync_angular_integral, x, 1);
   ylog;
   plot (x, abs(1.0 - yi/y));
#endif
}

main ();
