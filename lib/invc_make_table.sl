require ("nonthermal");

% jch "borrowed" this bisection routine from John Davis
%!%+
%\function{bisection}
%\synopsis{Find a root using bisection}
%\usage{root = bisection(Ref_Type f, Double_Type a, Double_Type b, ...)}
%\description
%  The \fun{bisection} function computes a root of the function \exmp{f(x,...)}
%  between in the interval \exmp{[a,b]} using the bisection algorithm.
%\notes
%  The interval \exmp{[a,b]} must be such that \exmp{f(a)*f(b)<0} to ensure that
%  a root exists in the interval.
%
%  This algorithm uses a combination bisection and the so-called false
%  projection technique to bracket a root.
%\seealso{newton}
%!%-
private define bisection ()
{
   variable f, a, b;
   variable client_data = __pop_args (_NARGS-3);
   (f,a,b) = ();

   variable count = 1;
   variable bisect_count = 5;
   variable fb, fa;
   variable x;

   if (a > b)
     (a,b) = (b,a);

   a *= 1.0;
   b *= 1.0;

   fa = @f(a, __push_args(client_data));
   fb = @f(b, __push_args(client_data));

   if (fa * fb > 0)
     verror ("bisection: the interval may not bracket a root: f(a)*f(b)>0");

   while (b > a)
     {
	if (fb == 0)
	  return b;
	if (fa == 0)
	  return a;

	if (count mod bisect_count)
	  {
	     x = (a*fb - b*fa) / (fb - fa);
	     if ((x <= a) or (x >= b))
	       x = 0.5 * (a + b);
	  }
	else			       %  bisect
	  x = 0.5 * (a + b);

	if ((x <= a) or (x >= b))
	  break;

	variable fx = @f(x, __push_args (client_data));
	count++;

	if (fx*fa < 0)
	  {
	     fb = fx;
	     b = x;
	  }
	else
	  {
	     fa = fx;
	     a = x;
	  }
     }

   return x;
}

private variable Special_Log= -50.0;
private variable Special_Val = 10.0^Special_Log;

private variable CBR_Temp = 2.725;
private variable ephoton_scale = Const_eV / (Const_m_e * Const_c^2);

% Maximum range of interest
private variable Log_X_Range = [1.0, 9.5];    % gamma
private variable Log_Y_Range = log10(ephoton_scale * 10.0^[2.0,15.0]);  % efinal
private variable X_Epsilon = 1.e-5 * (Log_X_Range[1] - Log_X_Range[0]);
private variable Sigma0 = 1.e-27;

private variable Xleft, Yleft;

private define log_to_lglg (xlg, ylg)
{
   variable xlg_left = interpol (ylg, Yleft, Xleft);
   variable xx = log (xlg - xlg_left);
   return xx;
}

private define lglg_to_log (xx, ylg)
{
   variable xlg_left = interpol (ylg, Yleft, Xleft);
   variable xlg = exp (xx) + xlg_left;
   return xlg;
}

private define safe_log10 (x)
{
   variable lg;

   if (typeof(x) != Array_Type)
     {
        if (x > Special_Val)
          {
             lg = log10(x);
             return lg;
          }
        return Special_Log;
     }

   lg = @x;
   variable i = where (x > Special_Val);
   lg[i] = log10(x[i]);

   % avoid comparing floats for equality
   variable j = [0:length(x)-1];
   j[i] = -1;
   lg[where (j>=0)] = Special_Log;

   return lg;
}

private define fcn (xlg, ylg)
{
   variable gam, ef;
   gam = 10.0^xlg;
   ef = 10.0^ylg;

   variable f;
   if (length(gam) == 1)
     {
        f = _invc_photon_integral(gam, ef, CBR_Temp, 0);
     }
   else
     {
        f = array_map (Double_Type, &_invc_photon_integral, gam, ef, CBR_Temp, 0);
     }

   return safe_log10 (f / Sigma0);
}

private define _fcn_lglg (xx, ylg)
{
   variable xlg = lglg_to_log (xx, ylg);
   variable f = fcn (xlg,ylg);
   return f;
}

private define fcn_lglg (xx, ylg)
{
   return array_map (Double_Type, &_fcn_lglg, xx, ylg);
}

private variable Y_value;
private variable Smallest_Useful_F;

private define rootfun (xlg)
{
   return fcn (xlg, Y_value) - Smallest_Useful_F;
}

private define find_root (ylg)
{
   variable _x = [Log_X_Range[0]:Log_X_Range[1]:0.1];
   variable _f = fcn(_x, ylg);
   variable imax = where (max(_f) == _f)[0];

   Y_value = ylg;
   Smallest_Useful_F = max(_f) - 12.0;  % log-scale

   variable ylg_max = ylg + 0.2 * (Log_Y_Range[1]-Log_Y_Range[0]);
   if (ylg_max > Log_Y_Range[1])
     ylg_max = Log_Y_Range[1];

   return bisection (&rootfun, Log_X_Range[0], _x[imax]);
}

private define save_boundary (fp, x, y)
{
   variable s = struct {x, y};
   s.x = x;
   s.y = y;

   fits_write_binary_table (fp, "BOUNDARY", s);
   fits_update_key (fp, "gammin", Log_X_Range[0], "log10(Gamma_min)");
   fits_update_key (fp, "gammax", Log_X_Range[1], "log10(Gamma_max)");
   fits_update_key (fp, "efnmin", Log_Y_Range[0], "log10(Efinal_min)");
   fits_update_key (fp, "efnmax", Log_Y_Range[1], "log10(Efinal_max)");
   fits_update_key (fp, "sigma0", Sigma0, "sigma0");
   fits_update_key (fp, "xepsilon", Sigma0, "xepsilon");
}

private define boundary_corners ()
{
   variable fp, file = "invc_bdry.fits";

   if (NULL == stat_file (file))
     {
        variable num = 1024;
        vmessage ("Finding min(gamma) vs. E_final...[%d points]", num);
        variable t = [0:num-1]/(num*1.0);  % don't want t[-1]=1
        Yleft = Log_Y_Range[0] + (Log_Y_Range[1]-Log_Y_Range[0])*t;
        Xleft = array_map (Double_Type, &find_root, Yleft);
        vmessage ("writing %s", file);
        fp = fits_open_file (file, "c");
        save_boundary (fp, Xleft, Yleft);
        fits_close_file (fp);
        vmessage ("... done");
        exit(0);
     }
   else
     {
        variable s = fits_read_table (file+"[BOUNDARY]");
        Xleft = s.x;
        Yleft = s.y;
     }

   variable xx, xlg, ylg;
#iffalse
   variable xright = ones(length(Yleft)) * (Log_X_Range[1] + 0.5);
   xlg = [Xleft + X_Epsilon, xright];
   ylg = [Yleft,             Yleft];
   xx = log_to_lglg (xlg, ylg);
#else   
   xlg = Xleft + X_Epsilon;
   ylg = Yleft;
   xx = log_to_lglg (xlg, ylg);
   xx = [xx, ones(length(Yleft))*3.0];
   ylg = [ylg, Yleft];
#endif   

   return (xx, ylg);
}

private define compute_table (file)
{
   variable x, y;
   (x,y) = boundary_corners ();

   variable l = [0:length(x)/2-1];
   variable r = l + length(x)/2;

   variable xleft, xright, yy;
   xleft = x[l];
   xright = x[r];
   yy = y[l];

   variable nx = 512;
   variable ny = 512;

   variable xmn, xmx, xgrid;
   xmn = min(xleft);
   xmx = max(xright);
   xgrid = xmn + (xmx - xmn) * [0:nx-1]/(nx-1.0);

   variable ymn, ymx, ygrid;
   ymn = Log_Y_Range[0];
   ymx = Log_Y_Range[1];
   ygrid = ymn + (ymx - ymn) * [0:ny-1]/(ny-1.0);

   variable xl, xr, i, j, xlg, ylg, xx;
   variable f = Double_Type[ny,nx];

   _for (0, ny-1, 1)
     {
        j = ();
        print(j);
        xl = interpol (ygrid[j], yy, xleft);
        xr = interpol (ygrid[j], yy, xright);
        i = where (xl <= xgrid and xgrid < xr);
        if (length(i) == 0)
          continue;
        f[j,i] = fcn_lglg (xgrid[i], ygrid[j]);
     }

   variable keys = struct
     {
        gammin, gammax, efnmin, efnmax, xepsilon, sigma0
     };
   keys.gammin = Log_X_Range[0];
   keys.gammax = Log_X_Range[1];
   keys.efnmin = Log_Y_Range[0];
   keys.efnmax = Log_Y_Range[1];
   keys.sigma0 = Sigma0;
   keys.xepsilon = X_Epsilon;

   variable bdry = struct
     {
        xleft, yleft
     };
   bdry.yleft = Yleft;
   bdry.xleft = Xleft;

   variable t = struct
     {
        xgrid, ygrid, f
     };
   t.xgrid = xgrid;
   t.ygrid = ygrid;
   t.f = f;

   variable fp = fits_open_file (file, "c");
   fits_write_binary_table (fp, "TABLE", t, keys, NULL);
   fits_write_binary_table (fp, "BOUNDARY", bdry, keys, NULL);
   fits_close_file(fp);
}

define _invc_make_table (file)
{
   compute_table(file);
}

provide ("invc_make_table");
