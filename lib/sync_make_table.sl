private define fcn (x)
{
   return _sync_angular_integral (x, 0);
}

private define make_entry (x)
{
   variable t = struct {x, y, next};
   t.x = x;
   t.y = fcn(x);
   %() = fprintf (stderr, "%17.8e %17.8e\r", log10(x), log10(t.y));
   return t;
}

private define init_list (a, b)
{
   variable t = make_entry (a);
   t.next = make_entry (b);
   return t;
}

private define add_entry_after (t)
{
   variable s = make_entry (sqrt (t.x * t.next.x));
   s.next = t.next;
   t.next = s;
   return t;
}

private define refine_table (list, x_tol, y_tol)
{
   variable t = list;
   variable added_entries = 0;

   forever
     {
        variable split = 0;

	if (orelse
            {t.next == NULL}
            {t.x >= t.next.x})
          break;

        if (abs(1.0 - t.x/t.next.x) > x_tol)
          {
             split = 1;
          }
        else
          {
             variable ya, yb;

             ya = fcn (sqrt(t.x * t.next.x));
             yb = sqrt(t.y * t.next.y);

             if (abs(1.0 - yb/ya) > y_tol)
               {
                  split = 1;
               }
          }

        if (split)
          {
             t = add_entry_after (t);
             added_entries = 1;
          }

        t = t.next;
     }

   return added_entries;
}

private define make_table (x_min, x_max, x_tol, y_tol)
{
   variable t = init_list (x_min, x_max);

   forever
     {
        if (0 == refine_table (t, x_tol, y_tol))
          break;
     }

   variable n, num = 0;
   foreach (t)
     {
        n = ();
        num++;
     }

   variable x, y;
   x = Double_Type[num];
   y = Double_Type[num];

   num = 1;
   x[0] = t.x;
   y[0] = t.y;
   foreach (t.next)
     {
        n = ();
        if (n.x == x[num-1])
          continue;
        x[num] = n.x;
        y[num] = n.y;
        num++;
     }

   variable s = struct {x,y};
   s.x = log10(x);
   s.y = log10(y);

   return s;
}

define _sync_make_table ()
{
   variable x_min = 1.0e-38;
   variable x_max = 100.0;
   variable x_tol = 0.05;
   variable y_tol = 1.25e-5;
   return make_table (x_min, x_max, x_tol, y_tol);
}

provide ("sync_make_table");
