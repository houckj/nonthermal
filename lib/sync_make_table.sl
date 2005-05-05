static define make_entry (x)
{
   variable t = struct {x, y, next, prev};
   t.x = x;
   t.y = _sync_angular_integral(x, 0);
   return t;
}

static define init_list (a, b)
{
   variable ta = make_entry (a);
   variable tb = make_entry (b);

   ta.next = tb;

   return ta;
}

static define add_entry_after (t)
{
   variable s = make_entry (sqrt (t.x * t.next.x));
   t.next.prev = s;
   s.next = t.next;
   s.prev = t;
   t.next = s;
   return t;
}

static define changing_too_fast (t, tol)
{
   return abs(t.next.y - t.y) > tol * abs(t.y);
}

static define need_new_entry (list, tol)
{
   foreach (list)
     {
	variable t = ();

	if (orelse
            {t.next == NULL}
            {t.x >= t.next.x})
	  return NULL;

        if (changing_too_fast (t, tol))
	  return t;
     }

   return NULL;
}

static define make_table (x_min, x_max, tol)
{
   variable t = init_list (x_min, x_max);
   variable n = t;

   variable num = 2;

   forever
     {
        n = need_new_entry (n, tol);
        if (n == NULL)
          break;
        n = add_entry_after (n);
%        () = fprintf (stderr, "%d %17.8e %17.8e\r", num, log10(n.x), log10(n.y));
        num++;
     }

   variable x, y;
   x = Double_Type[num];
   y = Double_Type[num];

   num = 0;
   foreach (t)
     {
        n = ();
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
   return make_table (1.e-38, 32.0, 1.e-2);
}

provide ("sync_make_table");
