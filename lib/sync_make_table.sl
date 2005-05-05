static define make_entry (x)
{
   variable t = struct {x, y, next, prev};
   t.x = x;
   t.y = _sync_angular_integral(x, 0);
   %() = fprintf (stderr, "%17.8e %17.8e\r", log10(x), log10(t.y));
   return t;
}

static define init_list (a, b)
{
   variable t = make_entry (a);
   t.next = make_entry (b);
   return t;
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

static define need_new_entry (list, tol)
{
   foreach (list)
     {
	variable t = ();

	if (orelse
            {t.next == NULL}
            {t.x >= t.next.x})
	  return NULL;

        % don't let y change too much
        if (abs(t.next.y - t.y) > tol * abs(t.y))
	  return t;

        % don't let x change too much
        if (abs(t.next.x - t.x) > tol * abs(t.next.x))
          return t;
     }

   return NULL;
}

static define make_table (x_min, x_max, tol)
{
   variable t = init_list (x_min, x_max);
   variable num = 2;

   variable n = t;
   forever
     {
        n = need_new_entry (n, tol);
        if (n == NULL)
          break;
        n = add_entry_after (n);
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
   variable x_min = 1.e-38;
   variable x_max = 32.0;
   variable y_tol = 0.05;
   return make_table (x_min, x_max, y_tol);
}

provide ("sync_make_table");
