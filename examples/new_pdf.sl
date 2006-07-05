require ("nonthermal");

variable lib = path_concat (path_dirname(__FILE__), "new_pdf.so");
add_user_pdf_intrin (lib, "simple", "");

private define pdf_param_defaults (i, val, freeze, min, max)
{
   return (val[i], freeze[i], min[i], max[i]);
}

private define pdf_simple (l,h,p) {return ("simple", p);}
private define pdf_simple_contin (x,p) {return ("simple", p);}

add_slang_function ("pdf_simple", [&pdf_simple, &pdf_simple_contin], ["index"]);
set_param_default_hook ("pdf_simple", &pdf_param_defaults,
                        [2.0], [0], [1.0], [3.0]);
