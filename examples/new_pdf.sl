require ("nonthermal");

variable lib = path_concat (path_dirname(__FILE__), "new_pdf.so");

add_pdf (lib, "simple", ["index"],   % file, name, param_names
         [2.0], [0], [1.0], [3.0]);  % value, freeze, min, max
