
_auto_declare=1;

(x,y,f) = readcol ("angfun.dat",1, 2, 3);

xlabel (latex2pg("log x"));
ylabel (latex2pg("log I(x)"));

id = plot_open ("angfun.ps/cps");
_pgslw(2);

multiplot([1,1]);
plot(x,y);

ylabel (latex2pg("I(x)/F(x)"));
yrange (,0.9);
plot(x, 10.0^(y - f));

print(max(10.0^(y-f)));

plot_close(id);

