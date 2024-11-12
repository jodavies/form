#ifndef `TEST'
  #message Use -D TEST=XXX
  #terminate
#else
  #include `NAME_' # `TEST'
#endif
.end

*--#[ factorin_ :
CFunction f;
Symbol x1,x2,x3;
Vector p,q;
Index mu;
Local F = 10*x1*p.q*p(mu)*f(x1)*(x1+x2+x3/3);
Local F2 = 10*x1*p.q*p(mu)*f(x1)*(x1+x2+x3/3);
#$f = 10*x1*p.q*p(mu)*f(x1)*(x1+x2+x3/3);
.sort
Drop F2;

Local G = factorin_(F);
Local G2 = factorin_(F2);
Local g = factorin_($f);

Print;
.end
assert succeeded?
assert result("G") =~ expr("10/3*f(x1)*p(mu)*p.q*x1")
assert result("G2") =~ expr("10/3*f(x1)*p(mu)*p.q*x1")
assert result("g") =~ expr("10/3*f(x1)*p(mu)*p.q*x1")
*--#] factorin_ :
