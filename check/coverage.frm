#ifndef `TEST'
  #message Use -D TEST=XXX
  #terminate
#else
  #include `NAME_' # `TEST'
#endif
.end

*--#[ tablebase_1 :
CTable,sparse f(1);
#do i = 1,10
	Fill f(`i') = `i';
#enddo
.sort
TableBase "test.tbl" create;
TableBase "test.tbl" addto f;
.end
TableBase "test.tbl" open;
Tablebase "test.tbl" load;
.sort
Local test = <f(1)>+...+<f(11)>;
TestUse;
.sort
TableBase "test.tbl" use;
Apply;
Print;
.end
#pend_if mpi?
assert succeeded?
assert result("test") =~ expr("55 + f(11)")
*--#] tablebase_1 :
*--#[ tablebase_2 :
CTable,sparse f(1);
#do i = 1,10
	Fill f(`i') = `i';
#enddo
.sort
TableBase "test2.tbl" create;
TableBase "test2.tbl" addto f;
.end
TableBase "test2.tbl" open;
Tablebase "test2.tbl" enter;
.sort
Local test = <f(1)>+...+<f(11)>;
Apply;
Print;
.end
#pend_if mpi?
assert succeeded?
assert result("test") =~ expr("55 + f(11)")
*--#] tablebase_2 :
*--#[ tablebase_3 :
CTable,sparse f(1);
#do i = 1,10
	Fill f(`i') = `i';
#enddo
.sort
TableBase "test3.tbl" create;
TableBase "test3.tbl" addto f;
.end
TableBase "test3.tbl" open;
Tablebase "test3.tbl" audit;
.end
#pend_if mpi?
assert succeeded?
assert stdout =~ exact_pattern(<<'EOF')
Table,sparse,f(1)
    f(1)
    f(2)
    f(3)
    f(4)
    f(5)
    f(6)
    f(7)
    f(8)
    f(9)
    f(10)
EOF
*--#] tablebase_3 :
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
