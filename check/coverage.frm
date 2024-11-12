#ifndef `TEST'
  #message Use -D TEST=XXX
  #terminate
#else
  #include `NAME_' # `TEST'
#endif
.end

*--#[ if :
Symbol x(:100),y,z;
CFunction f,g(symmetric),h;
Vector p,q;
Tensor T,G(symmetric);
Index mu,i1,...,i8;

Set hset: h;

Local test = x * (
	+ z^2
	+ z^3
	+ f(x)
	+ f(y)
	+ f(y)*x
	+ p(mu)
	+ p
	+ p.q
	+ T(p)*z
	+ f(p)
	+ g(p)
	+ h(p)
	+ g(i1,i2)*g(i1,i2)
	+ g(i1,i2)*g(i1,i3)
	+ g(i1,i2,i3)*g(i1,i2,i3)
	+ <g(i1,i2)>*...*<g(i7,i8)>*g(i8,i1)
	+ G(i1,i2)*G(i1,i2)
	+ G(i1,i2)*G(i1,i3)
	+ G(i1,i2,i3)*G(i1,i2,i3)
	+ <G(i1,i2)>*...*<G(i7,i8)>*G(i8,i1)
	);
Local test2 = x;
.sort

#$one = 1;
#$two = 2;
#$pdotq = p.q;

#if ( (maxpowerof(x) > 100) || (minpowerof(x) > -100) )
	Multiply z;
#endif

#if ( termsin(test) == 15 )
	Multiply 2;
#endif

#if ( exists(test3) )
	Multiply 1000;
#endif

#if (isfactorized(test2) == 0)
	If (Expression(test2));
		Identify x = 100;
	EndIf;
#endif

If (Match(f(x?)*x?));
	Identify x = 1;
EndIf;

If (Count(z,1) == 3);
	Identify x = 1;
ElseIf (Count(z,1) == $two);
	Identify x = 2;
EndIf;

If (Count(p+v,1) == 1);
	Identify x = 1;
EndIf;
If (Count(p+d,1) == 1);
	Identify x = 2;
EndIf;
If (Count(p+f,1) == 1);
	Identify x = 3;
EndIf;
If (Count(p+?{f,g},1) == 1);
	Identify x = 4;
EndIf;
If (Count(p+?hset,1) == 1);
	Identify x = 5;
EndIf;
If (Occurs(p.q));
	Multiply 6;
EndIf;
If (Match($pdotq));
	Multiply 6;
EndIf;

If (Occurs(f) || Match(g(p)));
	Multiply 2;
	Identify x = 2;
EndIf;

If ( (Coefficient == MultipleOf(3)) && (Count(z,1) == $one) );
	Multiply z;
EndIf;

If (FindLoop(g, arguments=2, loopsize=2));
	Identify x = 2;
EndIf;
If (FindLoop(g, arguments=3, loopsize=2));
	Identify x = 3;
EndIf;
If (FindLoop(g, arguments=3, loopsize<3, include=i3));
	Multiply 4;
EndIf;
If (FindLoop(g, arguments=2, loopsize=8));
	Multiply 5;
EndIf;

If (FindLoop(G, arguments=2, loopsize=2));
	Identify x = 2;
EndIf;
If (FindLoop(G, arguments=3, loopsize=2));
	Identify x = 3;
EndIf;
If (FindLoop(G, arguments=3, loopsize<3, include=i3));
	Multiply 4;
EndIf;
If (FindLoop(G, arguments=2, loopsize=8));
	Multiply 5;
EndIf;

If ($two);
	Multiply 2;
EndIf;

If (1 >= $two);
	Multiply 1000;
EndIf;

If (1 < 2);
	Multiply 2;
EndIf;
If (1 <= 2);
	Multiply 2;
EndIf;

If (1 > 2);
	Multiply 1000;
EndIf;
If (1 >= 2);
	Multiply 1000;
EndIf;

If (0);
	Multiply 1000;
EndIf;

Print;
.end
assert succeeded?
assert result("test") =~ expr("16*z^2 + 8*z^3 + 576*p.q + 8*p(mu) + 8*p + 16*f(x) + 96*f(y) + 64*f(p) + 64*g(p) + 16*g(i1,i2)^2 + 8*g(i1,i2)*g(i1,i3)*x + 40*g(i1,i2)*g(i1,i8)*g(i2,i3)*g(i3,i4)*g(i4,i5)*g(i5,i6)*g(i6,i7)*g(i7,i8)*x + 96*g(i1,i2,i3)^2 + 40*h(p) + 24*T(p)*z^2 + 96*G(i1,i2,i3)*G(i1,i2,i3) + 16*G(i1,i2)*G(i1,i2) + 8*G(i1,i2)*G(i1,i3)*x + 40*G(i1,i2)*G(i1,i8)*G(i2,i3)*G(i3,i4)*G(i4,i5)*G(i5,i6)*G(i6,i7)*G(i7,i8)*x")
assert result("test2") =~ expr("800")
*--#] if :
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
*--#[ hideskip :
#: scratchsize 1
#: hidesize 1
Symbol x,y,z;
Local f = (1+x+y+z)^70;
Local g = (1+x+y+z)^70;
.sort
Hide f;
Skip g;
Identify x = 2;
Identify y = 2;
Identify z = 2;
.sort
UnHide f;
Identify x = 1;
Identify y = 1;
Identify z = 1;
Print;
.end
assert succeeded?
assert result("f") =~ expr("1393796574908163946345982392040522594123776")
assert result("g") =~ expr("1393796574908163946345982392040522594123776")
*--#] hideskip :
*--#[ saveload :
Symbol x;
Function f;
CFunction g;
Index mu;
Vector p,q;
Tensor T;
Global F = (1+x)^20 + f(1) + g(2) + p.q + p(mu) + T(q) + T(mu) + e_(p,q);
.store
Save test.sav F;
.end
Load test.sav;
Local G = F;
.sort
Delete storage;
Identify x = 1;
Print;
.end
assert succeeded?
assert result("G") =~ expr("1048576 + p.q + p(mu) + e_(p,q) + g(2) + T(q) + T(mu) + f(1)")
*--#] saveload :
*--#[ symm :
* pattern matching of functions with symmtry properties
Indices x1,x2,x3,x4,x5;

CFunction fs(symmetric), fas(antisymmetric), fcs(cyclesymmetric);
Function gs(symmetric), gas(antisymmetric), gcs(cyclesymmetric);
CTensor Ts(symmetric), Tas(antisymmetric), Tcs(cyclesymmetric);
Tensor Us(symmetric), Uas(antisymmetric), Ucs(cyclesymmetric);

Local test =
	+ fs(1,2,3,4,5) + fas(1,2,3,4,5) + fcs(1,2,3,4,5)
	+ gs(1,2,3,4,5) + gas(1,2,3,4,5) + gcs(1,2,3,4,5)
	+ Ts(1,2,3,4,5) + Tas(1,2,3,4,5) + Tcs(1,2,3,4,5)
	+ Us(1,2,3,4,5) + Uas(1,2,3,4,5) + Ucs(1,2,3,4,5)
	;

Identify fs(x1?,x2?,4,3,x5?) = fs;
Identify fas(x1?,x2?,4,3,x5?) = - fas;
Identify fcs(5,1,x3?,x4?,x5?) = fcs;

Identify gs(x1?,x2?,4,3,x5?) = gs;
Identify gas(x1?,x2?,4,3,x5?) = - gas;
Identify gcs(5,1,x3?,x4?,x5?) = gcs;

Identify Ts(x1?,x2?,4,3,x5?) = Ts;
Identify Tas(x1?,x2?,4,3,x5?) = - Tas;
Identify Tcs(5,1,x3?,x4?,x5?) = Tcs;

Identify Us(x1?,x2?,4,3,x5?) = Us;
Identify Uas(x1?,x2?,4,3,x5?) = - Uas;
Identify Ucs(5,1,x3?,x4?,x5?) = Ucs;

Print;
.end
assert succeeded?
assert result("test") =~ expr("fs + fas + fcs + Ts + Tas + Tcs + Us + Uas + Ucs + gs + gas + gcs")
*--#] symm :
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
*--#[ transform :
CFunction Z,f,g,h,i;
Symbol inf,x,y,y1,y2,n1,n2;
Set funs: f,g;
#$one = 1;
#$two = 2;

Local test1 = f(1,2,3,4) + g(1,2,3,4) + h(1,2,3,4) + i(1,2,3,4);
Local test2 = f(1,2,3,4);
Local test3 = f(1,2,3,4,5,6,7);
Local test4 = f(1,2,3,4);
Local test5 = f(1,2,3,4) + g(1,2,3,4);
Local test6 = f(1) + f(1,1) + f(1,2) + f(2,1) + f(1,3,2) + f(2,3,1);
Local test7 = test6;
Local test8 =
	+ g(1)*f(1,1,1)
	+ g(2)*f(-1,1,1)
	+ g(3)*f(2,1,1)
	+ g(4)*f(1,-2,1)
	+ g(5)*f(-1,-1,-1)
*	Check against HarmonicSums.m
	- (
		+ g(1)*Z(1, 1, 1, inf)
		- g(2)*Z(-1, -1, 1, inf)
		+ g(3)*Z(2, 1, 1, inf)
		+ g(4)*(-Z(1, -2, -1, inf))
		- g(5)*Z(-1, 1, 1, inf)
	)
	;
Local test9 = test8;

InExpression test1;
	Transform funs addargs($one,last);
	Transform {h,i}, mulargs($one,last);
EndInExpression;

InExpression test2;
	Transform f explode(1,last);
	Identify f(?a) = f(?a)*g(nargs_(?a));
	Identify g(y?$len) = 1;
	Transform f encode(1,last):base=$two;
	Transform f decode(1,$len):base=2, implode(1,last);
EndInExpression;

InExpression test3;
	Transform f permute($one,3,5)($two,6);
EndInExpression;

InExpression test4;
	Transform f reverse($two,last);
EndInExpression;

InExpression test5;
	Transform f cycle($two,last)=+1;
	Transform g cycle(2,last)=-$one;
EndInExpression;

InExpression test6;
	Transform f islyndon(1,last)=(y1,n1);
	Transform f islyndon<(1,last)=(y1,n1);
	Transform f islyndon-(1,last)=(y1,n1);
	Transform f islyndon>(1,last)=(y2,n2);
	Transform f islyndon+(1,last)=(y2,n2);
EndInExpression;

InExpression test7;
	Transform f tolyndon(1,last)=(y1,n1);
	Transform f tolyndon<(1,last)=(y1,n1);
	Transform f tolyndon-(1,last)=(y1,n1);
	Transform f tolyndon>(1,last)=(y2,n2);
	Transform f tolyndon+(1,last)=(y2,n2);
EndInExpression;

InExpression test8;
	Transform f HtoZ(1,last);
	Identify f(?a) = Z(?a,inf);
EndInExpression;

InExpression test9;
	Transform Z ZtoH(1,last-1);
	Identify Z(?a,inf) = f(?a);
EndInExpression;

Print;
.end
assert succeeded?
assert result("test1") =~ expr("f(10) + g(10) + h(24) + i(24)")
assert result("test2") =~ expr("f(1,2,3,4)")
assert result("test3") =~ expr("f(3,6,5,4,1,2,7)")
assert result("test4") =~ expr("f(1,4,3,2)")
assert result("test5") =~ expr("f(1,4,2,3) + g(1,3,4,2)")
assert result("test6") =~ expr("f(1)*y1^3*y2^2 + f(1,1)*n1^3*n2^2 + f(1,2)*y1^3*n2^2 + f(1,3,2)*y1^3*n2^2 + f(2,1)*y2^2*n1^3 + f(2,3,1)*n1^3*n2^2")
assert result("test7") =~ expr("f(1)*y1^3*y2^2 + f(1,1)*n1^3*n2^2 + 2*f(2,1)*y1^3*y2^2 + f(3,1,2)*y1^3*y2^2 + f(3,2,1)*y1^3*y2^2")
assert result("test8") =~ expr("0")
assert result("test9") =~ expr("0")
*--#] transform :
