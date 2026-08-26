* Coverage-oriented tests for normal paths in sources/if.c.

*--#[ coverage_if_dollar_factor_number :
Symbol x;

* A factored dollar can provide its factor count to a run-time If.
#$factorcount = 1;
#factdollar $factorcount
Local factorcount = 1;
If ( $factorcount[0] == 1 );
  Multiply 2;
EndIf;
Print;
.end
assert succeeded?
assert result("factorcount") =~ expr("2")
*--#] coverage_if_dollar_factor_number :

*--#[ coverage_if_switch_dense :
Symbol x;

* Consecutive cases use the dense switch table and Break jumps to EndSwitch.
#$densecase = 2;
Local densecase = 1;
Switch $densecase;
Case 1;
  Multiply 10;
  Break;
Case 2;
  Multiply 20;
  Break;
Case 3;
  Multiply 30;
  Break;
Default;
  Multiply 100;
EndSwitch;
Print;
.end
assert succeeded?
assert result("densecase") =~ expr("20")
*--#] coverage_if_switch_dense :

*--#[ coverage_if_switch_sparse :
Symbol x;

* Unordered, widely-spaced cases use the sorted sparse table.  The second
* switch also checks the default lookup path.
#$sparsecase = -7;
#$defaultcase = 99;
Local sparsecase = 1;
Local defaultcase = 1;
.sort
Skip;
NSkip sparsecase;
Switch $sparsecase;
Case 500;
  Multiply 500;
  Break;
Case -7;
  Multiply 7;
  Break;
Case 10;
  Multiply 10;
  Break;
Default;
  Multiply 100;
EndSwitch;
.sort
Skip;
NSkip defaultcase;
Switch $defaultcase;
Case 500;
  Multiply 500;
  Break;
Case -7;
  Multiply 7;
  Break;
Case 10;
  Multiply 10;
  Break;
Default;
  Multiply 99;
EndSwitch;
.sort
Print;
.end
assert succeeded?
assert result("sparsecase") =~ expr("7")
assert result("defaultcase") =~ expr("99")
*--#] coverage_if_switch_sparse :

*--#[ coverage_if_occurs_variants :
Symbol x;
CFunction f;
Vector p,q;
Index mu;

* Exercise Occurs searches through symbols, functions, vector components,
* and dot products, including the reversed dot-product form.
Local occursvariants = f(x,p(mu))*p.q;
If ( Occurs(x) );
  Multiply 2;
EndIf;
If ( Occurs(f) );
  Multiply 3;
EndIf;
If ( Occurs(p) );
  Multiply 5;
EndIf;
If ( Occurs(q.p) );
  Multiply 7;
EndIf;
If ( Occurs(p.q) );
  Multiply 11;
EndIf;
Print;
.end
assert succeeded?
assert result("occursvariants") =~ expr("2310*f(x,p(mu))*p.q")
*--#] coverage_if_occurs_variants :

*--#[ coverage_if_occurs_nested_and_tensor :
CFunction f,g;
Tensor T;
Vector p;
Index mu,nu;

* A function used as a short argument and indices in vector/tensor objects
* take the recursive Occurs paths that a top-level function match skips.
*Local deepoccurs = g(f)*p(mu)*T(mu,nu);
Local deepoccurs = g(f)*p(mu)*T(mu,nu);
If ( Occurs(f) );
  Multiply 2;
EndIf;
* Note that here we will not find mu, since p is contracted with T and it is
* written in Schoonschip notation.
If ( Occurs(mu) );
  Multiply 3;
EndIf;
If ( Occurs(nu) );
  Multiply 5;
EndIf;
Print;
.end
assert succeeded?
assert result("deepoccurs") =~ expr("10*g(f)*T(p,nu)")
*--#] coverage_if_occurs_nested_and_tensor :
