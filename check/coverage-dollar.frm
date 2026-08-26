* Coverage-oriented tests for normal paths in sources/dollar.c.

*--#[ coverage_dollar_typed_names :
Function ff;
Tensor TT;
Vector vv;
Symbol x,y;
Index ii,jj;

* Dollar variables used in command parameters exercise the function, tensor
* and vector name conversions.
#$fun = ff;
#$ten = TT;
#$vec = vv;
Local chained = ff(x)*ff(y);
Local tensorized = vv(ii)*vv(jj);
ChainIn $fun;
ToTensor $vec,$ten;
Print;
.end
assert succeeded?
assert result("chained") =~ expr("ff(x,y)")
assert result("tensorized") =~ expr("TT(ii,jj)")
*--#] coverage_dollar_typed_names :

*--#[ coverage_dollar_preif :
Symbol x,y,z;

* The $(...) form evaluates complete dollar expressions before comparing them.
#$sum = x+y;
#$different = x+z;
#if $($sum == x+y)
  #$equal = 1;
#else
  #$equal = 0;
#endif
#if $($sum != $different)
  #$differentresult = 1;
#else
  #$differentresult = 0;
#endif
Local comparisons = $equal + 10*$differentresult;
Print;
.end
assert succeeded?
assert result("comparisons") =~ expr("11")
*--#] coverage_dollar_preif :

*--#[ coverage_factdollar_content :
Symbol x,y;

* FactDollar removes ordinary common content before factorization.
#$ordinary = x*(x+y) + x*(x+2*y);
#factdollar $ordinary
Local ordinary = $ordinary[1]*$ordinary[2];
Print;
.end
assert succeeded?
assert result("ordinary") =~ expr("3*x*y + 2*x^2")
*--#] coverage_factdollar_content :

*--#[ coverage_dollar_inside_and_queries :
Symbol x,y;

* Inside accepts a dollar expression and executes its normal processing path.
#$insidevalue = x+y;
Local trigger = 1;
Inside $insidevalue;
  Identify x = 2*x;
EndInside;
* TODO Local means the result is lost. Use of sum with a non-numeric $ is dodgy with tform.
ModuleOption sum $insidevalue;
.sort
Local insidevalue = $insidevalue;
Print;
.end
assert succeeded?
assert result("insidevalue") =~ expr("y+2*x")
*--#] coverage_dollar_inside_and_queries :

*--#[ coverage_dollar_output_and_preif_specials :
Symbol x,y;
Set choices: x,y;

* TODO These are not described in the manual!
#$member = x;
#$multiple = 2*x+2*y;
#$base = x+y;
#$factorprint = (1+x)*(1+y);
#if $($member == set_(choices))
  #$isset = 1;
#else
  #$isset = 0;
#endif
#if $($multiple == multipleof_($base))
  #$ismultiple = 1;
#else
  #$ismultiple = 0;
#endif
#if termsin($base) == 2
  #$termcountok = 1;
#else
  #$termcountok = 0;
#endif
#if sizeof($base) > 0
  #$sizeok = 1;
#else
  #$sizeok = 0;
#endif
* Dollar output also serializes a complete value and an individual factor.
#write "member is %$",$member
#factdollar $factorprint
#write "first factor is %$",$factorprint[1]
Local specialchecks = $isset + 10*$ismultiple + 100*$termcountok + 1000*$sizeok;
Print;
.end
assert succeeded?
assert result("specialchecks") =~ expr("1111")
assert stdout =~ exact_pattern("member is x")
assert stdout =~ exact_pattern("first factor is 1+y")
*--#] coverage_dollar_output_and_preif_specials :

*--#[ coverage_dollar_runtime_do :
Symbol x;

* A runtime Do loop evaluates dollar-valued bounds and its loop variable.
#$last = 3;
Local loopproduct = 1;
Do $i = 1,$last;
  Multiply $i;
EndDo;
ModuleOption local $i;
Print;
.end
assert succeeded?
assert result("loopproduct") =~ expr("6")
*--#] coverage_dollar_runtime_do :

*--#[ coverage_dollar_preprocessor_values :
Symbol x;

* The preprocessor loop reads its bounds directly from dollar variables.
* This exercises the integer conversion used by #do (rather than Do).
#$first = 2;
#$last = 4;
#$step = 1;
Local preprocessorloop = 1;
#do i = $first,$last,$step
  Multiply `i';
#enddo
Print;
.end
assert succeeded?
assert result("preprocessorloop") =~ expr("24")
*--#] coverage_dollar_preprocessor_values :

*--#[ coverage_dollar_exchange :
Symbol x;

* Exchanging two dollar names is a normal preprocessor operation; their
* values stay in place while their names are swapped.
#$left = 2;
#$right = 3;
#exchange $left,$right
Local exchanged = $left + 10*$right;
Print;
.end
assert succeeded?
assert result("exchanged") =~ expr("23")
*--#] coverage_dollar_exchange :

*--#[ coverage_dollar_typed_symbol_and_index :
#-
Off statistics;

Symbol x,y;
Index mu;
Vector p;

* Format Optimize accepts a symbol-valued dollar in a Horner scheme, while
* Trace4 accepts an index-valued dollar as its trace index.
#$scheme = x;
#$traceindex = mu;
Format O1,scheme=($scheme,y);
Local optimized = x*y + x^2;
Local traced = p(mu);
Trace4,$traceindex;
Print;
.sort
#optimize optimized
#clearoptimize
Print;
.end
assert succeeded?
assert stdout =~ exact_pattern(<<'EOF')
      Z1_=y + x;
      optimized=x*Z1_;

    Z1_=p(mu);
   traced=Z1_;


    Z1_=p(mu);
   traced=Z1_;
EOF
*--#] coverage_dollar_typed_symbol_and_index :
