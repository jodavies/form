* Coverage-oriented tests for normal paths in sources/argument.c.

*--#[ coverage_argument_implode_explode :
CFunction g,h;
Symbol x,y;

* Exercise the compact representations used for a function and a general
* expression after a run of zero arguments. Imploding and exploding again
* must preserve the FORM-level argument list.
Local imploded = g(0,0,h(x)) + h(0,0,x+y);
ArgImplode g,h;
ArgExplode g,h;
Print;
.end
assert succeeded?
assert result("imploded") =~ expr("g(0,0,h(x)) + h(0,0,y+x)")
*--#] coverage_argument_implode_explode :

*--#[ coverage_factarg_content :
CFunction f,g;
Symbol x,y,z;
Vector p,q;
Index mu;

* Common symbol, vector, dot-product, function and inverse factors take the
* normal content-extraction paths before the remaining argument is factored.
Local symbols = f(x*(x+y) + x*(x+2*y));
Local vector = f(p(mu)*x + p(mu)*y);
Local dotproduct = f(p.q*x + p.q*y);
Local inversesymbol = f(x^-1*y + x^-1*z);
Local inversedotproduct = f(p.q^-1*x + p.q^-1*y);
Local functionfactor = f(g(x)*y + g(x)*z);
.sort
FactArg f;
Print;
.end
assert succeeded?
assert result("symbols") =~ expr("f(x,3*y + 2*x)")
assert result("vector") =~ expr("f(y + x,p(mu))")
assert result("dotproduct") =~ expr("f(y + x,p.q)")
assert result("inversesymbol") =~ expr("f(x^-1,z + y)")
assert result("inversedotproduct") =~ expr("f(y + x,p.q^-1)")
assert result("functionfactor") =~ expr("f(z + y,g(x))")
*--#] coverage_factarg_content :

*--#[ coverage_factarg_modulus :
CFunction f,g;
Symbol x,y;

* FactArg under a modulus uses MakeMod rather than the integer GCD/LCM path.
Modulus 17;
Local modularf = f(3*x + 6*y);
Local modularg = g(3*x + 6*y);
.sort
FactArg f,g;
Print;
.end
assert succeeded?
assert result("modularf") =~ expr("f(3,2*y + x)")
assert result("modularg") =~ expr("g(3,2*y + x)")
*--#] coverage_factarg_modulus :

*--#[ coverage_argument_term_environment :
Symbol x,y;

* Term environments sort and reprocess each term independently.
Local termwise = x + y;
Term;
  Multiply 2;
EndTerm;
Print;
.end
assert succeeded?
assert result("termwise") =~ expr("2*y + 2*x")
*--#] coverage_argument_term_environment :

*--#[ coverage_argument_normalize_scales :
CFunction f;
Symbol x,y;

* Exercise Normalize's non-default scale forms.  Positive, negative and
* zero scale all take distinct normal (non-error) execution paths.
Local normalizedpositive = f(2*x + 3*y);
Local normalizednegative = f(2*x + 3*y);
Local normalizedzero = f(2*x + 3*y);
InExpression normalizedpositive;
	Normalize,^2,f;
EndInExpression;
InExpression normalizednegative;
	Normalize,^-1,f;
EndInExpression;
InExpression normalizedzero;
	Normalize,(0),f;
EndInExpression;
Print;
.end
assert succeeded?
assert result("normalizedpositive") =~ expr("9*f(y+2/3*x)")
assert result("normalizednegative") =~ expr("1/3*f(y+2/3*x)")
assert result("normalizedzero") =~ expr("f(y+2/3*x)")
*--#] coverage_argument_normalize_scales :

*--#[ coverage_argument_makeinteger :
CFunction f;
Symbol x,y;

* MakeInteger normalizes both numerator GCDs and denominator LCMs in a
* multi-term argument.
Local integralargument = f(2/3*x + 4/5*y);
MakeInteger f;
Print;
.end
assert succeeded?
assert result("integralargument") =~ expr("2/15*f(6*y+5*x)")
*--#] coverage_argument_makeinteger :

*--#[ coverage_factarg_tensor_and_delta :
CFunction f;
Tensor T;
Symbol x,y;
Index mu,nu;

* Factor extraction has separate representations for tensor functions and
* Kronecker deltas.
Local tensorfactor = f(T(mu,nu)*x + T(mu,nu)*y);
Local deltafactor = f(d_(mu,nu)*x + d_(mu,nu)*y);
FactArg f;
Print;
.end
assert succeeded?
assert result("tensorfactor") =~ expr("f(y+x,T(mu,nu))")
assert result("deltafactor") =~ expr("f(y+x,d_(mu,nu))")
*--#] coverage_factarg_tensor_and_delta :

*--#[ coverage_oldfactarg_tensor :
CFunction f;
Tensor T;
Symbol x,y;
Index mu,nu;

* OldFactArg remains a supported compatibility mode and has a separate
* normal implementation for tensor-function factors.
On OldFactArg;
Local oldtensorfactor = f(T(mu,nu)*x + T(mu,nu)*y);
FactArg f;
Print;
.end
assert succeeded?
assert result("oldtensorfactor") =~ expr("f(T(mu,nu),y+x)")
*--#] coverage_oldfactarg_tensor :

*--#[ coverage_oldfactarg_factor_kinds :
CFunction f,g;
Symbol x,y,z;
Vector p,q;
Index mu;

* The legacy factorizer has independent normal paths for common symbols,
* inverse symbols, dot products, inverse dot products, vectors and ordinary
* functions.
On OldFactArg;
Local oldsymbols = f(x*(x+y) + x*(x+2*y));
Local oldinversesymbol = f(x^-1*y + x^-1*z);
Local olddotproduct = f(p.q*x + p.q*y);
Local oldinversedotproduct = f(p.q^-1*x + p.q^-1*y);
Local oldvector = f(p(mu)*x + p(mu)*y);
Local oldfunction = f(g(x)*y + g(x)*z);
FactArg f;
Print;
.end
assert succeeded?
assert result("oldsymbols") =~ expr("f(x,3*y+2*x)")
assert result("oldinversesymbol") =~ expr("f(x^-1,z+y)")
assert result("olddotproduct") =~ expr("f(p.q,y+x)")
assert result("oldinversedotproduct") =~ expr("f(p.q^-1,y+x)")
assert result("oldvector") =~ expr("f(p(mu),y+x)")
assert result("oldfunction") =~ expr("f(g(x),z+y)")
*--#] coverage_oldfactarg_factor_kinds :

*--#[ coverage_argument_dollar_selected_normalize :
CFunction f;
Symbol x,y,a;

* A dollar in the argument selector is evaluated at execution time before
* normalization chooses the reference term.
* TODO this is not documented
#$reference = x;
Local dollarnormalize = f(2*x + 3*y);
Normalize,($reference),f;
Print;
.end
assert succeeded?
assert result("dollarnormalize") =~ expr("2*f(3/2*y + x)")
*--#] coverage_argument_dollar_selected_normalize :

*--#[ coverage_argument_makeinteger_long_coefficients :
CFunction f;
Symbol x,y;

* Multi-word numerator GCDs and denominator LCMs exercise the long-number
* branches of MakeInteger's argument normalizer.
Local longgcd = f(1180591620717411303424/3*x + 1180591620717411303424/5*y);
Local longlcm = f(x/1180591620717411303424 + y/3);
MakeInteger f;
Print;
.end
assert succeeded?
assert result("longgcd") =~ expr("1180591620717411303424/15*f(3*y + 5*x)")
assert result("longlcm") =~ expr("1/3541774862152233910272*f(1180591620717411303424*y + 3*x)")
*--#] coverage_argument_makeinteger_long_coefficients :

*--#[ coverage_factarg_option_and_short_argument_paths :
CFunction fzero,fplus,fminus,fvector;
Vector p;

* Coefficient options and negative vectors use distinct short-argument
* paths before the general factorization code.
Local factzero = fzero(-2);
Local factplus = fplus(-2);
Local factminus = fminus(-2);
Local factvector = fvector(-p);
FactArg,(0),fzero;
FactArg,(1),fplus;
FactArg,(-1),fminus;
FactArg fvector;
Print;
.end
assert succeeded?
assert result("factzero") =~ expr("fzero(1)")
assert result("factplus") =~ expr("fplus(-1,2)")
assert result("factminus") =~ expr("fminus(-2)")
assert result("factvector") =~ expr("fvector(-1,p)")
*--#] coverage_factarg_option_and_short_argument_paths :

*--#[ coverage_factarg_option_multiterm_coefficients :
CFunction fzero,fplus,fminus;
Symbol x,y;

* The three coefficient modes also have normal multi-term paths, including
* sign extraction and rational coefficient normalization.
Local multizero = fzero(-2*x - 4*y);
Local multiplus = fplus(-2*x - 4*y);
Local multiminus = fminus(-2*x - 4*y);
FactArg,(0),fzero;
FactArg,(1),fplus;
FactArg,(-1),fminus;
Print;
.end
assert succeeded?
assert result("multizero") =~ expr("fzero(2*y+x)")
assert result("multiplus") =~ expr("fplus(2*y+x,-1,2)")
assert result("multiminus") =~ expr("fminus(2*y+x,-2)")
*--#] coverage_factarg_option_multiterm_coefficients :

*--#[ coverage_argument_short_normalize :
CFunction fpos,fneg,fvec;
Vector p;

* Short numeric and negative-vector arguments take a different Normalize
* route from the general polynomial arguments.
Local directpositive = fpos(2);
Local directnegative = fneg(-2);
Local directvector = fvec(-p);
Normalize,^2,fpos;
Normalize,^-1,fneg;
Normalize,^2,fvec;
Print;
.end
assert succeeded?
assert result("directpositive") =~ expr("4*fpos(1)")
assert result("directnegative") =~ expr("-1/2*fneg(1)")
assert result("directvector") =~ expr("fvec(p)")
*--#] coverage_argument_short_normalize :

*--#[ coverage_factarg_high_power_dotproduct :
CFunction f;
Symbol x,y;
Vector p,q;

* A dot product with power greater than one passes through FactArg's
* non-simple-factor path before factorization.
Local dotpower = f(p.q^2*x + p.q^2*y);
FactArg f;
Print;
.end
assert succeeded?
assert result("dotpower") =~ expr("f(y+x,p.q,p.q)")
*--#] coverage_factarg_high_power_dotproduct :

*--#[ coverage_oldfactarg_short_and_one_term :
CFunction fzero,fplus,fminus,fvector;
Symbol x;
Vector p;

* The legacy path handles a negative vector and one-term general arguments
* separately from both short numerical arguments and multi-term arguments.
On OldFactArg;
Local oldzeroone = fzero(2*x);
Local oldplusone = fplus(-2*x);
Local oldminusone = fminus(-2*x);
Local oldvector = fvector(-p);
FactArg,(0),fzero;
FactArg,(1),fplus;
FactArg,(-1),fminus;
FactArg fvector;
Print;
.end
assert succeeded?
assert result("oldzeroone") =~ expr("fzero(x)")
assert result("oldplusone") =~ expr("fplus(x,-1,2)")
assert result("oldminusone") =~ expr("fminus(x,-2)")
assert result("oldvector") =~ expr("fvector(p,-1,1)")
*--#] coverage_oldfactarg_short_and_one_term :

*--#[ coverage_argtoextrasymbol_fast_and_general :
CFunction f,g;
Symbol x,y;
Vector p;
Index mu;

* ArgToExtraSymbol converts both compact arguments and full expressions;
* ToNumber uses the corresponding numeric replacement path.
Local extrasymbols = f(1) + f(x) + f(-p) + f(mu) + f(g) + f(x+y);
Local extranumbers = f(1) + f(x) + f(-p) + f(mu) + f(g) + f(x+y);
InExpression extrasymbols;
	ArgToExtraSymbol f;
EndInExpression;
InExpression extranumbers;
  ArgToExtraSymbol,ToNumber,f;
EndInExpression;
Print;
.sort
#write "%X"
.end
assert succeeded?
assert result("extrasymbols") =~ expr("f(Z6_) + f(Z5_) + f(Z4_) + f(Z3_) + f(Z2_) + f(Z1_)")
assert result("extranumbers") =~ expr("f(1) + f(2) + f(3) + f(4) + f(5) + f(6)")
assert stdout =~ exact_pattern(<<'EOF')
    Z1_=1;
    Z2_=x;
    Z3_= - p;
    Z4_=mu;
    Z5_=g;
    Z6_=y + x;
EOF
*--#] coverage_argtoextrasymbol_fast_and_general :

*--#[ coverage_argument_reprocess_short_and_sum :
CFunction f;
Symbol x,y;

* Argument reprocesses both a transformed short argument and a transformed
* multi-term argument through its normal sorting path.
Local argumentshort = f(x);
Local argumentsum = f(x+y);
Argument f;
  Identify x = 2;
EndArgument;
Print;
.end
assert succeeded?
assert result("argumentshort") =~ expr("f(2)")
assert result("argumentsum") =~ expr("f(2+y)")
*--#] coverage_argument_reprocess_short_and_sum :

*--#[ coverage_argument_normalize_unity :
CFunction fplus,fminus;
Symbol x;

* Unit-coefficient positive and negative terms use Normalize's no-op
* short-circuit cases.
Local positiveunity = fplus(x);
Local negativeunity = fminus(-x);
Normalize fplus;
Normalize,^-1,fminus;
Print;
.end
assert succeeded?
assert result("positiveunity") =~ expr("fplus(x)")
assert result("negativeunity") =~ expr("-fminus(x)")
*--#] coverage_argument_normalize_unity :

*--#[ coverage_factarg_option_negative_vector :
CFunction f;
Vector p;

* TODO is this correct?
* With an explicit coefficient mode, a negative vector follows the
* non-numeric short-argument branch of FactArg.
Local optionvector = f(-p);
FactArg,(1),f;
Print;
.end
assert succeeded?
assert result("optionvector") =~ expr("f(-p,1)")
*--#] coverage_factarg_option_negative_vector :

*--#[ coverage_makeinteger_zero_mode :
CFunction f;
Symbol x,y;

* The (0) selector mode keeps MakeInteger's argument normalization while
* omitting the extracted overall coefficient.
Local makeintegerzero = f(2/3*x + 4/5*y);
MakeInteger,(0),f;
Print;
.end
assert succeeded?
assert result("makeintegerzero") =~ expr("f(6*y+5*x)")
*--#] coverage_makeinteger_zero_mode :

*--#[ coverage_factarg_parenthesized_one_term :
CFunction f;
Symbol x;

* TODO this and the following "parenthesized" tests don't appear to do anything
* meaningful. The code which parses (x), ((x)) etc is common with SplitArg,
* but these are not documented for FactArg. (x) and (-x) are not equivalent to
* (1) and (-1). (x) appears to just pull the sign out, but not for multi-term
* arguments.

* A double-parenthesized argument specification compiles the TYPEFACTARG2
* form, exercising FactArg's one-term handling with a compiled factor.
Local selectedfactor = f(-2*x);
FactArg ((x)) f;
Print;
.end
assert succeeded?
assert result("selectedfactor") =~ expr("f(-1,2*x)")
*--#] coverage_factarg_parenthesized_one_term :

*--#[ coverage_factarg_parenthesized_multiple_arguments :
CFunction f;
Symbol x,y;

* TODO ???
* FactArg rebuilds a function with both untouched and nontrivial arguments
* when it has a double-parenthesized compiled argument specification.
Local selectedunmatched = f(1,x+y);
FactArg ((x)) f;
Print;
.end
assert succeeded?
assert result("selectedunmatched") =~ expr("f(1,y+x,1)")
*--#] coverage_factarg_parenthesized_multiple_arguments :

*--#[ coverage_factarg_parenthesized_one_term_coefficients :
CFunction f;
Symbol x;

* TODO ???
* Parenthesized one-term arguments have special paths for coefficients which
* do not fit FORM's compact one-word rational representation.
Local selectedlarge = f(-1180591620717411303424*x);
Local selectedordinary = f(2*x);
FactArg ((x)) f;
Print;
.end
assert succeeded?
assert result("selectedlarge") =~ expr("f(-1,1180591620717411303424*x)")
assert result("selectedordinary") =~ expr("f(2*x)")
*--#] coverage_factarg_parenthesized_one_term_coefficients :

*--#[ coverage_factarg_parenthesized_one_term_large_coefficients :
CFunction f;
Symbol x;

* TODO ???
* The single-parenthesized compiled argument specification exercises the
* same one-term factor path with large coefficient representations.
* These values have a most-significant limb whose sign bit is set, avoiding
* the compact rational shortcut used for small positive integers.
Local parenthesizedlarge = f(-170141183460469231731687303715884105728/3*x);
Local parenthesizedordinary = f(170141183460469231731687303715884105728*x);
FactArg (x) f;
Print;
.end
assert succeeded?
assert result("parenthesizedlarge") =~ expr("f(-1,170141183460469231731687303715884105728/3*x)")
assert result("parenthesizedordinary") =~ expr("f(170141183460469231731687303715884105728*x)")
*--#] coverage_factarg_parenthesized_one_term_large_coefficients :

*--#[ coverage_factarg_parenthesized_general_factor :
CFunction f;
Symbol x,y;

* TODO ???
* A composite parenthesized factor takes the general numerical-content branch,
* rather than the optimized single-symbol factor path.
Local selectedgeneral = f(-2*x*y-4*x*y^2);
FactArg ((x*y)) f;
Print;
.end
assert succeeded?
assert result("selectedgeneral") =~ expr("f(x*y + 2*x*y^2,-1,2)")
*--#] coverage_factarg_parenthesized_general_factor :

*--#[ coverage_argument_reprocess_fast_arguments :
CFunction f;
Symbol x,y;

* Argument must expand compact numeric and symbolic arguments before running
* the contained statement and compact them again afterwards.
Local reprocessfast = f(2) + f(x) + f(x+y);
Argument f;
  Multiply 2;
EndArgument;
Print;
.end
assert succeeded?
assert result("reprocessfast") =~ expr("f(2*y+2*x) + f(2*x) + f(4)")
*--#] coverage_argument_reprocess_fast_arguments :
