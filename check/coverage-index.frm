* Coverage-oriented tests for normal paths in sources/index.c.

*--#[ coverage_index_putinside_all :
Symbol x,y;
CFunction wrap;

* With no bracket selector, PutInside wraps each complete term.
Local allinside = 2*x + y;
PutInside wrap;
Print;
.end
assert succeeded?
assert result("allinside") =~ expr("wrap(2*x) + wrap(y)")
*--#] coverage_index_putinside_all :

*--#[ coverage_index_putinside_bracketed :
Symbol x,y,z;
CFunction wrap;

* Select only the x-dependent part of each term for the function argument.
Local selected = 1 + x + 2*x + y^2 + z^3;
PutInside wrap,x;
Print;
.end
assert succeeded?
assert result("selected") =~ expr("wrap(2*x) + wrap(x) + wrap(1) + wrap(1)*z^3 + wrap(1)*y^2")
*--#] coverage_index_putinside_bracketed :

*--#[ coverage_index_bracket_lookup :
Symbol x,y;
CFunction f;

* Build a bracket index, hide the expression, and retrieve a single bracket.
Local stored = f(1) + f(2)*x + f(3)*y;
B+ f;
.sort
Hide;
Local lookup = stored[f(2)];
Print;
.end
assert succeeded?
assert result("lookup") =~ expr("x")
*--#] coverage_index_bracket_lookup :

*--#[ coverage_index_many_and_missing_brackets :
CFunction f;

* More than the initial index capacity exercises index growth. A lookup of
* a bracket that is not present is a normal empty-result case.
Local indexed = <f(1)>+...+<f(80)>;
B+ f;
.sort
Hide;
Local found = indexed[f(60)];
Local missing = indexed[f(81)];
Print;
.end
assert succeeded?
assert result("found") =~ expr("1")
assert result("missing") =~ expr("0")
*--#] coverage_index_many_and_missing_brackets :
