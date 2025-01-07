/** @file flintwrap.cc
 *
 *   Contains methods to call flint polynomial interface routines from the rest of FORM.
 */

extern "C" {
#include "form3.h"
}

#include "flintinterface.h"

/*
	#[ flint_factorize_argument :
*/
int flint_factorize_argument(PHEAD WORD *argin, WORD *argout) {
	flint::factorize(BHEAD argin, argout, true, true);
	return 0;
}
/*
	#] flint_factorize_argument :
	#[ flint_factorize_dollar :
*/
WORD* flint_factorize_dollar(PHEAD WORD *argin) {
	return flint::factorize(BHEAD argin, NULL, false, false);
}
/*
	#] flint_factorize_dollar :
	#[ flint_gcd :
*/
WORD* flint_gcd(PHEAD WORD *a, WORD *b, const WORD must_fit) {
	// Extract expressions
	vector<WORD *> e;
	e.push_back(a);
	e.push_back(b);
	const map<unsigned,unsigned> var_map = flint::get_variables(e, false, false);

	if ( var_map.size() > 1 ) {
		return flint::gcd_mpoly(BHEAD a, b, must_fit, var_map);
	}
	else {
		return flint::gcd_poly(BHEAD a, b, must_fit, var_map);
	}
}
/*
	#] flint_gcd :
	#[ flint_ratfun_add :
*/
WORD* flint_ratfun_add(PHEAD WORD *t1, WORD *t2) {

	if ( AR.PolyFunExp == 1 ) {
		MesPrint("flint_ratfun_add: PolyFunExp unimplemented.");
		Terminate(-1);
	}

	WORD *oldworkpointer = AT.WorkPointer;

	// Extract expressions: the num and den of both prf
	vector<WORD *> e;
	for (WORD *t=t1+FUNHEAD; t<t1+t1[1];) {
		e.push_back(t);
		NEXTARG(t);
	}
	for (WORD *t=t2+FUNHEAD; t<t2+t2[1];) {
		e.push_back(t);
		NEXTARG(t);
	}
	const map<unsigned,unsigned> var_map = flint::get_variables(e, true, true);

	if ( var_map.size() > 1 ) {
		flint::ratfun_add_mpoly(BHEAD t1, t2, oldworkpointer, var_map);
	}
	else {
		flint::ratfun_add_poly(BHEAD t1, t2, oldworkpointer, var_map);
	}

	return oldworkpointer;
}
/*
	#] flint_ratfun_add :
	#[ flint_ratfun_normalize :
*/
int flint_ratfun_normalize(PHEAD WORD *term) {

	// The length of the coefficient
	const int ncoeff = (term + *term)[-1];
	// The end of the term data, before the coefficient:
	const WORD *tstop = term + *term - ABS(ncoeff);

	// Search the term for multiple PolyFun or one dirty one.
	int num_polyratfun = 0;
	for (WORD *t = term+1; t < tstop; t += t[1]) {
		if (*t == AR.PolyFun) {
			// Found one!
			num_polyratfun++;
			if ((t[2] & MUSTCLEANPRF) != 0) {
				// Dirty, increment again to force normalisation for single PolyFun
				num_polyratfun++;
			}
			if (num_polyratfun > 1) {
				// We're not counting all occurrences, just determinining if there
				// is something to be done.
				break;
			}
		}
	}
	if (num_polyratfun <= 1) {
		// There is nothing to do, return early
		return 0;
	}

	// When there are polyratfun with only one argument: rename them temporarily to TMPPOLYFUN.
	for (WORD *t = term+1; t < tstop; t += t[1]) {
		if (*t == AR.PolyFun && (t[1] == FUNHEAD+t[FUNHEAD] || t[1] == FUNHEAD+2 ) ) {
			*t = TMPPOLYFUN;
		}
	}


	// Extract all variables in the polyfuns
	// Collect pointers to each relevant argument
	vector<WORD *> e;
	for (WORD *t=term+1; t<tstop; t+=t[1]) {
		if (*t == AR.PolyFun) {
			for (WORD *t2 = t+FUNHEAD; t2<t+t[1];) {
				e.push_back(t2);
				NEXTARG(t2);
			}
		}
	}
	const map<unsigned,unsigned> var_map = flint::get_variables(e, true, true);

	if ( var_map.size() > 1 ) {
		flint::ratfun_normalize_mpoly(BHEAD term, var_map);
	}
	else {
		flint::ratfun_normalize_poly(BHEAD term, var_map);
	}


	// Undo renaming of single-argument PolyFun
	const WORD *new_tstop = term + *term - ABS((term + *term)[-1]);
	for (WORD *t=term+1; t<new_tstop; t+=t[1]) {
		if (*t == TMPPOLYFUN ) *t = AR.PolyFun;
	}

	return 0;
}
/*
	#] flint_ratfun_normalize :
*/
