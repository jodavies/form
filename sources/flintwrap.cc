/** @file flintwrap.cc
 *
 *   Contains methods to call flint polynomial routines from the rest of FORM.
 */

extern "C" {
#include "form3.h"
}

#include <flint/fmpq_mpoly.h>

#include <iostream>
#include <vector>
#include <map>
#include <cassert>

// The bits of std that are needed:
using std::cout;
using std::endl;
using std::map;
using std::vector;

// Prototypes
WORD flint_fmpz_get_form(fmpz_t, WORD *);
void flint_fmpz_set_form(fmpz_t, UWORD *, WORD);
map<unsigned,unsigned> flint_get_variables(const vector <WORD *> &, bool, bool);
void flint_mpoly_from_argument(fmpz_mpoly_t, fmpz_mpoly_t, const WORD *, const map<unsigned,unsigned>, const fmpz_mpoly_ctx_t);
ULONG flint_mpoly_to_argument(WORD *, ULONG, const fmpz_mpoly_t, const map<unsigned,unsigned>, const fmpz_mpoly_ctx_t);
int flint_ratfun_normalize(PHEAD WORD *);
void flint_ratfun_read(const WORD *, fmpz_mpoly_t, fmpz_mpoly_t, const map<unsigned,unsigned> &, fmpz_mpoly_ctx_t);

/*
	#[ flint_fmpz_get_form :
*/
// Write FORM's long integer representation of an fmpz at a, and put the number of WORDs at na.
WORD flint_fmpz_get_form(fmpz_t z, WORD *a) {

	WORD na = 0;
	const int sgn = fmpz_sgn(z);
	if ( sgn == -1 ) {
		fmpz_neg(z, z);
	}
	const long nlimbs = fmpz_size(z);

	// This works but is UB?
	//fmpz_get_ui_array(reinterpret_cast<ULONG*>(a), nlimbs, z);

	unsigned long limb_data[nlimbs] = {0};
	fmpz_get_ui_array(limb_data, nlimbs, z);
	for (long i = 0; i < nlimbs; i++) {
		a[2*i] = (WORD)(limb_data[i] & 0xFFFFFFFF);
		na++;
		a[2*i+1] = (WORD)(limb_data[i] >> BITSINWORD);
		if ( a[2*i+1] != 0 || i < (nlimbs-1) ) {
			// The final limb might fit in a single 32bit WORD. Only
			// increment na if the final WORD is non zero.
			na++;
		}
	}

	// And now put the sign in the number of limbs
	if ( sgn == -1 ) {
		na = -na;
	}

	return na;
}
/*
	#] flint_fmpz_get_form :
	#[ flint_fmpz_set_form :
*/
// Create an fmpz directly from FORM's long integer representation. fmpz uses 64bit
// unsigned limbs, but FORM uses 32bit UWORDs on 64bit architectures so we can't use
// fmpz_set_ui_array directly.
void flint_fmpz_set_form(fmpz_t z, UWORD *a, WORD na) {

	if ( na == 0 ) {
		fmpz_zero(z);
		return;
	}

	// Negative na represenents a negative number
	int sgn = 1;
	if ( na < 0 ) {
		sgn = -1;
		na = -na;
	}

	// Remove padding. FORM stores numerators and denominators with equal numbers
	// of limbs but we don't need to do this within the fmpz. It is not necessary
	// to do this really, the fmpz doesn't add zero limbs unnecessarily, but we
	// might be able to avoid creating the limb_data array below.
	while ( a[na-1] == 0 ) {
		na--;
	}

	// If the number fits in fixed-size  fmpz_set functions, we don't need to use
	// additional memory to convert to ULONG. These probably cover most real cases.
	if ( na == 1 ) {
		fmpz_set_ui(z, (ULONG)a[0]);
	}
	else if ( na == 2 ) {
		fmpz_set_ui(z, (((ULONG)a[1])<<BITSINWORD) + (ULONG)a[0]);
	}
	else if ( na == 3 ) {
		fmpz_set_uiui(z, (ULONG)a[2], (((ULONG)a[1])<<BITSINWORD) + (ULONG)a[0]);
	}
	else if ( na == 4 ) {
		fmpz_set_uiui(z, (((ULONG)a[3])<<BITSINWORD) + (ULONG)a[2],
			(((ULONG)a[1])<<BITSINWORD) + (ULONG)a[0]);
	}
	else {
		const LONG nlimbs = (na+1)/2;
		ULONG limb_data[nlimbs] = {0};
		for ( int i = 0; i < nlimbs; i++ ) {
			if ( 2*i+1 <= na-1 ) {
				limb_data[i] = (ULONG)a[2*i] + (((ULONG)a[2*i+1])<<BITSINWORD);
			}
			else {
				limb_data[i] = (ULONG)a[2*i];
			}
		}
		fmpz_set_ui_array(z, limb_data, nlimbs);
	}

	// Finally set the sign.
	if ( sgn == -1 ) {
		fmpz_neg(z, z);
	}

	return;
}
/*
	#] flint_fmpz_set_form :
	#[ flint_get_variables :
*/
// TODO sort vars
map<unsigned,unsigned> flint_get_variables(const vector <WORD *> &es, const bool with_arghead, const bool sort_vars) {
	DUMMYUSE(sort_vars);

	// This does the same job as poly::get_variables. We need to get
	// the list of symbols which appear in the vector of expressions.
	// These correspond to the numerator and denominator of one or
	// more polyratfun.

	int num_vars = 0;
	// To be used if we sort by highest degree, as the polu code does.
	vector<int> degrees;
	map<unsigned,unsigned> var_map;

	// extract all variables
	for (unsigned ei=0; ei < es.size(); ei++) {
		WORD *e = es[ei];

		// fast notation
		if (*e == -SNUMBER) {
		}
		else if (*e == -SYMBOL) {
			if (!var_map.count(e[1])) {
				var_map[e[1]] = num_vars++;
				degrees.push_back(1);
			}
		}
		// JD: Here we need to check for non-symbol/number terms in fast notation.
		else if (*e < 0) {
			MLOCK(ErrorMessageLock);
			MesPrint ((char*)"ERROR: polynomials and polyratfuns must contain symbols only");
			MUNLOCK(ErrorMessageLock);
			Terminate(1);
		}
		else {
			for (int i=with_arghead ? ARGHEAD : 0; with_arghead ? i<e[0] : e[i]!=0; i+=e[i]) {
				if (i+1<i+e[i]-ABS(e[i+e[i]-1]) && e[i+1]!=SYMBOL) {
					MLOCK(ErrorMessageLock);
					MesPrint ((char*)"ERROR: polynomials and polyratfuns must contain symbols only");
					MUNLOCK(ErrorMessageLock);
					Terminate(1);
				}

				for (int j=i+3; j<i+e[i]-ABS(e[i+e[i]-1]); j+=2) {
					if (!var_map.count(e[j])) {
						var_map[e[j]] = num_vars++;
						degrees.push_back(e[j+1]);
					}
					else {
						degrees[var_map[e[j]]] = MaX(degrees[var_map[e[j]]], e[j+1]);
					}
				}
			}
		}
	}

//	cout << "vars: ";
//	for (auto x: var_map) {
//		cout << x.second << " " << x.first << endl;
//	}

	return var_map;
}
/*
	#] flint_get_variables :
	#[ flint_mpoly_from_argument :
*/
// TODO with arghead
void flint_mpoly_from_argument(fmpz_mpoly_t poly, fmpz_mpoly_t denpoly, const WORD *args, const map<unsigned,unsigned> var_map, const fmpz_mpoly_ctx_t ctx) {

	// First check for "fast notation" arguments:
	if ( *args == -SNUMBER ) {
		fmpz_mpoly_set_si(poly, *(args+1), ctx);
		fmpz_mpoly_set_si(denpoly, 1, ctx);
		return;
	}

	if ( *args == -SYMBOL ) {
		// A "fast notation" SYMBOL has a power and coefficient of 1:
		fmpz_t coeff;
		fmpz_init(coeff);
		fmpz_set_si(coeff, (long)1);

		unsigned long exponents[var_map.size()] = {0};
		exponents[var_map.at(*(args+1))] = 1;

		fmpz_mpoly_set_coeff_fmpz_ui(poly, coeff, exponents, ctx);
		fmpz_mpoly_set_si(denpoly, 1, ctx);
		fmpz_clear(coeff);
		return;
	}


	// Now we can iterate through the terms of the argument
	const WORD* arg_stop = args+args[0];
	args += ARGHEAD;


	// Search for numerical or symbol denominators to create "denpoly".
	fmpz_t den_coeff;
	fmpz_init(den_coeff);
	fmpz_set_si(den_coeff, (long)1);
	unsigned long neg_exponents[var_map.size()] = {0};
	for (const WORD* term = args; term < arg_stop; term += term[0]) {

		const WORD* term_stop = term+term[0];
		const WORD  coeff_size = (term_stop)[-1];
		const WORD* symbol_stop = term_stop - ABS(coeff_size);
		const WORD* t = term;

		t++;
		if (t == symbol_stop) {
			// Just a number, no symbols
		}
		else {
			if (*t != SYMBOL) {
				// get_variables already determined that this won't happen
				cout << "    error, found a non-symbol: " << *t << endl;
			}
			t++; // this entry is SYMBOL
			t++; // this entry just has the size of the symbol array, but we can use symbol_stop

			for (const WORD* s = t; s < symbol_stop; s+=2) {
				if ( *(s+1) < 0 ) {
					neg_exponents[var_map.at(*s)] = MaX(neg_exponents[var_map.at(*s)], (unsigned long)(-(*(s+1))) );
				}
			}
		}

		// Now check for a denominator in the coefficient:
		if (*(symbol_stop+ABS(coeff_size/2)) != 1) {
			fmpz_t tmp;
			fmpz_init(tmp);
			flint_fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// Record the LCM of the coefficient denominators:
			fmpz_lcm(den_coeff, den_coeff, tmp);
			fmpz_clear(tmp);
		}
	}
	// Assemble denpoly.
	fmpz_mpoly_set_coeff_fmpz_ui(denpoly, den_coeff, neg_exponents, ctx);


	for (const WORD* term = args; term < arg_stop; term += term[0]) {

		const WORD* term_stop = term+term[0];
		const WORD  coeff_size = (term_stop)[-1];
		const WORD* symbol_stop = term_stop - ABS(coeff_size);
		const WORD* t = term;

		// Create fmpz_t coeff and exponent array from term data:
		fmpz_t coeff;
		fmpz_init(coeff);
		unsigned long exponents[var_map.size()] = {0};

		t++; // skip over the total size entry
		if (t == symbol_stop) {
			// Just a number, no symbols
		}
		else {
			if (*t != SYMBOL) {
				// get_variables already determined that this won't happen
				cout << "    error, found a non-symbol: " << *t << endl;
			}
			t++; // this entry is SYMBOL
			t++; // this entry just has the size of the symbol array, but we can use symbol_stop
			for (const WORD* s = t; s < symbol_stop; s+=2) {
				exponents[var_map.at(*s)] = *(s+1);
			}
		}
		// Now read the coefficient
		flint_fmpz_set_form(coeff, (UWORD*)symbol_stop, coeff_size/2);

		// Multiply by denominator LCM
		fmpz_mul(coeff, coeff, den_coeff);

		// Shift by neg_exponents
		for (unsigned i = 0; i < var_map.size(); i++) {
			exponents[i] += neg_exponents[i];
		}

		// Read the denominator if there is one, and divide it out of the coefficient
		if (*(symbol_stop+ABS(coeff_size/2)) != 1) {
			fmpz_t tmp;
			fmpz_init(tmp);
			flint_fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// By construction, this is an exact division
			fmpz_divexact(coeff, coeff, tmp);
			fmpz_clear(tmp);
		}

		// Add the term to the poly
		fmpz_mpoly_set_coeff_fmpz_ui(poly, coeff, exponents, ctx);
		fmpz_clear(coeff);
	}

	fmpz_clear(den_coeff);
}
/*
	#] flint_mpoly_from_argument :
	#[ flint_mpoly_to_argument :
*/
// TODO with arghead
ULONG flint_mpoly_to_argument(WORD *out, ULONG prev_size, const fmpz_mpoly_t poly, const map<unsigned,unsigned> var_map, const fmpz_mpoly_ctx_t ctx) {

	fmpz_t coeff;
	fmpz_init(coeff);
	LONG exponents[var_map.size()];

	const LONG n_terms = fmpz_mpoly_length(poly, ctx);

	if ( n_terms == 0 ) {
// TODO arghead
		*out++ = -SNUMBER;
		*out++ = 0;
		return 2;
	}

	// TODO fast notation

	WORD* arg_size = out++; // total arg size
	WORD* arg_flag = out++;
//	*arg_flag = var_map.size() > 1 ? DIRTYFLAG : 0;
	*arg_flag = 0; // TODO dirty arg ?

	for (LONG i = 0; i < n_terms; i++) {
		WORD* term_size = out++;

		fmpz_mpoly_get_term_exp_si(exponents, poly, i, ctx);
		fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

		bool has_symbols = 0;
		for ( unsigned i = 0; i < var_map.size(); i++ ) {
			if ( exponents[i] != 0 ) { has_symbols = 1; }
		}
		if ( has_symbols ) {
			*out++ = SYMBOL;
			WORD* symbol_size = out++;
			*symbol_size = 2;

			for ( unsigned i = 0; i < var_map.size(); i++ ) {
				if ( exponents[i] != 0 ) {
					// TODO this seems bad
					for (auto x: var_map) {
						if ( x.second == i ) {
							*out++ = x.first;
						}
					}
					*out++ = exponents[i];
					*symbol_size += 2;
				}
			}
		}

		int coeff_size = flint_fmpz_get_form(coeff, out);
		out += ABS(coeff_size);
		*out++ = 1; // the denominator
		for (int i = 0; i < ABS(coeff_size)-1; i++) {
			*out++ = 0;
		}
		*out = 2*ABS(coeff_size) + 1; // the size of the coefficient
		if ( coeff_size < 0 ) { *out = -(*out); }
		out++;

		*term_size = out - term_size;

		if ( sizeof(WORD)*(prev_size + out - arg_size) >= (size_t)AM.MaxTer ) {
			MLOCK(ErrorMessageLock);
			MesPrint("flint_mpoly_to_argument: output exceeds MaxTermSize");
			MUNLOCK(ErrorMessageLock);
			Terminate(-1);
		}

	}

	*arg_size = out - arg_size;

	return *arg_size;
}
/*
	#] flint_mpoly_to_argument :
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
	const map<unsigned,unsigned> var_map = flint_get_variables(e, true, true);


	// Now we know how many variables appear in all polyratfuns. Prepare flint context.
	fmpz_mpoly_ctx_t ctx;
	fmpz_mpoly_ctx_init(ctx, var_map.size(), ORD_LEX);

	fmpz_mpoly_t gcd, num1, den1, num2, den2;
	fmpz_mpoly_init(gcd, ctx);
	fmpz_mpoly_init(num1, ctx);
	fmpz_mpoly_init(den1, ctx);
	fmpz_mpoly_init(num2, ctx);
	fmpz_mpoly_init(den2, ctx);

	flint_ratfun_read(t1, num1, den1, var_map, ctx);
	flint_ratfun_read(t2, num2, den2, var_map, ctx);

	if ( fmpz_mpoly_cmp(den1, den2, ctx) != 0 ) {
		fmpz_mpoly_gcd_cofactors(gcd, den1, den2, den1, den2, ctx);

		fmpz_mpoly_mul(num1, num1, den2, ctx);
		fmpz_mpoly_mul(num2, num2, den1, ctx);

		fmpz_mpoly_add(num1, num1, num2, ctx);
		fmpz_mpoly_mul(den1, den1, den2, ctx);
		fmpz_mpoly_mul(den1, den1, gcd,  ctx);
	}
	else {
		fmpz_mpoly_add(num1, num1, num2, ctx);
	}

	// Finally divide out any common factors between the resulting num1, den1:
	fmpz_mpoly_gcd_cofactors(gcd, num1, den1, num1, den1, ctx);

	// Fix sign: the leading denominator term should have a positive coeff.
//	fmpz_t leading_coeff;
//	fmpz_init(leading_coeff);
//	fmpz_mpoly_get_term_coeff_fmpz(leading_coeff, den1, 0, ctx);
//	if ( fmpz_sgn(leading_coeff) == -1 ) {
//		fmpz_mpoly_neg(num1, num1, ctx);
//		fmpz_mpoly_neg(den1, den1, ctx);
//	}
//	fmpz_clear(leading_coeff);

	// Result in FORM notation:
	WORD* out = oldworkpointer;
	*out++ = AR.PolyFun;
	WORD* args_size = out++;
	WORD* args_flag = out++;
//	*args_flag &= ~MUSTCLEANPRF;
	*args_flag = 0; // clean prf
	FILLFUN3(out); // TODO ?

	out += flint_mpoly_to_argument(out, out-args_size, num1, var_map, ctx);
	out += flint_mpoly_to_argument(out, out-args_size, den1, var_map, ctx);

	*args_size = out - args_size + 1; // The +1 is to include the function ID
	AT.WorkPointer = out;


	fmpz_mpoly_clear(num1, ctx);
	fmpz_mpoly_clear(den1, ctx);
	fmpz_mpoly_clear(num2, ctx);
	fmpz_mpoly_clear(den2, ctx);
	fmpz_mpoly_clear(gcd, ctx);
	fmpz_mpoly_ctx_clear(ctx);

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
	const map<unsigned,unsigned> var_map = flint_get_variables(e, true, true);


	// Now we know how many variables appear in all polyratfuns. Prepare flint context.
	fmpz_mpoly_ctx_t ctx;
	fmpz_mpoly_ctx_init(ctx, var_map.size(), ORD_LEX);

	fmpz_mpoly_t num1, den1;
	// TODO might want to use init2, if we have some idea of the number of terms which appear
	// some lower bound can be determined in get_variables?
	fmpz_mpoly_init(num1, ctx);
	fmpz_mpoly_init(den1, ctx);
	// Start with "trivial" polynomials, and multiply in the num and den of the prf which appear.
	fmpz_t tmpNum, tmpDen;
	fmpz_init(tmpNum);
	fmpz_init(tmpDen);
	flint_fmpz_set_form(tmpNum, (UWORD*)tstop, ncoeff/2);
	flint_fmpz_set_form(tmpDen, (UWORD*)tstop+ABS(ncoeff/2), ABS(ncoeff/2));
	fmpz_mpoly_set_fmpz(num1, tmpNum, ctx);
	fmpz_mpoly_set_fmpz(den1, tmpDen, ctx);
	fmpz_clear(tmpNum);
	fmpz_clear(tmpDen);

	// Loop over the occurrences of PolyFun in the term, and multiply in to num1, den1.
	// s tracks where we are writing the non-PolyFun term data. The final PolyFun will
	// go at the end.
	WORD* term_size = term;
	WORD* s = term + 1;
	for (WORD *t=term+1; t<tstop;) {
		if (*t == AR.PolyFun) {

			fmpz_mpoly_t num2, den2;
			fmpz_mpoly_init(num2, ctx);
			fmpz_mpoly_init(den2, ctx);
			flint_ratfun_read(t, num2, den2, var_map, ctx);

			if ((t[2] & MUSTCLEANPRF) != 0) { // first normalize TODO do this inside read?
				fmpz_mpoly_t gcd;
				fmpz_mpoly_init(gcd, ctx);
				if (!fmpz_mpoly_gcd_cofactors(gcd, num2, den2, num2, den2, ctx)) {
					cout << "flint_ratfun_normalize: error in fmpz_mpoly_gcd_cofactors" << endl;
					Terminate(-1);
				}
				fmpz_mpoly_clear(gcd, ctx);
			}

			// get gcd of num1,den2 and num2,den1 and then assemble
			fmpz_mpoly_t gcd;
			fmpz_mpoly_init(gcd, ctx);
			fmpz_mpoly_gcd_cofactors(gcd, num1, den2, num1, den2, ctx);
			fmpz_mpoly_gcd_cofactors(gcd, num2, den1, num2, den1, ctx);
			fmpz_mpoly_clear(gcd, ctx);

			fmpz_mpoly_mul(num1, num1, num2, ctx);
			fmpz_mpoly_mul(den1, den1, den2, ctx);
			fmpz_mpoly_clear(num2, ctx);
			fmpz_mpoly_clear(den2, ctx);

			t += t[1];
		}

		else {
			// Not a PolyFun, just copy or skip over
			int i = t[1];
			if (s != t) { NCOPY(s,t,i); }
			else { t += i; s += i; }
		}
	}

	// Fix sign: leading term of den should be positive. Maybe not necessary,
	// does gcd_cofactors already arrange for this somehow?
	// Fix sign
//	fmpz_t leading_coeff;
//	fmpz_init(leading_coeff);
//	fmpz_mpoly_get_term_coeff_fmpz(leading_coeff, den1, 0, ctx);
//	if ( fmpz_sgn(leading_coeff) == -1 ) {
//		fmpz_mpoly_neg(num1, num1, ctx);
//		fmpz_mpoly_neg(den1, den1, ctx);
//	}
//	fmpz_clear(leading_coeff);
	

	// Result in FORM notation:
	WORD* out = s;
	*out++ = AR.PolyFun;
	WORD* args_size = out++;
	WORD* args_flag = out++;
	*args_flag &= ~MUSTCLEANPRF;

	out += flint_mpoly_to_argument(out, out-args_size, num1, var_map, ctx);
	out += flint_mpoly_to_argument(out, out-args_size, den1, var_map, ctx);

	*args_size = out - args_size + 1; // The +1 is to include the function ID

	// +3 for the coefficient of 1/1, which is added after the check
	if ( sizeof(WORD)*(*args_size+3) >= (size_t)AM.MaxTer ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint_ratfun_normalize: output exceeds MaxTermSize");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	*out++ = 1;
	*out++ = 1;
	*out++ = 3; // the term's coefficient is now 1/1

	*term_size = out - term_size;


	fmpz_mpoly_clear(num1, ctx);
	fmpz_mpoly_clear(den1, ctx);
	fmpz_mpoly_ctx_clear(ctx);


	// Undo renaming of single-argument PolyFun
	const WORD *new_tstop = term + *term - ABS((term + *term)[-1]);
	for (WORD *t=term+1; t<new_tstop; t+=t[1]) {
		if (*t == TMPPOLYFUN ) *t = AR.PolyFun;
	}

	return 0;
}
/*
	#] flint_ratfun_normalize :
	#[ flint_ratfun_read :
*/
void flint_ratfun_read(const WORD *a, fmpz_mpoly_t num, fmpz_mpoly_t den, const map<unsigned,unsigned> &var_map, fmpz_mpoly_ctx_t ctx) {

//	MesPrint("flint_ratfun_read: %a", *(a+1), a);

	// The end of the arguments:
	const WORD* arg_stop = a+a[1];

	// TODO Do we need to normalize? This is done outside of read currently, but maybe it is better here?
//	const bool clean = (a[2] & MUSTCLEANPRF) == 0;

	a += FUNHEAD;
	if (a >= arg_stop) {
		MLOCK(ErrorMessageLock);
		MesPrint ((char*)"ERROR: PolyRatFun cannot have zero arguments");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// Polys to collect the "de of the num" and "den of the den".
	// Input can arrive like this when enabling the PolyRatFun or moving things into it.
	fmpz_mpoly_t den_num;
	fmpz_mpoly_init(den_num, ctx);
	fmpz_mpoly_t den_den;
	fmpz_mpoly_init(den_den, ctx);

	// Read the numerator
	flint_mpoly_from_argument(num, den_num, a, var_map, ctx);

	NEXTARG(a);
	if ( a < arg_stop ) {
		// Read the denominator
		flint_mpoly_from_argument(den, den_den, a, var_map, ctx);
		NEXTARG(a);
	}
	else {
		// The denominator is 1
		cout << "implement this" << endl;
		Terminate(-1);
		// TODO I am not sure how to get here... the poly code seems buggy?
	}

	// Multiply the num by den_den and den by den_num:
	fmpz_mpoly_mul(num, num, den_den, ctx);
	fmpz_mpoly_mul(den, den, den_num, ctx);


	fmpz_mpoly_clear(den_num, ctx);
	fmpz_mpoly_clear(den_den, ctx);

	if (a < arg_stop) {
		MLOCK(ErrorMessageLock);
		MesPrint ((char*)"ERROR: PolyRatFun cannot have more than two arguments");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
}
/*
	#] flint_ratfun_read :
*/
