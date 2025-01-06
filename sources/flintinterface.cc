/** @file flintwrap.cc
 *
 *   Contains methods to call flint polynomial routines from the rest of FORM.
 */

extern "C" {
#include "form3.h"
}

#include "flintinterface.h"

/*
	#[ flint::fmpz_get_form :
*/
// Write FORM's long integer representation of an fmpz at a, and put the number of WORDs at na.
// na carries the sign of the integer.
WORD flint::fmpz_get_form(fmpz_t z, WORD *a) {

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
	#] flint::fmpz_get_form :
	#[ flint::fmpz_set_form :
*/
// Create an fmpz directly from FORM's long integer representation. fmpz uses 64bit
// unsigned limbs, but FORM uses 32bit UWORDs on 64bit architectures so we can't use
// fmpz_set_ui_array directly.
void flint::fmpz_set_form(fmpz_t z, UWORD *a, WORD na) {

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
	#] flint::fmpz_set_form :
	#[ flint::gcd_mpoly :
*/
// Return a pointer to a buffer containing the GCD of the 0-terminated term lists at a and b.
// If must_fit_term, this should be a TermMalloc buffer. Otherwise Malloc1 the buffer.
// For multi-variate cases.
WORD* flint::gcd_mpoly(PHEAD const WORD *a, const WORD *b, const WORD must_fit_term, const map<unsigned,unsigned> &var_map) {

	fmpz_mpoly_ctx_t ctx;
	fmpz_mpoly_ctx_init(ctx, var_map.size(), ORD_LEX);

	fmpz_mpoly_t pa, pb, denpa, denpb, gcd;
	fmpz_mpoly_init(pa, ctx);
	fmpz_mpoly_init(pb, ctx);
	fmpz_mpoly_init(denpa, ctx);
	fmpz_mpoly_init(denpb, ctx);
	fmpz_mpoly_init(gcd, ctx);

	const unsigned size_a = flint::mpoly_from_argument(pa, denpa, a, false, var_map, ctx);
	const unsigned size_b = flint::mpoly_from_argument(pb, denpb, b, false, var_map, ctx);
	// denpa, denpb should be 1:
	if ( fmpz_mpoly_is_one(denpa, ctx) != 1 ) {
		cout << "gcd error: denpa != 1";
		Terminate(-1);
	}
	if ( fmpz_mpoly_is_one(denpb, ctx) != 1 ) {
		cout << "gcd error: denpb != 1";
		Terminate(-1);
	}

	fmpz_mpoly_gcd(gcd, pa, pb, ctx);

	// This is freed by the caller
	WORD *res;
	if ( must_fit_term ) {
		res = TermMalloc("flint::gcd_mpoly");
	}
	else {
		// The GCD result can not exceed the size of the smaller of the input polynomials.
		// TODO This will often be an over-allocation, but we avoid determining the exact size.
		res = (WORD *)Malloc1(sizeof(WORD)*(size_a < size_b ? size_a : size_b), "flint::gcd_mpoly");
	}

	flint::mpoly_to_argument(BHEAD res, false, false, 0, gcd, var_map, ctx);

	fmpz_mpoly_clear(pa, ctx);
	fmpz_mpoly_clear(pb, ctx);
	fmpz_mpoly_clear(denpa, ctx);
	fmpz_mpoly_clear(denpb, ctx);
	fmpz_mpoly_clear(gcd, ctx);
	fmpz_mpoly_ctx_clear(ctx);

	return res;
}
/*
	#] flint::gcd_mpoly :
	#[ flint::gcd_poly :
*/
// Return a pointer to a buffer containing the GCD of the 0-terminated term lists at a and b.
// If must_fit_term, this should be a TermMalloc buffer. Otherwise Malloc1 the buffer.
// For uni-variate cases.
WORD* flint::gcd_poly(PHEAD const WORD *a, const WORD *b, const WORD must_fit_term, const map<unsigned,unsigned> &var_map) {

	fmpz_poly_t pa, pb, denpa, denpb, gcd;
	fmpz_poly_init(pa);
	fmpz_poly_init(pb);
	fmpz_poly_init(denpa);
	fmpz_poly_init(denpb);
	fmpz_poly_init(gcd);

	const unsigned size_a = flint::poly_from_argument(pa, denpa, a, false);
	const unsigned size_b = flint::poly_from_argument(pb, denpb, b, false);
	// denpa, denpb should be 1:
	if ( fmpz_poly_is_one(denpa) != 1 ) {
		cout << "gcd error: denpa != 1";
		Terminate(-1);
	}
	if ( fmpz_poly_is_one(denpb) != 1 ) {
		cout << "gcd error: denpb != 1";
		Terminate(-1);
	}

	fmpz_poly_gcd(gcd, pa, pb);

	// This is freed by the caller
	WORD *res;
	if ( must_fit_term ) {
		res = TermMalloc("flint::gcd_poly");
	}
	else {
		// The GCD result can not exceed the size of the smaller of the input polynomials.
		// TODO This will often be an over-allocation, but we avoid determining the exact size.
		res = (WORD *)Malloc1(sizeof(WORD)*(size_a < size_b ? size_a : size_b), "flint::gcd_poly");
	}

	flint::poly_to_argument(BHEAD res, false, false, 0, gcd, var_map);

	fmpz_poly_clear(pa);
	fmpz_poly_clear(pb);
	fmpz_poly_clear(denpa);
	fmpz_poly_clear(denpb);
	fmpz_poly_clear(gcd);

	return res;
}
/*
	#] flint::gcd_poly :
	#[ flint::get_variables :
*/
// Get the list of symbols which appear in the vector of expressions. These are polyratfun
// numerators, denominators or expressions from calls to gcd_ etc. Return this list as a
// map between indices and symbol codes.
// TODO sort vars
map<unsigned,unsigned> flint::get_variables(const vector <WORD *> &es, const bool with_arghead, const bool sort_vars) {
	DUMMYUSE(sort_vars);

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

	return var_map;
}
/*
	#] flint::get_variables :
	#[ flint::mpoly_from_argument :
*/
// Convert a FORM argument (or 0-terminated list of terms: with_arghead == false) to
// a (multi-variate) fmpz_mpoly_t poly. The "denominator" is return in denpoly,
// and contains the overall negative-power numeric and symbolic factor.
unsigned flint::mpoly_from_argument(fmpz_mpoly_t poly, fmpz_mpoly_t denpoly, const WORD *args, const bool with_arghead,
	const map<unsigned,unsigned> &var_map, const fmpz_mpoly_ctx_t ctx) {

	// First check for "fast notation" arguments:
	if ( *args == -SNUMBER ) {
		fmpz_mpoly_set_si(poly, *(args+1), ctx);
		fmpz_mpoly_set_si(denpoly, 1, ctx);
		return 2;
	}

	if ( *args == -SYMBOL ) {
		// A "fast notation" SYMBOL has a power and coefficient of 1:
		unsigned long exponents[var_map.size()] = {0};
		exponents[var_map.at(*(args+1))] = 1;

		fmpz_mpoly_set_coeff_ui_ui(poly, (unsigned long)1, exponents, ctx);
		fmpz_mpoly_set_ui(denpoly, (unsigned long)1, ctx);
		return 2;
	}


	// Now we can iterate through the terms of the argument. If we have
	// an ARGHEAD, we already know where to terminate. Otherwise we'll have
	// to loop until the terminating 0.
	const WORD* arg_stop = with_arghead ? args+args[0] : (WORD*)ULONG_MAX;
	unsigned arg_size = 0;
	if ( with_arghead ) {
		arg_size = args[0];
		args += ARGHEAD;
	}


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
			flint::fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// Record the LCM of the coefficient denominators:
			fmpz_lcm(den_coeff, den_coeff, tmp);
			fmpz_clear(tmp);
		}

		if ( *term_stop == 0 ) {
			// This should only ever happen when with_arghead is false
			if ( with_arghead ) {
				cout << "flint::mpoly_from_argument: arghead error" << endl;
				Terminate(-1);
			}
			// + 1 for the terminating 0
			arg_size = term_stop - args + 1;
			break;
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
			t++; // this entry is SYMBOL
			t++; // this entry just has the size of the symbol array, but we can use symbol_stop
			for (const WORD* s = t; s < symbol_stop; s+=2) {
				exponents[var_map.at(*s)] = *(s+1);
			}
		}
		// Now read the coefficient
		flint::fmpz_set_form(coeff, (UWORD*)symbol_stop, coeff_size/2);

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
			flint::fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// By construction, this is an exact division
			fmpz_divexact(coeff, coeff, tmp);
			fmpz_clear(tmp);
		}

		// Add the term to the poly
		fmpz_mpoly_set_coeff_fmpz_ui(poly, coeff, exponents, ctx);
		fmpz_clear(coeff);

		if ( *term_stop == 0 ) {
			// This should only ever happen when with_arghead is false
			if ( with_arghead ) {
				cout << "flint::mpoly_from_argument: arghead error" << endl;
				Terminate(-1);
			}
			break;
		}

	}

	fmpz_clear(den_coeff);

	return arg_size;
}
/*
	#] flint::mpoly_from_argument :
	#[ flint::mpoly_to_argument :
*/
// Convert a fmpz_mpoly_t to a FORM argument (or 0-terminated list of terms: with_arghead == false).
// If the caller is building an output term, prev_size contains the size of the term so far, to
// check that the output fits if must_fit_term.
ULONG flint::mpoly_to_argument(PHEAD WORD *out, const bool with_arghead, const bool must_fit_term, const ULONG prev_size,
	const fmpz_mpoly_t poly, const map<unsigned,unsigned> &var_map, const fmpz_mpoly_ctx_t ctx) {

	// Check there is at least space for ARGHEAD WORDs (the arghead or a fast-notation number)
	if ( must_fit_term && (sizeof(WORD)*(prev_size + ARGHEAD) > (size_t)AM.MaxTer) ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::mpoly_to_argument: output exceeds MaxTermSize");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	const LONG n_terms = fmpz_mpoly_length(poly, ctx);

	if ( n_terms == 0 ) {
		if ( with_arghead ) {
			*out++ = -SNUMBER;
			*out = 0;
			return 2;
		}
		else {
			*out = 0;
			return 1;
		}
	}

	fmpz_t coeff; // There is no fmpz_mpoly_get_term_coeff_ptr
	fmpz_init(coeff);
	WORD *tmp_coeff = (WORD *)NumberMalloc("flint::mpoly_to_argument");
	LONG exponents[var_map.size()];

	// Create the inverse of var_map, so we don't have to search it for each symbol written
	map<unsigned,unsigned> var_map_inv;
	for (auto x: var_map) {
		var_map_inv[x.second] = x.first;
	}

	// TODO fast notation for single symbol arg

	WORD* arg_size = 0;
	WORD* arg_flag = 0;
	if ( with_arghead ) {
		arg_size = out++; // total arg size
		arg_flag = out++;
		*arg_flag = 0; // clean argument
	}
	else {
		// If we are not writing an ARGHEAD, we still keep track of the
		// size of the sum of terms, but must not write to this location!
		arg_size = out;
	}

	for (LONG i = 0; i < n_terms; i++) {

		fmpz_mpoly_get_term_exp_si(exponents, poly, i, ctx);
		fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);

		int num_symbols = 0;
		for ( unsigned i = 0; i < var_map.size(); i++ ) {
			if ( exponents[i] != 0 ) { num_symbols += 1; }
		}

		// Convert the coefficient, write in temporary space
		const int coeff_size = flint::fmpz_get_form(coeff, tmp_coeff);

		// Now we have the number of symbols and the coeff size, we can determine the output size.
		// Check it fits if necessary: term size, num,den of "coeff_size", +1 for the total coeff size
		unsigned current_size = prev_size + (out - arg_size) + 1 + 2*ABS(coeff_size) + 1;
		if ( num_symbols ) {
			// symbols header, code,power of each symbol:
			current_size += 2 + 2*num_symbols;
		}
		if ( must_fit_term && (sizeof(WORD)*current_size > (size_t)AM.MaxTer) ) {
			MLOCK(ErrorMessageLock);
			MesPrint("flint::mpoly_to_argument: output exceeds MaxTermSize");
			MUNLOCK(ErrorMessageLock);
			Terminate(-1);
		}

		WORD* term_size = out++;
		if ( num_symbols ) {
			*out++ = SYMBOL;
			WORD* symbol_size = out++;
			*symbol_size = 2;

			for ( unsigned i = 0; i < var_map.size(); i++ ) {
				if ( exponents[i] != 0 ) {
					*out++ = var_map_inv[i];
					*out++ = exponents[i];
					*symbol_size += 2;
				}
			}
		}

		// Copy numerator
		for ( int i = 0; i < ABS(coeff_size); i++ ) {
			*out++ = tmp_coeff[i];
		}
		*out++ = 1; // the denominator
		for (int i = 0; i < ABS(coeff_size)-1; i++) {
			*out++ = 0;
		}
		*out = 2*ABS(coeff_size) + 1; // the size of the coefficient
		if ( coeff_size < 0 ) { *out = -(*out); }
		out++;

		*term_size = out - term_size;
	}

	if ( with_arghead ) {
		*arg_size = out - arg_size;
	}
	else {
		// with no arghead, we write a terminating zero
		*out++ = 0;
	}

	fmpz_clear(coeff);
	NumberFree(tmp_coeff, "flint::mpoly_to_argument");

	return out - arg_size;
}
/*
	#] flint::mpoly_to_argument :
	#[ flint::poly_from_argument :
*/
// Convert a FORM argument (or 0-terminated list of terms: with_arghead == false) to
// a (uni-variate) fmpz_poly_t poly. The "denominator" is return in denpoly,
// and contains the overall negative-power numeric and symbolic factor.
unsigned flint::poly_from_argument(fmpz_poly_t poly, fmpz_poly_t denpoly, const WORD *args, const bool with_arghead) {

	// First check for "fast notation" arguments:
	if ( *args == -SNUMBER ) {
		fmpz_poly_set_si(poly, *(args+1));
		fmpz_poly_set_si(denpoly, 1);
		return 2;
	}
	
	if ( *args == -SYMBOL ) {
		// A "fast notation" SYMBOL has a power and coefficient of 1:
		fmpz_poly_set_coeff_si(poly, 1, 1);
		fmpz_poly_set_si(denpoly, 1);
		return 2;
	}

	// Now we can iterate through the terms of the argument. If we have
	// an ARGHEAD, we already know where to terminate. Otherwise we'll have
	// to loop until the terminating 0.
	const WORD* arg_stop = with_arghead ? args+args[0] : (WORD*)ULONG_MAX;
	unsigned arg_size = 0;
	if ( with_arghead ) {
		arg_size = args[0];
		args += ARGHEAD;
	}

	// Search for numerical or symbol denominators to create "denpoly".
	fmpz_t den_coeff;
	fmpz_init(den_coeff);
	fmpz_set_si(den_coeff, 1);
	unsigned long neg_exponent = 0;
	for (const WORD* term = args; term < arg_stop; term += term[0]) {

		const WORD* term_stop = term+term[0];
		const WORD  coeff_size = (term_stop)[-1];
		const WORD* symbol_stop = term_stop - ABS(coeff_size);
		const WORD* t = term;

		t++; // skip over the total size entry
		if (t == symbol_stop) {
			// Just a number, no symbols
		}
		else {
			t++; // this entry is SYMBOL
			t++; // this entry is the size of the symbol array
			t++; // this is the first (and only) symbol code
			if ( *t < 0 ) {
				neg_exponent = MaX(neg_exponent, (unsigned long)(-(*t)) );
			}
			t++;
		}

		// Now check for a denominator in the coefficient:
		if (*(symbol_stop+ABS(coeff_size/2)) != 1) {
			fmpz_t tmp;
			fmpz_init(tmp);
			flint::fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// Record the LCM of the coefficient denominators:
			fmpz_lcm(den_coeff, den_coeff, tmp);
			fmpz_clear(tmp);
		}

		if ( *term_stop == 0 ) {
			// This should only ever happen when with_arghead is false
			if ( with_arghead ) {
				cout << "flint::poly_from_argument: arghead error" << endl;
				Terminate(-1);
			}
			// + 1 for the terminating 0
			arg_size = term_stop - args + 1;
			break;
		}
	}
	// Assemble denpoly.
	fmpz_poly_set_coeff_fmpz(denpoly, neg_exponent, den_coeff);


	for (const WORD* term = args; term < arg_stop; term += term[0]) {

		const WORD* term_stop = term+term[0];
		const WORD  coeff_size = (term_stop)[-1];
		const WORD* symbol_stop = term_stop - ABS(coeff_size);
		const WORD* t = term;

		// Create fmpz_t coeff and exponent array from term data:
		fmpz_t coeff;
		fmpz_init(coeff);
		unsigned long exponent = 0;

		t++; // skip over the total size entry
		if (t == symbol_stop) {
			// Just a number, no symbols
		}
		else {
			t++; // this entry is SYMBOL
			t++; // this entry is the size of the symbol array
			t++; // this is the first (and only) symbol code
			exponent = *t++;
		}

		// Now read the coefficient
		flint::fmpz_set_form(coeff, (UWORD*)symbol_stop, coeff_size/2);

		// Multiply by denominator LCM
		fmpz_mul(coeff, coeff, den_coeff);

		// Shift by neg_exponent
		exponent += neg_exponent;

		// Read the denominator if there is one, and divide it out of the coefficient
		if (*(symbol_stop+ABS(coeff_size/2)) != 1) {
			fmpz_t tmp;
			fmpz_init(tmp);
			flint::fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// By construction, this is an exact division
			fmpz_divexact(coeff, coeff, tmp);
			fmpz_clear(tmp);
		}

		// Add the term to the poly
		fmpz_poly_set_coeff_fmpz(poly, exponent, coeff);
		fmpz_clear(coeff);

		if ( *term_stop == 0 ) {
			// This should only ever happen when with_arghead is false
			if ( with_arghead ) {
				cout << "flint::poly_from_argument: arghead error" << endl;
				Terminate(-1);
			}
			break;
		}
	}

	fmpz_clear(den_coeff);

	return arg_size;
}
/*
	#] flint::poly_from_argument :
	#[ flint::poly_to_argument :
*/
// Convert a fmpz_poly_t to a FORM argument (or 0-terminated list of terms: with_arghead == false).
// If the caller is building an output term, prev_size contains the size of the term so far, to
// check that the output fits if must_fit_term.
ULONG flint::poly_to_argument(PHEAD WORD *out, const bool with_arghead, const bool must_fit_term, const ULONG prev_size,
	const fmpz_poly_t poly, const map<unsigned,unsigned> &var_map) {

	// Check there is at least space for ARGHEAD WORDs (the arghead or a fast-notation number)
	if ( must_fit_term && (sizeof(WORD)*(prev_size + ARGHEAD) > (size_t)AM.MaxTer) ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::poly_to_argument: output exceeds MaxTermSize");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	const LONG n_terms = fmpz_poly_length(poly);

	// The poly is zero
	if ( n_terms == 0 ) {
		if ( with_arghead ) {
			*out++ = -SNUMBER;
			*out = 0;
			return 2;
		}
		else {
			*out = 0;
			return 1;
		}
	}

	// The poly is constant, use fast notation if the coefficient is small enough
	if ( with_arghead && n_terms == 1 ) {
		if ( fmpz_fits_si(fmpz_poly_get_coeff_ptr(poly, 0)) ) {
			long fast_coeff = fmpz_poly_get_coeff_si(poly, 0);
			// While ">=", could work here, FORM does not use fast notation for INT_MIN
			if ( fast_coeff > INT_MIN && fast_coeff <= INT_MAX ) {
				*out++ = -SNUMBER;
				*out++ = (WORD)fast_coeff;
				return 2;
			}
		}
	}

	fmpz *coeff; // No init, we can get pointers directly to the poly's coefficients
	WORD *tmp_coeff = (WORD *)NumberMalloc("flint::poly_to_argument");

	// Create the inverse of var_map, so we don't have to search it for each symbol written
	map<unsigned,unsigned> var_map_inv;
	for (auto x: var_map) {
		var_map_inv[x.second] = x.first;
	}

	// TODO fast notation for single symbol arg

	WORD* arg_size = 0;
	WORD* arg_flag = 0;
	if ( with_arghead ) {
		arg_size = out++; // total arg size
		arg_flag = out++;
		*arg_flag = 0; // clean argument
	}
	else {
		// If we are not writing an ARGHEAD, we can still keep track of the
		// size of the sum of terms, but must not write to this location!
		arg_size = out;
	}

	// In reverse, since we want a "highfirst" output
	for (LONG i = n_terms-1; i >= 0; i--) {

		coeff = fmpz_poly_get_coeff_ptr(poly, i);

		// fmpz_poly is dense, there might be many zero coefficients:
		if ( !fmpz_is_zero(coeff) ) {

			// Convert the coefficient, write in temporary space
			const int coeff_size = flint::fmpz_get_form(coeff, tmp_coeff);

			// Now we have the coeff size, we can determine the output size
			// Check it fits if necessary: symbol code,power, num,den of "coeff_size", +1 for total coeff size
			unsigned current_size = prev_size + (out - arg_size) + 1 + 2*ABS(coeff_size) + 1;
			if ( i > 0 ) {
				// and also symbols header, code,power of the symbol
				current_size += 4;
			}
			if ( must_fit_term && (sizeof(WORD)*current_size > (size_t)AM.MaxTer) ) {
				MLOCK(ErrorMessageLock);
				MesPrint("flint::poly_to_argument: output exceeds MaxTermSize");
				MUNLOCK(ErrorMessageLock);
				Terminate(-1);
			}

			WORD* term_size = out++;

			if ( i > 0 ) {
				*out++ = SYMBOL;
				*out++ = 4; // The symbol array size, it is univariate
				*out++ = var_map_inv[0];
				*out++ = i;
			}

			// Copy numerator
			for ( int i = 0; i < ABS(coeff_size); i++ ) {
				*out++ = tmp_coeff[i];
			}
			*out++ = 1; // the denominator
			for (int i = 0; i < ABS(coeff_size)-1; i++) {
				*out++ = 0;
			}
			*out = 2*ABS(coeff_size) + 1; // the size of the coefficient
			if ( coeff_size < 0 ) { *out = -(*out); }
			out++;

			*term_size = out - term_size;
		}

	}

	if ( with_arghead ) {
		*arg_size = out - arg_size;
	}
	else {
		// with no arghead, we write a terminating zero
		*out++ = 0;
	}

	NumberFree(tmp_coeff, "flint::poly_to_argument");

	return out - arg_size;
}
/*
	#] flint::mpoly_to_argument :
	#[ flint::ratfun_add_mpoly :
*/
// Add the multi-variate FORM rational polynomials at t1 and t2. The result is written at out.
void flint::ratfun_add_mpoly(PHEAD WORD *t1, WORD *t2, WORD *out, const map<unsigned,unsigned> &var_map) {

	fmpz_mpoly_ctx_t ctx;
	fmpz_mpoly_ctx_init(ctx, var_map.size(), ORD_LEX);

	fmpz_mpoly_t gcd, num1, den1, num2, den2;
	fmpz_mpoly_init(gcd, ctx);
	fmpz_mpoly_init(num1, ctx);
	fmpz_mpoly_init(den1, ctx);
	fmpz_mpoly_init(num2, ctx);
	fmpz_mpoly_init(den2, ctx);

	flint::ratfun_read_mpoly(t1, num1, den1, var_map, ctx);
	flint::ratfun_read_mpoly(t2, num2, den2, var_map, ctx);

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
	*out++ = AR.PolyFun;
	WORD* args_size = out++;
	WORD* args_flag = out++;
	*args_flag = 0; // clean prf
	FILLFUN3(out); // Remainder of funhead, if it is larger than 3

	out += flint::mpoly_to_argument(BHEAD out, true, false, out-args_size, num1, var_map, ctx);
	out += flint::mpoly_to_argument(BHEAD out, true, false, out-args_size, den1, var_map, ctx);

	*args_size = out - args_size + 1; // The +1 is to include the function ID
	AT.WorkPointer = out;

	fmpz_mpoly_clear(num1, ctx);
	fmpz_mpoly_clear(den1, ctx);
	fmpz_mpoly_clear(num2, ctx);
	fmpz_mpoly_clear(den2, ctx);
	fmpz_mpoly_clear(gcd, ctx);
	fmpz_mpoly_ctx_clear(ctx);
}
/*
	#] flint::ratfun_add_mpoly :
	#[ flint::ratfun_add_poly :
*/
// Add the uni-variate FORM rational polynomials at t1 and t2. The result is written at out.
void flint::ratfun_add_poly(PHEAD WORD *t1, WORD *t2, WORD *out, const map<unsigned,unsigned> &var_map) {

	fmpz_poly_t gcd, num1, den1, num2, den2;
	fmpz_poly_init(gcd);
	fmpz_poly_init(num1);
	fmpz_poly_init(den1);
	fmpz_poly_init(num2);
	fmpz_poly_init(den2);

	flint::ratfun_read_poly(t1, num1, den1);
	flint::ratfun_read_poly(t2, num2, den2);

	if ( fmpz_poly_equal(den1, den2) == 0 ) {
		fmpz_poly_gcd(gcd, den1, den2);
		fmpz_poly_div(den1, den1, gcd);
		fmpz_poly_div(den2, den2, gcd);

		fmpz_poly_mul(num1, num1, den2);
		fmpz_poly_mul(num2, num2, den1);

		fmpz_poly_add(num1, num1, num2);
		fmpz_poly_mul(den1, den1, den2);
		fmpz_poly_mul(den1, den1, gcd);
	}
	else {
		fmpz_poly_add(num1, num1, num2);
	}

	// Finally divide out any common factors between the resulting num1, den1:
	fmpz_poly_gcd(gcd, num1, den1);
	fmpz_poly_div(num1, num1, gcd);
	fmpz_poly_div(den1, den1, gcd);

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
	*out++ = AR.PolyFun;
	WORD* args_size = out++;
	WORD* args_flag = out++;
	*args_flag = 0; // clean prf
	FILLFUN3(out); // Remainder of funhead, if it is larger than 3

	out += flint::poly_to_argument(BHEAD out, true, false, out-args_size, num1, var_map);
	out += flint::poly_to_argument(BHEAD out, true, false, out-args_size, den1, var_map);

	*args_size = out - args_size + 1; // The +1 is to include the function ID
	AT.WorkPointer = out;

	fmpz_poly_clear(num1);
	fmpz_poly_clear(den1);
	fmpz_poly_clear(num2);
	fmpz_poly_clear(den2);
	fmpz_poly_clear(gcd);
}
/*
	#] flint::ratfun_add_poly :
	#[ flint::ratfun_normalize_mpoly :
*/
// Multiply and simplify occurrences of the multi-variate FORM rational polynomials found in term.
// The final term is written in place, with the rational polynomial at the end.
void flint::ratfun_normalize_mpoly(PHEAD WORD *term, const map<unsigned,unsigned> &var_map) {

	// The length of the coefficient
	const int ncoeff = (term + *term)[-1];
	// The end of the term data, before the coefficient:
	const WORD *tstop = term + *term - ABS(ncoeff);

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
	flint::fmpz_set_form(tmpNum, (UWORD*)tstop, ncoeff/2);
	flint::fmpz_set_form(tmpDen, (UWORD*)tstop+ABS(ncoeff/2), ABS(ncoeff/2));
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
			flint::ratfun_read_mpoly(t, num2, den2, var_map, ctx);

			if ((t[2] & MUSTCLEANPRF) != 0) { // first normalize TODO do this inside read?
				fmpz_mpoly_t gcd;
				fmpz_mpoly_init(gcd, ctx);
				if (!fmpz_mpoly_gcd_cofactors(gcd, num2, den2, num2, den2, ctx)) {
					cout << "flint::ratfun_normalize: error in fmpz_mpoly_gcd_cofactors" << endl;
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

	// Fix sign: leading term of den should be positive.
	// Maybe not necessary, does gcd_cofactors already arrange for this somehow?
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

	out += flint::mpoly_to_argument(BHEAD out, true, false, out-args_size, num1, var_map, ctx);
	out += flint::mpoly_to_argument(BHEAD out, true, false, out-args_size, den1, var_map, ctx);

	*args_size = out - args_size + 1; // The +1 is to include the function ID

	// +3 for the coefficient of 1/1, which is added after the check
	if ( sizeof(WORD)*(*args_size+3) >= (size_t)AM.MaxTer ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::ratfun_normalize: output exceeds MaxTermSize");
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
}
/*
	#] flint::ratfun_normalize_mpoly :
	#[ flint::ratfun_normalize_poly :
*/
// Multiply and simplify occurrences of the uni-variate FORM rational polynomials found in term.
// The final term is written in place, with the rational polynomial at the end.
void flint::ratfun_normalize_poly(PHEAD WORD *term, const map<unsigned,unsigned> &var_map) {

	// The length of the coefficient
	const int ncoeff = (term + *term)[-1];
	// The end of the term data, before the coefficient:
	const WORD *tstop = term + *term - ABS(ncoeff);

	fmpz_poly_t num1, den1;
	// TODO might want to use init2, if we have some idea of the number of terms which appear
	// some lower bound can be determined in get_variables?
	fmpz_poly_init(num1);
	fmpz_poly_init(den1);
	// Start with "trivial" polynomials, and multiply in the num and den of the prf which appear.
	fmpz_t tmpNum, tmpDen;
	fmpz_init(tmpNum);
	fmpz_init(tmpDen);
	flint::fmpz_set_form(tmpNum, (UWORD*)tstop, ncoeff/2);
	flint::fmpz_set_form(tmpDen, (UWORD*)tstop+ABS(ncoeff/2), ABS(ncoeff/2));
	fmpz_poly_set_fmpz(num1, tmpNum);
	fmpz_poly_set_fmpz(den1, tmpDen);
	fmpz_clear(tmpNum);
	fmpz_clear(tmpDen);

	// Loop over the occurrences of PolyFun in the term, and multiply in to num1, den1.
	// s tracks where we are writing the non-PolyFun term data. The final PolyFun will
	// go at the end.
	WORD* term_size = term;
	WORD* s = term + 1;
	for (WORD *t=term+1; t<tstop;) {
		if (*t == AR.PolyFun) {

			fmpz_poly_t num2, den2;
			fmpz_poly_init(num2);
			fmpz_poly_init(den2);
			flint::ratfun_read_poly(t, num2, den2);

			if ((t[2] & MUSTCLEANPRF) != 0) { // first normalize TODO do this inside read?
				fmpz_poly_t gcd;
				fmpz_poly_init(gcd);
				fmpz_poly_gcd(gcd, num2, den2);
				fmpz_poly_div(num2, num2, gcd);
				fmpz_poly_div(den2, den2, gcd);
				fmpz_poly_clear(gcd);
			}

			// get gcd of num1,den2 and num2,den1 and then assemble
			fmpz_poly_t gcd;
			fmpz_poly_init(gcd);
			fmpz_poly_gcd(gcd, num1, den2);
			fmpz_poly_div(num1, num1, gcd);
			fmpz_poly_div(den2, den2, gcd);
			fmpz_poly_gcd(gcd, num2, den1);
			fmpz_poly_div(num2, num2, gcd);
			fmpz_poly_div(den1, den1, gcd);
			fmpz_poly_clear(gcd);

			fmpz_poly_mul(num1, num1, num2);
			fmpz_poly_mul(den1, den1, den2);
			fmpz_poly_clear(num2);
			fmpz_poly_clear(den2);

			t += t[1];
		}

		else {
			// Not a PolyFun, just copy or skip over
			int i = t[1];
			if (s != t) { NCOPY(s,t,i); }
			else { t += i; s += i; }
		}
	}

	// Fix sign: leading term of den should be positive. Since poly_to_argument writes
	// out the terms in high-first order, the "leading term" is the final term:
	if ( fmpz_sgn(fmpz_poly_get_coeff_ptr(den1, fmpz_poly_degree(den1))) == -1 ) {
		fmpz_poly_neg(num1, num1);
		fmpz_poly_neg(den1, den1);
	}

	// Result in FORM notation:
	WORD* out = s;
	*out++ = AR.PolyFun;
	WORD* args_size = out++;
	WORD* args_flag = out++;
	*args_flag &= ~MUSTCLEANPRF;

	out += flint::poly_to_argument(BHEAD out, true, false, out-args_size, num1, var_map);
	out += flint::poly_to_argument(BHEAD out, true, false, out-args_size, den1, var_map);

	*args_size = out - args_size + 1; // The +1 is to include the function ID

	// +3 for the coefficient of 1/1, which is added after the check
	if ( sizeof(WORD)*(*args_size+3) >= (size_t)AM.MaxTer ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::ratfun_normalize: output exceeds MaxTermSize");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	*out++ = 1;
	*out++ = 1;
	*out++ = 3; // the term's coefficient is now 1/1

	*term_size = out - term_size;

	fmpz_poly_clear(num1);
	fmpz_poly_clear(den1);
}
/*
	#] flint::ratfun_normalize_poly :
	#[ flint::ratfun_read_mpoly :
*/
// Read the multi-variate FORM rational polynomial at a and create fmpz_mpoly_t numerator and denominator
void flint::ratfun_read_mpoly(const WORD *a, fmpz_mpoly_t num, fmpz_mpoly_t den, const map<unsigned,unsigned> &var_map,
	fmpz_mpoly_ctx_t ctx) {

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

	// Polys to collect the "den of the num" and "den of the den".
	// Input can arrive like this when enabling the PolyRatFun or moving things into it.
	fmpz_mpoly_t den_num;
	fmpz_mpoly_init(den_num, ctx);
	fmpz_mpoly_t den_den;
	fmpz_mpoly_init(den_den, ctx);

	// Read the numerator
	flint::mpoly_from_argument(num, den_num, a, true, var_map, ctx);

	NEXTARG(a);
	if ( a < arg_stop ) {
		// Read the denominator
		flint::mpoly_from_argument(den, den_den, a, true, var_map, ctx);
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
	#] flint::ratfun_read_mpoly :
	#[ flint::ratfun_read_poly :
*/
// Read the uni-variate FORM rational polynomial at a and create fmpz_mpoly_t numerator and denominator
void flint::ratfun_read_poly(const WORD *a, fmpz_poly_t num, fmpz_poly_t den) {

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

	// Polys to collect the "den of the num" and "den of the den".
	// Input can arrive like this when enabling the PolyRatFun or moving things into it.
	fmpz_poly_t den_num;
	fmpz_poly_init(den_num);
	fmpz_poly_t den_den;
	fmpz_poly_init(den_den);

	// Read the numerator
	flint::poly_from_argument(num, den_num, a, true);

	NEXTARG(a);
	if ( a < arg_stop ) {
		// Read the denominator
		flint::poly_from_argument(den, den_den, a, true);
		NEXTARG(a);
	}
	else {
		// The denominator is 1
		cout << "implement this" << endl;
		Terminate(-1);
		// TODO I am not sure how to get here... the poly code seems buggy?
	}

	// Multiply the num by den_den and den by den_num:
	fmpz_poly_mul(num, num, den_den);
	fmpz_poly_mul(den, den, den_num);


	fmpz_poly_clear(den_num);
	fmpz_poly_clear(den_den);

	if (a < arg_stop) {
		MLOCK(ErrorMessageLock);
		MesPrint ((char*)"ERROR: PolyRatFun cannot have more than two arguments");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
}
/*
	#] flint::ratfun_read_poly :
*/
