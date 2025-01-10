/** @file flintinterface.cc
 *
 *   Contains functions which interface with FLINT functions to perform arithmetic with uni and
 *   multivariate polynomials and perform factorization.
 */

extern "C" {
#include "form3.h"
}

#include "flintinterface.h"

/*
	#[ flint::factorize_mpoly :
*/
WORD* flint::factorize_mpoly(PHEAD WORD *argin, WORD *argout, const bool with_arghead,
	const bool is_fun_arg, const var_map_t &var_map) {

	fmpz_mpoly_ctx_t ctx;
	fmpz_mpoly_ctx_init(ctx, var_map.size(), ORD_LEX);

	fmpz_mpoly_t arg, den, base;
	fmpz_mpoly_init(arg, ctx);
	fmpz_mpoly_init(den, ctx);
	fmpz_mpoly_init(base, ctx);

	flint::from_argument_mpoly(arg, den, argin, with_arghead, var_map, ctx);
	// The denominator must be 1:
	if ( fmpz_mpoly_is_one(den, ctx) != 1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::factorize_mpoly error: den != 1");
		MUNLOCK(ErrorMessageLock);
		fmpz_mpoly_clear(den, ctx);
		Terminate(-1);
	}
	fmpz_mpoly_clear(den, ctx);


	// Now we can factor the mpoly:
	fmpz_mpoly_factor_t arg_fac;
	fmpz_mpoly_factor_init(arg_fac, ctx);
	fmpz_mpoly_factor(arg_fac, arg, ctx);
	const long num_factors = fmpz_mpoly_factor_length(arg_fac, ctx);


	// Construct the output. If argout is not NULL, we write the result there.
	// Otherwise, allocate memory.
	// The output is zero-terminated list of factors. If with_arghead, each has an arghead which
	// contains its size. Otherwise, the factors are zero separated.

	// We only need to determine the size of the output if we are allocating memory, but we need to
	// loop through the factors to check for an overall sign anyway. Do both together:
	int overall_sign = 1;

	// Initially 1, for the final trailing 0.
	unsigned long output_size = 1;

	for ( long i = 0; i < num_factors; i++ ) {
		fmpz_mpoly_factor_get_base(base, arg_fac, i, ctx);
		const long exponent = fmpz_mpoly_factor_get_exp_si(arg_fac, i, ctx);

		// poly_factorize makes sure that the highest power of the highest symbol has a positive
		// coefficient. This should be the last term of the mpoly. Check the sign:
		if ( fmpz_sgn(fmpz_mpoly_term_coeff_ref(base, fmpz_mpoly_length(base, ctx)-1, ctx)) == -1 ) {
			if ( exponent % 2 == 1 ) {
				overall_sign *= -1;
			}
		}

		const bool write = false;
		for ( long j = 0; j < exponent; j++ ) {
			output_size += (unsigned long)flint::to_argument_mpoly(BHEAD NULL, with_arghead,
				is_fun_arg, write, 0, base, var_map, ctx);
		}
	}
	if ( overall_sign == -1 ) {
		// Add space for a number and an ARGHEAD or zero separator
		output_size += 4 + with_arghead ? ARGHEAD : 1;
	}

	// Now make the allocation of necessary:
	if ( argout == NULL ) {
		argout = (WORD*)Malloc1(sizeof(WORD)*output_size, "flint::factorize_mpoly");
	}


	WORD* old_argout = argout;

	// If the overall sign is negative, first write a full-notation -1. It will be absorbed into the
	// overall factor in the content by the caller.
	if ( with_arghead ) {
		*argout++ = 4 + ARGHEAD;
		*argout++ = 0;
		for ( int i = 2; i < ARGHEAD; i++ ) {
			*argout++ = 0;
		}
	}
	*argout++ = 4; // term size
	*argout++ = 1; // numerator
	*argout++ = 1; // denominator
	*argout++ = -3; // coeff size, negative number
	if ( !with_arghead ) {
		*argout++ = 0; // factor separator
	}

	for ( long i = 0; i < num_factors; i++ ) {
		fmpz_mpoly_factor_get_base(base, arg_fac, i, ctx);
		const long exponent = fmpz_mpoly_factor_get_exp_si(arg_fac, i, ctx);

		// poly_factorize makes sure that the highest power of the highest symbol has a positive
		// coefficient. This should be the last term of the mpoly. Fix the sign:
		if ( fmpz_sgn(fmpz_mpoly_term_coeff_ref(base, fmpz_mpoly_length(base, ctx)-1, ctx)) == -1 ) {
			fmpz_mpoly_neg(base, base, ctx);
		}

		const bool write = true;
		for ( long j = 0; j < exponent; j++ ) {
			argout += flint::to_argument_mpoly(BHEAD argout, with_arghead, is_fun_arg, write,
				argout-old_argout, base, var_map, ctx);
		}
	}
	// Final trailing zero to denote the end of the factors.
	*argout++ = 0;

	fmpz_mpoly_factor_clear(arg_fac, ctx);
	fmpz_mpoly_clear(base, ctx);
	fmpz_mpoly_clear(arg, ctx);
	fmpz_mpoly_ctx_clear(ctx);

	return old_argout;
}
/*
	#] flint::factorize_mpoly :
	#[ flint::factorize_poly :
*/
WORD* flint::factorize_poly(PHEAD WORD *argin, WORD *argout, const bool with_arghead,
	const bool is_fun_arg, const var_map_t &var_map) {

	fmpz_poly_t arg, den;
	fmpz_poly_init(arg);
	fmpz_poly_init(den);

	flint::from_argument_poly(arg, den, argin, with_arghead);
	// The denominator must be 1:
	if ( fmpz_poly_is_one(den) != 1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::factorize_poly error: den != 1");
		MUNLOCK(ErrorMessageLock);
		fmpz_poly_clear(den);
		Terminate(-1);
	}
	fmpz_poly_clear(den);


	// Now we can factor the poly:
	fmpz_poly_factor_t arg_fac;
	fmpz_poly_factor_init(arg_fac);
	fmpz_poly_factor(arg_fac, arg);
	// fmpz_poly_factor_t lacks some convenience functions which fmpz_mpoly_factor_t has.
	// I have worked out how to get the factors by looking at how fmpz_poly_factor_print works.
	const long num_factors = arg_fac->num;


	// Construct the output. If argout is not NULL, we write the result there.
	// Otherwise, allocate memory.
	// The output is zero-terminated list of factors. If with_arghead, each has an arghead which
	// contains its size. Otherwise, the factors are zero separated.
	if ( argout == NULL ) {
		// First we need to determine the size of the output. This is the same procedure as the
		// loop below, but we don't write anything in to_argument_poly (write arg: false).
		// Initially 1, for the final trailing 0.
		unsigned long output_size = 1;

		for ( long i = 0; i < num_factors; i++ ) {
			fmpz_poly_struct* base = arg_fac->p + i;
			
			const long exponent = arg_fac->exp[i];

			const bool write = false;
			for ( long j = 0; j < exponent; j++ ) {
				output_size += (unsigned long)flint::to_argument_poly(BHEAD NULL, with_arghead,
					is_fun_arg, write, 0, base, var_map);
			}
		}

		argout = (WORD*)Malloc1(sizeof(WORD)*output_size, "flint::factorize_poly");
	}

	WORD* old_argout = argout;

	for ( long i = 0; i < num_factors; i++ ) {
		fmpz_poly_struct* base = arg_fac->p + i;
		
		const long exponent = arg_fac->exp[i];

		const bool write = true;
		for ( long j = 0; j < exponent; j++ ) {
			argout += flint::to_argument_poly(BHEAD argout, with_arghead, is_fun_arg, write,
				argout-old_argout, base, var_map);
		}
	}
	*argout = 0;


	fmpz_poly_factor_clear(arg_fac);
	fmpz_poly_clear(arg);

	return old_argout;
}
/*
	#] flint::factorize_poly :

	#[ flint::form_sort :
*/
// Sort terms using form's sorting routines. Uses a custom (faster) compare routine, since here
// only symbols can appear.
// This is a modified poly_sort from polywrap.cc.
void flint::form_sort(PHEAD WORD *terms) {

	if ( terms[0] < 0 ) {
		// Fast notation, there is nothing to do
		return;
	}

	const WORD oldsorttype = AR.SortType;
	AR.SortType = SORTHIGHFIRST;

	const WORD in_size = terms[0] - ARGHEAD;
	WORD out_size;

	if ( NewSort(BHEAD0) ) {
		Terminate(-1);
	}
	AR.CompareRoutine = (COMPAREDUMMY)(&CompareSymbols);

	// Make sure the symbols are in the right order within the terms
	for ( int i = ARGHEAD; i < terms[0]; i += terms[i] ) {
		if ( SymbolNormalize(terms+i) < 0 || StoreTerm(BHEAD terms+i) ) {
			AR.SortType = oldsorttype;
			AR.CompareRoutine = (COMPAREDUMMY)(&Compare1);
			LowerSortLevel();
			Terminate(-1);
		}
	}

	if ( ( out_size = EndSort(BHEAD terms+ARGHEAD, 1) ) < 0 ) {
		AR.SortType = oldsorttype;
		AR.CompareRoutine = (COMPAREDUMMY)(&Compare1);
		Terminate(-1);
	}

	// Check the final size
	if ( in_size != out_size  ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::form_sort: error: unexpected sorted arg length change %d->%d", in_size,
			out_size);
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	AR.SortType = oldsorttype;
	AR.CompareRoutine = (COMPAREDUMMY)(&Compare1);
	terms[1] = 0; // set dirty flag to zero
}
/*
	#] flint::form_sort :

	#[ flint::from_argument_mpoly :
*/
// Convert a FORM argument (or 0-terminated list of terms: with_arghead == false) to a
// (multi-variate) fmpz_mpoly_t poly. The "denominator" is return in denpoly, and contains the
// overall negative-power numeric and symbolic factor.
unsigned flint::from_argument_mpoly(fmpz_mpoly_t poly, fmpz_mpoly_t denpoly, const WORD *args,
	const bool with_arghead, const var_map_t &var_map, const fmpz_mpoly_ctx_t ctx) {

	// Some callers re-use their poly, denpoly to avoid calling init/clear unnecessarily.
	// Make sure they are 0 to begin.
	fmpz_mpoly_set_si(poly, 0, ctx);
	fmpz_mpoly_set_si(denpoly, 0, ctx);

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
	fmpz_t tmp;
	fmpz_init(tmp);
	unsigned long neg_exponents[var_map.size()] = {0};

	for ( const WORD* term = args; term < arg_stop; term += term[0] ) {
		const WORD* term_stop = term+term[0];
		const WORD  coeff_size = (term_stop)[-1];
		const WORD* symbol_stop = term_stop - ABS(coeff_size);
		const WORD* t = term;

		t++;
		if ( t == symbol_stop ) {
			// Just a number, no symbols
		}
		else {
			t++; // this entry is SYMBOL
			t++; // this entry just has the size of the symbol array, but we can use symbol_stop

			for ( const WORD* s = t; s < symbol_stop; s += 2 ) {
				if ( *(s+1) < 0 ) {
					neg_exponents[var_map.at(*s)] =
						MaX(neg_exponents[var_map.at(*s)], (unsigned long)(-(*(s+1))) );
				}
			}
		}

		// Now check for a denominator in the coefficient:
		if ( *(symbol_stop+ABS(coeff_size/2)) != 1 ) {
			flint::fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// Record the LCM of the coefficient denominators:
			fmpz_lcm(den_coeff, den_coeff, tmp);
		}

		if ( !with_arghead && *term_stop == 0 ) {
			// + 1 for the terminating 0
			arg_size = term_stop - args + 1;
			break;
		}

	}
	// Assemble denpoly.
	fmpz_mpoly_set_coeff_fmpz_ui(denpoly, den_coeff, neg_exponents, ctx);


	// For the term coefficients
	fmpz_t coeff;
	fmpz_init(coeff);

	for ( const WORD* term = args; term < arg_stop; term += term[0] ) {
		const WORD* term_stop = term+term[0];
		const WORD  coeff_size = (term_stop)[-1];
		const WORD* symbol_stop = term_stop - ABS(coeff_size);
		const WORD* t = term;

		unsigned long exponents[var_map.size()] = {0};

		t++; // skip over the total size entry
		if ( t == symbol_stop ) {
			// Just a number, no symbols
		}
		else {
			t++; // this entry is SYMBOL
			t++; // this entry just has the size of the symbol array, but we can use symbol_stop
			for ( const WORD* s = t; s < symbol_stop; s += 2 ) {
				exponents[var_map.at(*s)] = *(s+1);
			}
		}
		// Now read the coefficient
		flint::fmpz_set_form(coeff, (UWORD*)symbol_stop, coeff_size/2);

		// Multiply by denominator LCM
		fmpz_mul(coeff, coeff, den_coeff);

		// Shift by neg_exponents
		for ( unsigned i = 0; i < var_map.size(); i++ ) {
			exponents[i] += neg_exponents[i];
		}

		// Read the denominator if there is one, and divide it out of the coefficient
		if ( *(symbol_stop+ABS(coeff_size/2)) != 1 ) {
			flint::fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// By construction, this is an exact division
			fmpz_divexact(coeff, coeff, tmp);
		}

		// Push the term to the mpoly, remember to sort when finished! This is much faster than using
		// fmpz_mpoly_set_coeff_fmpz_ui when the terms arrive in the "wrong order".
		fmpz_mpoly_push_term_fmpz_ui(poly, coeff, exponents, ctx);

		if ( !with_arghead && *term_stop == 0 ) {
			break;
		}

	}
	// And now sort the mpoly
	fmpz_mpoly_sort_terms(poly, ctx);


	fmpz_clear(tmp);
	fmpz_clear(den_coeff);
	fmpz_clear(coeff);

	return arg_size;
}
/*
	#] flint::from_argument_mpoly :
	#[ flint::from_argument_poly :
*/
// Convert a FORM argument (or 0-terminated list of terms: with_arghead == false) to a
// (uni-variate) fmpz_poly_t poly. The "denominator" is return in denpoly, and contains the
// overall negative-power numeric and symbolic factor.
unsigned flint::from_argument_poly(fmpz_poly_t poly, fmpz_poly_t denpoly, const WORD *args,
	const bool with_arghead) {

	// Some callers re-use their poly, denpoly to avoid calling init/clear unnecessarily.
	// Make sure they are 0 to begin.
	fmpz_poly_set_si(poly, 0);
	fmpz_poly_set_si(denpoly, 0);

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
	fmpz_t tmp;
	fmpz_init(tmp);
	unsigned long neg_exponent = 0;
	for ( const WORD* term = args; term < arg_stop; term += term[0] ) {

		const WORD* term_stop = term+term[0];
		const WORD  coeff_size = (term_stop)[-1];
		const WORD* symbol_stop = term_stop - ABS(coeff_size);
		const WORD* t = term;

		t++; // skip over the total size entry
		if ( t == symbol_stop ) {
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
		if ( *(symbol_stop+ABS(coeff_size/2)) != 1 ) {
			flint::fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// Record the LCM of the coefficient denominators:
			fmpz_lcm(den_coeff, den_coeff, tmp);
		}

		if ( *term_stop == 0 ) {
			// + 1 for the terminating 0
			arg_size = term_stop - args + 1;
			break;
		}
	}
	// Assemble denpoly.
	fmpz_poly_set_coeff_fmpz(denpoly, neg_exponent, den_coeff);


	// For the term coefficients
	fmpz_t coeff;
	fmpz_init(coeff);

	for ( const WORD* term = args; term < arg_stop; term += term[0] ) {

		const WORD* term_stop = term+term[0];
		const WORD  coeff_size = (term_stop)[-1];
		const WORD* symbol_stop = term_stop - ABS(coeff_size);
		const WORD* t = term;

		unsigned long exponent = 0;

		t++; // skip over the total size entry
		if ( t == symbol_stop ) {
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
		if ( *(symbol_stop+ABS(coeff_size/2)) != 1 ) {
			flint::fmpz_set_form(tmp, (UWORD*)(symbol_stop+ABS(coeff_size/2)), ABS(coeff_size/2));
			// By construction, this is an exact division
			fmpz_divexact(coeff, coeff, tmp);
		}

		// Add the term to the poly
		fmpz_poly_set_coeff_fmpz(poly, exponent, coeff);

		if ( *term_stop == 0 ) {
			break;
		}
	}


	fmpz_clear(tmp);
	fmpz_clear(den_coeff);
	fmpz_clear(coeff);

	return arg_size;
}
/*
	#] flint::from_argument_poly :

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

	// Use fixed-size functions to get limb data where possible. These probably cover most real
	// cases. 
	if ( nlimbs == 1 ) {
		const unsigned long limb = fmpz_get_ui(z);
		a[0] = (WORD)(limb & 0xFFFFFFFF);
		na++;
		a[1] = (WORD)(limb >> BITSINWORD);
		if ( a[1] != 0 ) {
			na++;
		}
	}
	else if ( nlimbs == 2 ) {
		unsigned long limb_hi = 0, limb_lo = 0;
		fmpz_get_uiui(&limb_hi, &limb_lo, z);
		a[0] = (WORD)(limb_lo & 0xFFFFFFFF);
			na++;
		a[1] = (WORD)(limb_lo >> BITSINWORD);
			na++;
		a[2] = (WORD)(limb_hi & 0xFFFFFFFF);
			na++;
		a[3] = (WORD)(limb_hi >> BITSINWORD);
		if ( a[3] != 0 ) {
			na++;
		}
	}
	else {
		unsigned long limb_data[nlimbs] = {0};
		fmpz_get_ui_array(limb_data, nlimbs, z);
		for ( long i = 0; i < nlimbs; i++ ) {
			a[2*i] = (WORD)(limb_data[i] & 0xFFFFFFFF);
			na++;
			a[2*i+1] = (WORD)(limb_data[i] >> BITSINWORD);
			if ( a[2*i+1] != 0 || i < (nlimbs-1) ) {
				// The final limb might fit in a single 32bit WORD. Only
				// increment na if the final WORD is non zero.
				na++;
			}
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
// Create an fmpz directly from FORM's long integer representation. fmpz uses 64bit unsigned limbs,
// but FORM uses 32bit UWORDs on 64bit architectures so we can't use fmpz_set_ui_array directly.
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

	// Remove padding. FORM stores numerators and denominators with equal numbers of limbs but we
	// don't need to do this within the fmpz. It is not necessary to do this really, the fmpz
	// doesn't add zero limbs unnecessarily, but we might be able to avoid creating the limb_data
	// array below.
	while ( a[na-1] == 0 ) {
		na--;
	}

	// If the number fits in fixed-size fmpz_set functions, we don't need to use additional memory
	// to convert to ULONG. These probably cover most real cases.
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
WORD* flint::gcd_mpoly(PHEAD const WORD *a, const WORD *b, const WORD must_fit_term,
	const var_map_t &var_map) {

	fmpz_mpoly_ctx_t ctx;
	fmpz_mpoly_ctx_init(ctx, var_map.size(), ORD_LEX);

	fmpz_mpoly_t pa, pb, denpa, denpb, gcd;
	fmpz_mpoly_init(pa, ctx);
	fmpz_mpoly_init(pb, ctx);
	fmpz_mpoly_init(denpa, ctx);
	fmpz_mpoly_init(denpb, ctx);
	fmpz_mpoly_init(gcd, ctx);

	flint::from_argument_mpoly(pa, denpa, a, false, var_map, ctx);
	flint::from_argument_mpoly(pb, denpb, b, false, var_map, ctx);

	// denpa, denpb must be 1:
	if ( fmpz_mpoly_is_one(denpa, ctx) != 1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::gcd_mpoly: error: denpa != 1");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	if ( fmpz_mpoly_is_one(denpb, ctx) != 1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::gcd_mpoly: error: denpb != 1");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	fmpz_mpoly_gcd(gcd, pa, pb, ctx);

	// This is freed by the caller
	WORD *res;
	if ( must_fit_term ) {
		res = TermMalloc("flint::gcd_mpoly");
	}
	else {
		// Determine the size of the GCD by passing write = false.
		const bool with_arghead = false;
		const bool write = false;
		const ULONG prev_size = 0;
		const unsigned long gcd_size = (unsigned long)flint::to_argument_mpoly(BHEAD NULL,
			with_arghead, must_fit_term, write, prev_size, gcd, var_map, ctx);
		
		res = (WORD *)Malloc1(sizeof(WORD)*gcd_size, "flint::gcd_mpoly");
	}

	const bool with_arghead = false;
	const bool write = true;
	const ULONG prev_size = 0;
	flint::to_argument_mpoly(BHEAD res, with_arghead, must_fit_term, write, prev_size, gcd, var_map,
		ctx);

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
WORD* flint::gcd_poly(PHEAD const WORD *a, const WORD *b, const WORD must_fit_term,
	const var_map_t &var_map) {

	fmpz_poly_t pa, pb, denpa, denpb, gcd;
	fmpz_poly_init(pa);
	fmpz_poly_init(pb);
	fmpz_poly_init(denpa);
	fmpz_poly_init(denpb);
	fmpz_poly_init(gcd);

	flint::from_argument_poly(pa, denpa, a, false);
	flint::from_argument_poly(pb, denpb, b, false);

	// denpa, denpb must be 1:
	if ( fmpz_poly_is_one(denpa) != 1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::gcd_poly: error: denpa != 1");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	if ( fmpz_poly_is_one(denpb) != 1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::gcd_poly: error: denpb != 1");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	fmpz_poly_gcd(gcd, pa, pb);

	// This is freed by the caller
	WORD *res;
	if ( must_fit_term ) {
		res = TermMalloc("flint::gcd_poly");
	}
	else {
		// Determine the size of the GCD by passing write = false.
		const bool with_arghead = false;
		const bool write = false;
		const ULONG prev_size = 0;
		const unsigned long gcd_size = (unsigned long)flint::to_argument_poly(BHEAD NULL,
			with_arghead, must_fit_term, write, prev_size, gcd, var_map);

		res = (WORD *)Malloc1(sizeof(WORD)*gcd_size, "flint::gcd_poly");
	}

	const bool with_arghead = false;
	const ULONG prev_size = 0;
	const bool write = true;
	flint::to_argument_poly(BHEAD res, with_arghead, must_fit_term, write, prev_size, gcd, var_map);

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
// numerators, denominators or expressions from calls to gcd_ etc. Return this list as a map
// between indices and symbol codes.
// TODO FACTORSYMBOL last?
flint::var_map_t flint::get_variables(const vector <WORD *> &es, const bool with_arghead,
	const bool sort_vars) {

	int num_vars = 0;
	// To be used if we sort by highest degree, as the polu code does.
	vector<int> degrees;
	var_map_t var_map;

	// extract all variables
	for ( unsigned ei = 0; ei < es.size(); ei++ ) {
		WORD *e = es[ei];

		// fast notation
		if ( *e == -SNUMBER ) {
		}
		else if ( *e == -SYMBOL ) {
			if ( !var_map.count(e[1]) ) {
				var_map[e[1]] = num_vars++;
				degrees.push_back(1);
			}
		}
		// JD: Here we need to check for non-symbol/number terms in fast notation.
		else if ( *e < 0 ) {
			MLOCK(ErrorMessageLock);
			MesPrint("ERROR: polynomials and polyratfuns must contain symbols only");
			MUNLOCK(ErrorMessageLock);
			Terminate(1);
		}
		else {
			for ( int i = with_arghead ? ARGHEAD:0; with_arghead ? i < e[0]:e[i] != 0; i += e[i] ) {
				if ( i+1 < i+e[i]-ABS(e[i+e[i]-1]) && e[i+1] != SYMBOL ) {
					MLOCK(ErrorMessageLock);
					MesPrint("ERROR: polynomials and polyratfuns must contain symbols only");
					MUNLOCK(ErrorMessageLock);
					Terminate(1);
				}

				for ( int j = i+3; j<i+e[i]-ABS(e[i+e[i]-1]); j += 2 ) {
					if ( !var_map.count(e[j]) ) {
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

	if ( sort_vars ) {
		// bubble sort variables in decreasing order of degree
		// (this seems better for factorization)
		for ( unsigned i = 0; i < var_map.size(); i++ ) {
			for ( unsigned j = 0; j+1 < var_map.size(); j++ ) {
				if ( degrees[j] < degrees[j+1] ) {
					swap(degrees[j], degrees[j+1]);

					// Find the map keys associated with the values we want to swap
					unsigned j0 = 0;
					unsigned j1 = 0;
					for ( auto x: var_map ) {
						if ( x.second == j ) {
							j0 = x.first;
						}
						else if ( x.second == j+1 ) {
							j1 = x.first;
						}
					}
					swap(var_map.at(j0), var_map.at(j1));
				}
			}
		}
	}
	// Otherwise, sort lexicographically in FORM's ordering
	else {
		for ( unsigned i = 0; i < var_map.size(); i++ ) {
			for ( unsigned j = 0; j+1 < var_map.size(); j++ ) {
				unsigned j0 = 0;
				unsigned j1 = 0;
				for ( auto x: var_map ) {
					if ( x.second == j ) {
						j0 = x.first;
					}
					else if ( x.second == j+1 ) {
						j1 = x.first;
					}
				}
				if ( j0 > j1 ) {
					swap(var_map.at(j0), var_map.at(j1));
				}
			}
		}
	}

	return var_map;
}
/*
	#] flint::get_variables :

	#[ flint::mul_mpoly :
*/
// Return a pointer to a buffer containing the product of the 0-terminated term lists at a and b.
// For multi-variate cases.
WORD* flint::mul_mpoly(PHEAD const WORD *a, const WORD *b, const var_map_t &var_map) {

	fmpz_mpoly_ctx_t ctx;
	fmpz_mpoly_ctx_init(ctx, var_map.size(), ORD_LEX);

	fmpz_mpoly_t pa, pb, denpa, denpb;
	fmpz_mpoly_init(pa, ctx);
	fmpz_mpoly_init(pb, ctx);
	fmpz_mpoly_init(denpa, ctx);
	fmpz_mpoly_init(denpb, ctx);

	flint::from_argument_mpoly(pa, denpa, a, false, var_map, ctx);
	flint::from_argument_mpoly(pb, denpb, b, false, var_map, ctx);

	// denpa, denpb must be integers. Negative symbol powers have been converted to extra symbols.
	if ( fmpz_mpoly_is_fmpz(denpa, ctx) != 1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::mul_mpoly: error: denpa is non-constant");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	if ( fmpz_mpoly_is_fmpz(denpb, ctx) != 1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::mul_mpoly: error: denpb is non-constant");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// Multiply numerators, store result in pa
	fmpz_mpoly_mul(pa, pa, pb, ctx);
	// Multiply denominators, store result in denpa, and convert to an fmpz:
	fmpz_mpoly_mul(denpa, denpa, denpb, ctx);
	fmpz_t den;
	fmpz_init(den);
	fmpz_mpoly_get_fmpz(den, denpa, ctx);

	WORD* res;
	// First determine the result size, and malloc. The result should have no arghead. Here we use
	// the "scale" argument of to_argument_mpoly since den might not be 1.
	const bool with_arghead = false;
	bool write = false;
	const bool must_fit_term = false;
	const ULONG prev_size = 0;
	const unsigned long mul_size = (unsigned long)flint::to_argument_mpoly(BHEAD NULL,
		with_arghead, must_fit_term, write, prev_size, pa, var_map, ctx, den);
	res = (WORD*)Malloc1(sizeof(WORD)*mul_size, "flint::mul_mpoly");

	write = true;
	flint::to_argument_mpoly(BHEAD res, with_arghead, must_fit_term, write, prev_size, pa,
		var_map, ctx, den);

	fmpz_clear(den);
	fmpz_mpoly_clear(pa, ctx);
	fmpz_mpoly_clear(pb, ctx);
	fmpz_mpoly_clear(denpa, ctx);
	fmpz_mpoly_clear(denpb, ctx);
	fmpz_mpoly_ctx_clear(ctx);

	return res;
}
/*
	#] flint::mul_mpoly :
	#[ flint::mul_poly :
*/
// Return a pointer to a buffer containing the product of the 0-terminated term lists at a and b.
// For uni-variate cases.
WORD* flint::mul_poly(PHEAD const WORD *a, const WORD *b, const var_map_t &var_map) {

	fmpz_poly_t pa, pb, denpa, denpb;
	fmpz_poly_init(pa);
	fmpz_poly_init(pb);
	fmpz_poly_init(denpa);
	fmpz_poly_init(denpb);

	flint::from_argument_poly(pa, denpa, a, false);
	flint::from_argument_poly(pb, denpb, b, false);

	// denpa, denpb must be integers. Negative symbol powers have been converted to extra symbols.
	if ( fmpz_poly_degree(denpa) != 0 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::mul_poly: error: denpa is non-constant");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	if ( fmpz_poly_degree(denpb) != 0 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::mul_poly: error: denpb is non-constant");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// Multiply numerators, store result in pa
	fmpz_poly_mul(pa, pa, pb);
	// Multiply denominators, store result in denpa, and convert to an fmpz:
	fmpz_poly_mul(denpa, denpa, denpb);
	fmpz_t den;
	fmpz_init(den);
	fmpz_poly_get_coeff_fmpz(den, denpa, 0);

	WORD* res;
	// First determine the result size, and malloc. The result should have no arghead. Here we use
	// the "scale" argument of to_argument_poly since den might not be 1.
	const bool with_arghead = false;
	bool write = false;
	const bool must_fit_term = false;
	const ULONG prev_size = 0;
	const unsigned long mul_size = (unsigned long)flint::to_argument_poly(BHEAD NULL,
		with_arghead, must_fit_term, write, prev_size, pa, var_map, den);
	res = (WORD*)Malloc1(sizeof(WORD)*mul_size, "flint::mul_poly");

	write = true;
	flint::to_argument_poly(BHEAD res, with_arghead, must_fit_term, write, prev_size, pa,
		var_map, den);

	fmpz_clear(den);
	fmpz_poly_clear(pa);
	fmpz_poly_clear(pb);
	fmpz_poly_clear(denpa);
	fmpz_poly_clear(denpb);

	return res;
}
/*
	#] flint::mul_poly :

	#[ flint::ratfun_add_mpoly :
*/
// Add the multi-variate FORM rational polynomials at t1 and t2. The result is written at out.
void flint::ratfun_add_mpoly(PHEAD WORD *t1, WORD *t2, WORD *out, const var_map_t &var_map) {

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

	const bool with_arghead = true;
	const bool must_fit_term = true;
	const bool write = true;
	out += flint::to_argument_mpoly(BHEAD out, with_arghead, must_fit_term, write, out-args_size,
		num1, var_map, ctx);
	out += flint::to_argument_mpoly(BHEAD out, with_arghead, must_fit_term, write, out-args_size,
		den1, var_map, ctx);

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
void flint::ratfun_add_poly(PHEAD WORD *t1, WORD *t2, WORD *out, const var_map_t &var_map) {

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

	const bool with_arghead = true;
	const bool must_fit_term = true;
	const bool write = true;
	out += flint::to_argument_poly(BHEAD out, with_arghead, must_fit_term, write, out-args_size,
		num1, var_map);
	out += flint::to_argument_poly(BHEAD out, with_arghead, must_fit_term, write, out-args_size,
		den1, var_map);

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
void flint::ratfun_normalize_mpoly(PHEAD WORD *term, const var_map_t &var_map) {

	// The length of the coefficient
	const int ncoeff = (term + *term)[-1];
	// The end of the term data, before the coefficient:
	const WORD *tstop = term + *term - ABS(ncoeff);

	fmpz_mpoly_ctx_t ctx;
	fmpz_mpoly_ctx_init(ctx, var_map.size(), ORD_LEX);

	fmpz_mpoly_t num1, den1;
	fmpz_mpoly_t num2, den2, gcd;
	// TODO might want to use init2, if we have some idea of the number of terms which appear
	// some lower bound can be determined in get_variables?
	fmpz_mpoly_init(num1, ctx);
	fmpz_mpoly_init(den1, ctx);
	fmpz_mpoly_init(num2, ctx);
	fmpz_mpoly_init(den2, ctx);
	fmpz_mpoly_init(gcd, ctx);

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
	for ( WORD *t = term + 1; t < tstop; ) {
		if ( *t == AR.PolyFun ) {
			flint::ratfun_read_mpoly(t, num2, den2, var_map, ctx);

			// get gcd of num1,den2 and num2,den1 and then assemble
			fmpz_mpoly_gcd_cofactors(gcd, num1, den2, num1, den2, ctx);
			fmpz_mpoly_gcd_cofactors(gcd, num2, den1, num2, den1, ctx);

			fmpz_mpoly_mul(num1, num1, num2, ctx);
			fmpz_mpoly_mul(den1, den1, den2, ctx);

			t += t[1];
		}

		else {
			// Not a PolyFun, just copy or skip over
			int i = t[1];
			if ( s != t ) { NCOPY(s,t,i); }
			else { t += i; s += i; }
		}
	}

	// Fix sign: leading term of den should be positive.
	if ( fmpz_sgn(fmpz_mpoly_term_coeff_ref(den1, 0, ctx)) == -1 ) {
		fmpz_mpoly_neg(num1, num1, ctx);
		fmpz_mpoly_neg(den1, den1, ctx);
	}

	// Result in FORM notation:
	WORD* out = s;
	*out++ = AR.PolyFun;
	WORD* args_size = out++;
	WORD* args_flag = out++;
	*args_flag &= ~MUSTCLEANPRF;

	const bool with_arghead = true;
	const bool must_fit_term = true;
	const bool write = true;
	out += flint::to_argument_mpoly(BHEAD out, with_arghead, must_fit_term, write, out-args_size,
		num1, var_map, ctx);
	out += flint::to_argument_mpoly(BHEAD out, with_arghead, must_fit_term, write, out-args_size,
		den1, var_map, ctx);

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
	fmpz_mpoly_clear(num2, ctx);
	fmpz_mpoly_clear(den2, ctx);
	fmpz_mpoly_clear(gcd, ctx);
	fmpz_mpoly_ctx_clear(ctx);
}
/*
	#] flint::ratfun_normalize_mpoly :
	#[ flint::ratfun_normalize_poly :
*/
// Multiply and simplify occurrences of the uni-variate FORM rational polynomials found in term.
// The final term is written in place, with the rational polynomial at the end.
void flint::ratfun_normalize_poly(PHEAD WORD *term, const var_map_t &var_map) {

	// The length of the coefficient
	const int ncoeff = (term + *term)[-1];
	// The end of the term data, before the coefficient:
	const WORD *tstop = term + *term - ABS(ncoeff);

	fmpz_poly_t num1, den1;
	fmpz_poly_t num2, den2, gcd;
	// TODO might want to use init2, if we have some idea of the number of terms which appear
	// some lower bound can be determined in get_variables?
	fmpz_poly_init(num1);
	fmpz_poly_init(den1);
	fmpz_poly_init(num2);
	fmpz_poly_init(den2);
	fmpz_poly_init(gcd);

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
	for ( WORD *t = term + 1; t < tstop; ) {
		if ( *t == AR.PolyFun ) {
			flint::ratfun_read_poly(t, num2, den2);

			// get gcd of num1,den2 and num2,den1 and then assemble
			fmpz_poly_gcd(gcd, num1, den2);
			fmpz_poly_div(num1, num1, gcd);
			fmpz_poly_div(den2, den2, gcd);
			fmpz_poly_gcd(gcd, num2, den1);
			fmpz_poly_div(num2, num2, gcd);
			fmpz_poly_div(den1, den1, gcd);

			fmpz_poly_mul(num1, num1, num2);
			fmpz_poly_mul(den1, den1, den2);

			t += t[1];
		}

		else {
			// Not a PolyFun, just copy or skip over
			int i = t[1];
			if ( s != t ) { NCOPY(s,t,i); }
			else { t += i; s += i; }
		}
	}

	// Fix sign: leading term of den should be positive. Since to_argument_poly writes
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

	const bool with_arghead = true;
	const bool must_fit_term = true;
	const bool write = true;
	out += flint::to_argument_poly(BHEAD out, with_arghead, must_fit_term, write, out-args_size,
		num1, var_map);
	out += flint::to_argument_poly(BHEAD out, with_arghead, must_fit_term, write, out-args_size,
		den1, var_map);

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
	fmpz_poly_clear(num2);
	fmpz_poly_clear(den2);
	fmpz_poly_clear(gcd);
}
/*
	#] flint::ratfun_normalize_poly :

	#[ flint::ratfun_read_mpoly :
*/
// Read the multi-variate FORM rational polynomial at a and create fmpz_mpoly_t numerator and
// denominator.
void flint::ratfun_read_mpoly(const WORD *a, fmpz_mpoly_t num, fmpz_mpoly_t den,
	const var_map_t &var_map, fmpz_mpoly_ctx_t ctx) {

	// The end of the arguments:
	const WORD* arg_stop = a+a[1];

	const bool must_normalize = (a[2] & MUSTCLEANPRF) != 0;

	a += FUNHEAD;
	if ( a >= arg_stop ) {
		MLOCK(ErrorMessageLock);
		MesPrint("ERROR: PolyRatFun cannot have zero arguments");
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
	flint::from_argument_mpoly(num, den_num, a, true, var_map, ctx);
	NEXTARG(a);

	if ( a < arg_stop ) {
		// Read the denominator
		flint::from_argument_mpoly(den, den_den, a, true, var_map, ctx);
		NEXTARG(a);
	}
	else {
		// The denominator is 1
		MLOCK(ErrorMessageLock);
		MesPrint("implement this");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	if ( a < arg_stop ) {
		MLOCK(ErrorMessageLock);
		MesPrint("ERROR: PolyRatFun cannot have more than two arguments");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// Multiply the num by den_den and den by den_num:
	fmpz_mpoly_mul(num, num, den_den, ctx);
	fmpz_mpoly_mul(den, den, den_num, ctx);

	if ( must_normalize ) {
		fmpz_mpoly_t gcd;
		fmpz_mpoly_init(gcd, ctx);
		fmpz_mpoly_gcd_cofactors(gcd, num, den, num, den, ctx);
		fmpz_mpoly_clear(gcd, ctx);
	}


	fmpz_mpoly_clear(den_num, ctx);
	fmpz_mpoly_clear(den_den, ctx);
}
/*
	#] flint::ratfun_read_mpoly :
	#[ flint::ratfun_read_poly :
*/
// Read the uni-variate FORM rational polynomial at a and create fmpz_mpoly_t numerator and
// denominator.
void flint::ratfun_read_poly(const WORD *a, fmpz_poly_t num, fmpz_poly_t den) {

	// The end of the arguments:
	const WORD* arg_stop = a+a[1];

	const bool must_normalize = (a[2] & MUSTCLEANPRF) != 0;

	a += FUNHEAD;
	if ( a >= arg_stop ) {
		MLOCK(ErrorMessageLock);
		MesPrint("ERROR: PolyRatFun cannot have zero arguments");
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
	flint::from_argument_poly(num, den_num, a, true);
	NEXTARG(a);

	if ( a < arg_stop ) {
		// Read the denominator
		flint::from_argument_poly(den, den_den, a, true);
		NEXTARG(a);
	}
	else {
		// The denominator is 1
		MLOCK(ErrorMessageLock);
		MesPrint("implement this");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}
	if ( a < arg_stop ) {
		MLOCK(ErrorMessageLock);
		MesPrint("ERROR: PolyRatFun cannot have more than two arguments");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// Multiply the num by den_den and den by den_num:
	fmpz_poly_mul(num, num, den_den);
	fmpz_poly_mul(den, den, den_num);

	if ( must_normalize ) {
		fmpz_poly_t gcd;
		fmpz_poly_init(gcd);
		fmpz_poly_gcd(gcd, num, den);
		fmpz_poly_div(num, num, gcd);
		fmpz_poly_div(den, den, gcd);
		fmpz_poly_clear(gcd);
	}


	fmpz_poly_clear(den_num);
	fmpz_poly_clear(den_den);
}
/*
	#] flint::ratfun_read_poly :

	#[ flint::to_argument_mpoly :
*/
// Convert a fmpz_mpoly_t to a FORM argument (or 0-terminated list of terms: with_arghead==false).
// If the caller is building an output term, prev_size contains the size of the term so far, to
// check that the output fits in AM.MaxTer if must_fit_term.
// All coefficients will be divided by denscale (which might just be 1).
// If write is false, we never write to out but only track the total would-be size. This lets this
// function be repurposed as a "size of FORM notation" function without duplicating the code.
#define IFW(x) { if ( write ) x; }
ULONG flint::to_argument_mpoly(PHEAD WORD *out, const bool with_arghead, const bool must_fit_term,
	const bool write, const ULONG prev_size, const fmpz_mpoly_t poly, const var_map_t &var_map,
	const fmpz_mpoly_ctx_t ctx, const fmpz_t denscale) {

	// denscale should be positive, check this:
	if ( fmpz_sgn(denscale) == -1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::to_argument_mpoly: error: denscale < 0");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// out is modified later, keep the pointer at entry
	const WORD* out_entry = out;

	// Track the total size written. We could do this with pointer differences, but if
	// write == false we don't write to or move out to be able to find the size that way.
	ULONG ws = 0;

	// Check there is at least space for ARGHEAD WORDs (the arghead or fast-notation number/symbol)
	if ( write && must_fit_term && (sizeof(WORD)*(prev_size + ARGHEAD) > (size_t)AM.MaxTer) ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::to_argument_mpoly: output exceeds MaxTermSize");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// Create the inverse of var_map, so we don't have to search it for each symbol written
	var_map_t var_map_inv;
	for ( auto x: var_map ) {
		var_map_inv[x.second] = x.first;
	}

	LONG exponents[var_map.size()];
	const LONG n_terms = fmpz_mpoly_length(poly, ctx);

	if ( n_terms == 0 ) {
		if ( with_arghead ) {
			IFW(*out++ = -SNUMBER); ws++;
			IFW(*out++ = 0); ws++;
			return ws;
		}
		else {
			IFW(*out++ = 0); ws++;
			return ws;
		}
	}

	// For dividing out denscale
	fmpz_t coeff, den, gcd;
	fmpz_init(coeff);
	fmpz_init(den);
	fmpz_init(gcd);

	// The mpoly might be constant or a single symbol with coeff 1. Use fast notation if possible.
	if ( with_arghead && n_terms == 1 ) {

		if ( fmpz_mpoly_is_fmpz(poly, ctx) ) {
			// The mpoly is constant. Use fast notation if the number is an integer and small enough:

			fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, 0, ctx);
			fmpz_set(den, denscale);
			flint::util::simplify_fmpz(coeff, den, gcd);

			if ( fmpz_is_one(den) && fmpz_fits_si(coeff) ) {
				const long fast_coeff = fmpz_get_si(coeff);
				// While ">=", could work here, FORM does not use fast notation for INT_MIN
				if ( fast_coeff > INT_MIN && fast_coeff <= INT_MAX ) {
					IFW(*out++ = -SNUMBER); ws++;
					IFW(*out++ = (WORD)fast_coeff); ws++;

					fmpz_clear(coeff);
					fmpz_clear(den);
					fmpz_clear(gcd);

					return ws;
				}
			}
		}

		else {
			fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, 0, ctx);
			fmpz_set(den, denscale);
			flint::util::simplify_fmpz(coeff, den, gcd);

			if ( fmpz_is_one(coeff) && fmpz_is_one(den) ) {
				// The coefficient is one. Now check the symbol powers:

				fmpz_mpoly_get_term_exp_si(exponents, poly, 0, ctx);
				int use_fast = 0;
				WORD fast_symbol = 0;

				for ( unsigned i = 0; i < var_map.size(); i++ ) {
					if ( exponents[i] == 1 ) fast_symbol = var_map_inv[i];
					use_fast += exponents[i];
				}

				// use_fast has collected the total degree. If it is 1, then fast_symbol holds the code
				if ( use_fast == 1 ) {
					IFW(*out++ = -SYMBOL); ws++;
					IFW(*out++ = fast_symbol); ws++;

					fmpz_clear(coeff);
					fmpz_clear(den);
					fmpz_clear(gcd);

					return ws;
				}
			}
		}
	}

	WORD *tmp_coeff = (WORD *)NumberMalloc("flint::to_argument_mpoly");
	WORD *tmp_den = (WORD *)NumberMalloc("flint::to_argument_mpoly");

	WORD* arg_size = 0;
	WORD* arg_flag = 0;
	if ( with_arghead ) {
		IFW(arg_size = out++); ws++; // total arg size
		IFW(arg_flag = out++); ws++;
		IFW(*arg_flag = 0); // clean argument
	}

	for ( LONG i = 0; i < n_terms; i++ ) {

		fmpz_mpoly_get_term_exp_si(exponents, poly, i, ctx);

		fmpz_mpoly_get_term_coeff_fmpz(coeff, poly, i, ctx);
		fmpz_set(den, denscale);
		flint::util::simplify_fmpz(coeff, den, gcd);

		int num_symbols = 0;
		for ( unsigned i = 0; i < var_map.size(); i++ ) {
			if ( exponents[i] != 0 ) { num_symbols += 1; }
		}

		// Convert the coefficient, write in temporary space
		const int num_size = flint::fmpz_get_form(coeff, tmp_coeff);
		const int den_size = flint::fmpz_get_form(den, tmp_den);
		assert(den_size >= 1);
		// It is important that this comparison is >=, or we loose the sign!
		const int coeff_size = ABS(num_size) >= ABS(den_size) ? num_size : den_size;

		// Now we have the number of symbols and the coeff size, we can determine the output size.
		// Check it fits if necessary: term size, num,den of "coeff_size", +1 for total coeff size
		unsigned current_size = prev_size + ws + 1 + 2*ABS(coeff_size) + 1;
		if ( num_symbols ) {
			// symbols header, code,power of each symbol:
			current_size += 2 + 2*num_symbols;
		}
		if ( write && must_fit_term && (sizeof(WORD)*current_size > (size_t)AM.MaxTer) ) {
			MLOCK(ErrorMessageLock);
			MesPrint("flint::to_argument_mpoly: output exceeds MaxTermSize");
			MUNLOCK(ErrorMessageLock);
			Terminate(-1);
		}

		WORD* term_size = 0;
		IFW(term_size = out++); ws++;
		if ( num_symbols ) {
			IFW(*out++ = SYMBOL); ws++;
			WORD* symbol_size = 0;
			IFW(symbol_size = out++); ws++;
			IFW(*symbol_size = 2);

			for ( unsigned i = 0; i < var_map.size(); i++ ) {
				if ( exponents[i] != 0 ) {
					IFW(*out++ = var_map_inv[i]); ws++;
					IFW(*out++ = exponents[i]); ws++;
					IFW(*symbol_size += 2);
				}
			}
		}

		// Copy numerator
		for ( int i = 0; i < ABS(num_size); i++ ) {
			IFW(*out++ = tmp_coeff[i]); ws++;
		}
		for ( int i = ABS(num_size); i < ABS(coeff_size); i++ ) {
			IFW(*out++ = 0); ws++;
		}
		// Copy denominator
		for ( int i = 0; i < ABS(den_size); i++ ) {
			IFW(*out++ = tmp_den[i]); ws++;
		}
		for ( int i = ABS(den_size); i < ABS(coeff_size); i++ ) {
			IFW(*out++ = 0); ws++;
		}

		IFW(*out = 2*ABS(coeff_size) + 1); // the size of the coefficient
		IFW(if ( coeff_size < 0 ) { *out = -(*out); });
		IFW(out++); ws++;

		IFW(*term_size = out - term_size);
	}

	if ( with_arghead ) {
		IFW(*arg_size = out - arg_size);
		if ( write ) {
			// Sort into form highfirst ordering
			flint::form_sort(BHEAD (WORD*)(out_entry));
		}
	}
	else {
		// with no arghead, we write a terminating zero
		IFW(*out++ = 0); ws++;
	}

	NumberFree(tmp_coeff, "flint::to_argument_mpoly");
	NumberFree(tmp_den, "flint::to_argument_mpoly");

	fmpz_clear(coeff);
	fmpz_clear(den);
	fmpz_clear(gcd);

	return ws;
}

// If no denscale argument is supplied, just set it to 1 and call the usual function
ULONG flint::to_argument_mpoly(PHEAD WORD *out, const bool with_arghead, const bool must_fit_term,
	const bool write, const ULONG prev_size, const fmpz_mpoly_t poly, const var_map_t &var_map,
	const fmpz_mpoly_ctx_t ctx) {

	fmpz_t tmp;
	fmpz_init_set_ui(tmp, 1);

	ULONG ret = flint::to_argument_mpoly(BHEAD out, with_arghead, must_fit_term, write, prev_size,
		poly, var_map, ctx, tmp);

	fmpz_clear(tmp);
	return ret;
}
/*
	#] flint::to_argument_mpoly :
	#[ flint::to_argument_poly :
*/
// Convert a fmpz_poly_t to a FORM argument (or 0-terminated list of terms: with_arghead==false).
// If the caller is building an output term, prev_size contains the size of the term so far, to
// check that the output fits in AM.MaxTer if must_fit_term.
// All coefficients will be divided by denscale (which might just be 1).
// If write is false, we never write to out but only track the total would-be size. This lets this
// function be repurposed as a "size of FORM notation" function without duplicating the code.
ULONG flint::to_argument_poly(PHEAD WORD *out, const bool with_arghead, const bool must_fit_term,
	const bool write, const ULONG prev_size, const fmpz_poly_t poly, const var_map_t &var_map,
	const fmpz_t denscale) {

	// denscale should be positive, check this:
	if ( fmpz_sgn(denscale) == -1 ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::to_argument_poly: error: denscale < 0");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// Track the total size written. We could do this with pointer differences, but if
	// write == false we don't write to or move out to be able to find the size that way.
	ULONG ws = 0;

	// Check there is at least space for ARGHEAD WORDs (the arghead or fast-notation number/symbol)
	if ( write && must_fit_term && (sizeof(WORD)*(prev_size + ARGHEAD) > (size_t)AM.MaxTer) ) {
		MLOCK(ErrorMessageLock);
		MesPrint("flint::to_argument_poly: output exceeds MaxTermSize");
		MUNLOCK(ErrorMessageLock);
		Terminate(-1);
	}

	// Create the inverse of var_map, so we don't have to search it for each symbol written
	var_map_t var_map_inv;
	for ( auto x: var_map ) {
		var_map_inv[x.second] = x.first;
	}

	const LONG n_terms = fmpz_poly_length(poly);

	// The poly is zero
	if ( n_terms == 0 ) {
		if ( with_arghead ) {
			IFW(*out++ = -SNUMBER); ws++;
			IFW(*out++ = 0); ws++;
			return ws;
		}
		else {
			IFW(*out++ = 0); ws++;
			return ws;
		}
	}

	// For dividing out denscale
	fmpz_t coeff, den, gcd;
	fmpz_init(coeff);
	fmpz_init(den);
	fmpz_init(gcd);

	// The poly is constant, use fast notation if the coefficient is integer and small enough
	if ( with_arghead && n_terms == 1 ) {

		fmpz_poly_get_coeff_fmpz(coeff, poly, 0);
		fmpz_set(den, denscale);
		flint::util::simplify_fmpz(coeff, den, gcd);

		if ( fmpz_is_one(den) && fmpz_fits_si(coeff) ) {
			const long fast_coeff = fmpz_get_si(coeff);
			// While ">=", could work here, FORM does not use fast notation for INT_MIN
			if ( fast_coeff > INT_MIN && fast_coeff <= INT_MAX ) {
				IFW(*out++ = -SNUMBER); ws++;
				IFW(*out++ = (WORD)fast_coeff); ws++;

				fmpz_clear(coeff);
				fmpz_clear(den);
				fmpz_clear(gcd);

				return ws;
			}
		}
	}

	// The poly might be a single symbol with coeff 1, use fast notation if so.
	if ( with_arghead && n_terms == 2 ) {
		if ( fmpz_is_zero(fmpz_poly_get_coeff_ptr(poly, 0)) ) {
			// The constant term is zero

			fmpz_poly_get_coeff_fmpz(coeff, poly, 1);
			fmpz_set(den, denscale);
			flint::util::simplify_fmpz(coeff, den, gcd);

			if ( fmpz_is_one(coeff) && fmpz_is_one(den) ) {
				// Single symbol with coeff 1. Use fast notation:
				IFW(*out++ = -SYMBOL); ws++;
				IFW(*out++ = var_map_inv[0]); ws++;

				fmpz_clear(coeff);
				fmpz_clear(den);
				fmpz_clear(gcd);

				return ws;
			}
		}
	}

	WORD *tmp_coeff = (WORD *)NumberMalloc("flint::to_argument_poly");
	WORD *tmp_den = (WORD *)NumberMalloc("flint::to_argument_mpoly");

	WORD* arg_size = 0;
	WORD* arg_flag = 0;
	if ( with_arghead ) {
		IFW(arg_size = out++); ws++; // total arg size
		IFW(arg_flag = out++); ws++;
		IFW(*arg_flag = 0); // clean argument
	}

	// In reverse, since we want a "highfirst" output
	for ( LONG i = n_terms-1; i >= 0; i-- ) {

		// fmpz_poly is dense, there might be many zero coefficients:
		if ( !fmpz_is_zero(fmpz_poly_get_coeff_ptr(poly, i)) ) {

			fmpz_poly_get_coeff_fmpz(coeff, poly, i);
			fmpz_set(den, denscale);
			flint::util::simplify_fmpz(coeff, den, gcd);

			// Convert the coefficient, write in temporary space
			const int num_size = flint::fmpz_get_form(coeff, tmp_coeff);
			const int den_size = flint::fmpz_get_form(den, tmp_den);
			assert(den_size >= 1);
			// It is important that this comparison is >=, or we loose the sign!
			const int coeff_size = ABS(num_size) >= ABS(den_size) ? num_size : den_size;

			// Now we have the coeff size, we can determine the output size
			// Check it fits if necessary: symbol code,power, num,den of "coeff_size",
			// +1 for total coeff size
			unsigned current_size = prev_size + ws + 1 + 2*ABS(coeff_size) + 1;
			if ( i > 0 ) {
				// and also symbols header, code,power of the symbol
				current_size += 4;
			}
			if ( write && must_fit_term && (sizeof(WORD)*current_size > (size_t)AM.MaxTer) ) {
				MLOCK(ErrorMessageLock);
				MesPrint("flint::to_argument_poly: output exceeds MaxTermSize");
				MUNLOCK(ErrorMessageLock);
				Terminate(-1);
			}

			WORD* term_size = 0;
			IFW(term_size = out++); ws++;

			if ( i > 0 ) {
				IFW(*out++ = SYMBOL); ws++;
				IFW(*out++ = 4); ws++; // The symbol array size, it is univariate
				IFW(*out++ = var_map_inv[0]); ws++;
				IFW(*out++ = i); ws++;
			}

			// Copy numerator
			for ( int i = 0; i < ABS(num_size); i++ ) {
				IFW(*out++ = tmp_coeff[i]); ws++;
			}
			for ( int i = ABS(num_size); i < ABS(coeff_size); i++ ) {
				IFW(*out++ = 0); ws++;
			}
			// Copy denominator
			for ( int i = 0; i < ABS(den_size); i++ ) {
				IFW(*out++ = tmp_den[i]); ws++;
			}
			for ( int i = ABS(den_size); i < ABS(coeff_size); i++ ) {
				IFW(*out++ = 0); ws++;
			}

			IFW(*out = 2*ABS(coeff_size) + 1); // the size of the coefficient
			IFW(if ( coeff_size < 0 ) { *out = -(*out); });
			IFW(out++); ws++;

			IFW(*term_size = out - term_size);
		}

	}

	if ( with_arghead ) {
		IFW(*arg_size = out - arg_size);
	}
	else {
		// with no arghead, we write a terminating zero
		IFW(*out++ = 0); ws++;
	}

	NumberFree(tmp_coeff, "flint::to_argument_poly");
	NumberFree(tmp_den, "flint::to_argument_poly");

	fmpz_clear(coeff);
	fmpz_clear(den);
	fmpz_clear(gcd);

	return ws;
}

// If no denscale argument is supplied, just set it to 1 and call the usual function
ULONG flint::to_argument_poly(PHEAD WORD *out, const bool with_arghead, const bool must_fit_term,
	const bool write, const ULONG prev_size, const fmpz_poly_t poly, const var_map_t &var_map) {

	fmpz_t tmp;
	fmpz_init_set_ui(tmp, 1);

	ULONG ret = flint::to_argument_poly(BHEAD out, with_arghead, must_fit_term, write, prev_size,
		poly, var_map, tmp);

	fmpz_clear(tmp);
	return ret;
}
/*
	#] flint::to_argument_mpoly :
*/

// Utility functions
/*
	#[ flint::util::simplify_fmpz :
*/
// Divide the GCD out of num and den
inline void flint::util::simplify_fmpz(fmpz_t num, fmpz_t den, fmpz_t gcd) {
	fmpz_gcd(gcd, num, den);
	fmpz_divexact(num, num, gcd);
	fmpz_divexact(den, den, gcd);
}
/*
	#] flint::util::simplify_fmpz :
*/
