#pragma once
/** @file flintinterface.h
 *
 *   Prototypes for functions in flintinterface.cc
 */

extern "C" {
#include "form3.h"
}


#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>


#include <iostream>
#include <vector>
#include <map>
#include <cassert>


// The bits of std that are needed:
using std::cout;
using std::endl;
using std::map;
using std::vector;
using std::swap;


namespace flint {

	typedef std::map<unsigned,unsigned> var_map_t;

	WORD* factorize_mpoly(PHEAD WORD *, WORD *, const bool, const bool, const var_map_t &);
	WORD* factorize_poly(PHEAD WORD *, WORD *, const bool, const bool, const var_map_t &);

	void form_sort(PHEAD WORD *);

	unsigned from_argument_mpoly(fmpz_mpoly_t, fmpz_mpoly_t, const WORD *, const bool,
		const var_map_t &, const fmpz_mpoly_ctx_t);
	unsigned from_argument_poly(fmpz_poly_t, fmpz_poly_t, const WORD *, const bool);

	WORD fmpz_get_form(fmpz_t, WORD *);
	void fmpz_set_form(fmpz_t, UWORD *, WORD);

	WORD* gcd_mpoly(PHEAD const WORD *, const WORD *, const WORD, const var_map_t &);
	WORD* gcd_poly(PHEAD const WORD *, const WORD *, const WORD, const var_map_t &);

	var_map_t get_variables(const vector <WORD *> &, const bool, const bool);

	WORD* mul_mpoly(PHEAD const WORD *, const WORD *, const var_map_t &);
	WORD* mul_poly(PHEAD const WORD *, const WORD *, const var_map_t &);

	void ratfun_add_mpoly(PHEAD WORD *, WORD *, WORD *, const var_map_t &);
	void ratfun_add_poly(PHEAD WORD *, WORD *, WORD *, const var_map_t &);

	void ratfun_normalize_mpoly(PHEAD WORD *, const var_map_t &);
	void ratfun_normalize_poly(PHEAD WORD *, const var_map_t &);

	void ratfun_read_mpoly(const WORD *, fmpz_mpoly_t, fmpz_mpoly_t, const var_map_t &,
		fmpz_mpoly_ctx_t);
	void ratfun_read_poly(const WORD *, fmpz_poly_t, fmpz_poly_t);

	ULONG to_argument_mpoly(PHEAD WORD *, const bool, const bool, const bool, const ULONG,
		const fmpz_mpoly_t, const var_map_t &, const fmpz_mpoly_ctx_t);
	ULONG to_argument_mpoly(PHEAD WORD *, const bool, const bool, const bool, const ULONG,
		const fmpz_mpoly_t, const var_map_t &, const fmpz_mpoly_ctx_t, const fmpz_t);
	ULONG to_argument_poly(PHEAD WORD *, const bool, const bool, const bool, const ULONG,
		const fmpz_poly_t, const var_map_t &);
	ULONG to_argument_poly(PHEAD WORD *, const bool, const bool, const bool, const ULONG,
		const fmpz_poly_t, const var_map_t &, const fmpz_t);


	namespace util {

		void simplify_fmpz(fmpz_t, fmpz_t, fmpz_t);

	}

}
