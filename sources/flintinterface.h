#pragma once
/** @file flintinterface.h
 *
 *   Prototypes for functions in flintinterface.cc
 */

extern "C" {
#include "form3.h"
}


#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_poly.h>


#include <iostream>
#include <vector>
#include <map>
#include <cassert>


// The bits of std that are needed:
using std::cout;
using std::endl;
using std::map;
using std::vector;


namespace flint {

	WORD fmpz_get_form(fmpz_t, WORD *);
	void fmpz_set_form(fmpz_t, UWORD *, WORD);

	WORD* gcd_mpoly(PHEAD const WORD *, const WORD *, const WORD, const map<unsigned,unsigned> &);
	WORD* gcd_poly(PHEAD const WORD *, const WORD *, const WORD, const map<unsigned,unsigned> &);

	map<unsigned,unsigned> get_variables(const vector <WORD *> &, const bool, const bool);

	unsigned mpoly_from_argument(fmpz_mpoly_t, fmpz_mpoly_t, const WORD *, const bool, const map<unsigned,unsigned> &, const fmpz_mpoly_ctx_t);
	ULONG mpoly_to_argument(PHEAD WORD *, const bool, const bool, const ULONG, const fmpz_mpoly_t, const map<unsigned,unsigned> &, const fmpz_mpoly_ctx_t);

	unsigned poly_from_argument(fmpz_poly_t, fmpz_poly_t, const WORD *, const bool);
	ULONG poly_to_argument(PHEAD WORD *, const bool, const bool, const ULONG, const fmpz_poly_t, const map<unsigned,unsigned> &);

	void ratfun_add_mpoly(PHEAD WORD *, WORD *, WORD *, const map<unsigned,unsigned> &);
	void ratfun_add_poly(PHEAD WORD *, WORD *, WORD *, const map<unsigned,unsigned> &);

	void ratfun_normalize_mpoly(PHEAD WORD *, const map<unsigned,unsigned> &);
	void ratfun_normalize_poly(PHEAD WORD *, const map<unsigned,unsigned> &);

	void ratfun_read_mpoly(const WORD *, fmpz_mpoly_t, fmpz_mpoly_t, const map<unsigned,unsigned> &, fmpz_mpoly_ctx_t);
	void ratfun_read_poly(const WORD *, fmpz_poly_t, fmpz_poly_t);
}
