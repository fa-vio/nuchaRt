
#pragma once

namespace Rcpp {
    // Foo
	template<> SEXP wrap(const SamData &s);
	template<> SamData as(SEXP s);
}