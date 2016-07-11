#pragma once

namespace Rcpp {
    template<> SEXP wrap(const Gene &e);
    template <> Gene as( SEXP e ) ;
}
