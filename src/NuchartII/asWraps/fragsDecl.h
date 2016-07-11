#pragma once

namespace Rcpp {
    template<> SEXP wrap(const Fragment &e);
    template <> Fragment as( SEXP e ) ;
}
