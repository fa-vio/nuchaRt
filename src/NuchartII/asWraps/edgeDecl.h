#pragma once

namespace Rcpp {
    template<> SEXP wrap(const Edge &e);
    template <> Edge as( SEXP e ) ;
}
