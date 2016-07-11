#pragma once

namespace Rcpp {
    template<> SEXP wrap(const miniEdge &e);
    template <> miniEdge as( SEXP e ) ;
}
