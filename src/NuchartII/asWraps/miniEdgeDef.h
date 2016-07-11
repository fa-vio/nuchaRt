#pragma once

namespace Rcpp {
    // define template specialisations for as and wrap
    template<> SEXP wrap(const miniEdge &e) {
        List ret = List::create( Rcpp::Named("Start1") = e.st1,
                                Rcpp::Named("Start2") = e.st2,
                                Rcpp::Named("Weight") = e.score,
                                Rcpp::Named("Prob") = e.prob,
                                Rcpp::Named("HS1") = e.HS1, 
                                Rcpp::Named("HS2") = e.HS2 );                              
        return wrap(ret);
    }

    template <> miniEdge as( SEXP e ) {
        Rcpp::List edL = Rcpp::as<Rcpp::List>(e);
        miniEdge ed;

        ed.st1   = Rcpp::as<long>( edL["Start1"]);
        ed.st2   = Rcpp::as<long>( edL["Start2"]);
        ed.score = Rcpp::as<double>(edL["Weight"]);
        ed.prob  = Rcpp::as<double>(edL["Prob"]);
        ed.HS1   = Rcpp::as<long>(edL["HS1"]);
        ed.HS2   =  Rcpp::as<long>(edL["HS2"]);
        
        return ed;
    }
}