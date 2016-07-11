#pragma once

namespace Rcpp {
    
    // define template specialisations for as and wrap
    template<> SEXP wrap(const Fragment &f) {
        List ret = Rcpp::List::create( Rcpp::Named("Chr") = f.getChromName(),
                                   Rcpp::Named("Start") = f.getStart(),
                                   Rcpp::Named("Stop") = f.getStop(),
                                   Rcpp::Named("Id") = f.getFragNum() );
        return wrap(ret);
    }

    template <> Fragment as( SEXP f ) {
        Rcpp::List frL = Rcpp::as<Rcpp::List>(f);
        Fragment fr;
        fr.setChromosome( Rcpp::as<std::string>(frL["Chr"]) );
        fr.setStart( Rcpp::as<long>( frL["Start"]) );
        fr.setStop( Rcpp::as<long>( frL["Stop"]) );
        fr.setFragNum( Rcpp::as<long>( frL["Id"]) );
        
        return fr;
    }

}