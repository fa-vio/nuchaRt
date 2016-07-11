
#pragma once

namespace Rcpp {
    
    // define template specialisations for as and wrap
    template<> SEXP wrap(const SamData &s) {
        List ret = Rcpp::List::create(  Rcpp::Named("Id")   = s.getId(),
                                        Rcpp::Named("Chr1") = s.getChr1(),
                                        Rcpp::Named("Start1") = s.getStart1(),
                                        Rcpp::Named("Chr2") = s.getChr2(),
                                        Rcpp::Named("Start2") = s.getStart2(),
                                        Rcpp::Named("Seq") = s.getSeq(),
                                        Rcpp::Named("HS1")  = s.getHS1(),
                                        Rcpp::Named("HS2")  = s.getHS2() 
                                        );
        return Rcpp::wrap(ret);
    }

    template <> SamData as( SEXP s ) {
        Rcpp::List samL = Rcpp::as<Rcpp::List>(s);
        SamData sam;
                        
        sam.setId( Rcpp::as<long>( samL["Id"]) );
        sam.setChr1( Rcpp::as<std::string>(samL["Chr1"]) );
        sam.setStart1( Rcpp::as<long>( samL["Start1"]) );
        sam.setChr2( Rcpp::as<std::string>(samL["Chr2"]) );
        sam.setStart2( Rcpp::as<long>( samL["Start2"]) );
        sam.setSeq( Rcpp::as<std::string>(samL["Seq"]) );
        sam.setHS1( Rcpp::as<size_t>(samL["HS1"]) );
        sam.setHS2( Rcpp::as<size_t>(samL["HS2"]) );
    
        return sam;
    }
}