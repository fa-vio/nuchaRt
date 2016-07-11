#pragma once

namespace Rcpp {
    
    // define template specialisations for as and wrap
    template<> SEXP wrap(const Gene& g) {
        List ret = Rcpp::List::create( Rcpp::Named("Id") = g.getPosition(), 
                                    Rcpp::Named("Chr") = g.getChr(),
                                   Rcpp::Named("Symbol") = g.getSymbol(),
                                   Rcpp::Named("Start") = g.getStart(),
                                   Rcpp::Named("Stop") = g.getStop(),
                                   Rcpp::Named("E_ID") = g.getEID(),
                                   Rcpp::Named("ExtStart") = g.getExtendedStart(),
                                   Rcpp::Named("ExtStopt") = g.getExtendedStop(),
                                   Rcpp::Named("NumIntrs") = g.getIntr(),
                                   Rcpp::Named("NumCis") = g.getNumCis(),
                                   Rcpp::Named("NumTrans") = g.getNumTrans(),
                                   Rcpp::Named("HS") = g.getHS() );
        return wrap(ret);
    }

    template <> Gene as( SEXP g ) {
        Rcpp::List gnL = Rcpp::as<Rcpp::List>(g);
        Gene gn;
        gn.setChromosomeName( Rcpp::as<std::string>(gnL["Chr"]) );
        gn.setSymbol( Rcpp::as<std::string>(gnL["Symbol"]) );
        gn.setStartCoordinate( Rcpp::as<long>(gnL["Start"]) );
        gn.setStopCoordinate( Rcpp::as<long>(gnL["Stop"]) );
        gn.setEID( Rcpp::as<long>(gnL["E_ID"]) );
        gn.setExtendedStart( Rcpp::as<long>(gnL["ExtStart"]) );
        gn.setExtendedStop( Rcpp::as<long>(gnL["ExtStop"]) );
        gn.setCis( Rcpp::as<long>(gnL["NumCis"]) );
        gn.setTrans( Rcpp::as<long>(gnL["NumTrans"]) );
        gn.setHS( Rcpp::as<long>(gnL["HS"]) );
        gn.setPosition( Rcpp::as<long>(gnL["Position"]) );
        
        long intrs = Rcpp::as<long>( gnL["NumIntrs"] );
        for(unsigned i=0; i<intrs; ++i)
            gn.setIntr();
            
        return gn;
    }
}