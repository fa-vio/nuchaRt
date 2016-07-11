#pragma once

namespace Rcpp {
    // define template specialisations for as and wrap
    template<> SEXP wrap(const Edge &e) {
        List ret = List::create( Rcpp::Named("Gene1") = e.getGeneSymbol1(),                              
                                Rcpp::Named("Chr1") = e.getChromosome1(),
                                Rcpp::Named("Start1") = e.getStart1(),
                                Rcpp::Named("End1") = e.getEnd1(),
                                Rcpp::Named("Seq1") = e.getSeq1(),
                                Rcpp::Named("Gene2") = e.getGeneSymbol2(),
                                Rcpp::Named("Chr2") = e.getChromosome2(),
                                Rcpp::Named("Start2") = e.getStart2(),
                                Rcpp::Named("End2") = e.getEnd2(),
                                Rcpp::Named("Seq2") = e.getSeq2(),
                                Rcpp::Named("Weight") = e.getWeight(),
                                Rcpp::Named("Prob") = e.getProb(),
                                Rcpp::Named("HS1") = e.getHS1(), 
                                Rcpp::Named("HS2") = e.getHS2() );                              
        return wrap(ret);
    }

    template <> Edge as( SEXP e ) {
        Rcpp::List edL = Rcpp::as<Rcpp::List>(e);
        Edge ed;
        ed.setChromosome1( Rcpp::as<std::string>(edL["Chr1"]) );
        ed.setGeneSymbol1( Rcpp::as<std::string>(edL["Gene1"]) );
        ed.setStart1( Rcpp::as<long>( edL["Start1"]) );
        ed.setEnd1( Rcpp::as<long>( edL["End1"]) );
        ed.setSeq1( Rcpp::as<std::string>(edL["Seq1"]) );
        ed.setChromosome2( Rcpp::as<std::string>(edL["Chr2"]) );
        ed.setGeneSymbol2( Rcpp::as<std::string>(edL["Gene2"]) );
        ed.setStart2( Rcpp::as<long>( edL["Start2"]) );
        ed.setEnd2( Rcpp::as<long>( edL["End2"]) );
        ed.setSeq2( Rcpp::as<std::string>(edL["Seq2"]) );
        ed.setWeight( Rcpp::as<double>(edL["Weight"]) );
        ed.setProb( Rcpp::as<double>(edL["Prob"]) );
        ed.setHS1( Rcpp::as<long>(edL["HS1"]) );
        ed.setHS2( Rcpp::as<long>(edL["HS2"]) );
        
        return ed;
    }
}