#include "NuchartII/common.hpp"
#include "NuchartII/Fragment.hpp"
#include "NuchartII/SamData.hpp"
#include "NuchartII/Gene.hpp"
#include "NuchartII/Edge.hpp"

#include <RcppCommon.h>
#include "NuchartII/asWraps/samDataDecl.h"
#include "NuchartII/asWraps/geneDecl.h"
#include "NuchartII/asWraps/edgeDecl.h"
#include "NuchartII/asWraps/fragsDecl.h"
#include "NuchartII/asWraps/miniEdgeDecl.h"
#include <Rcpp.h>
#include "NuchartII/asWraps/samDataDef.h"
#include "NuchartII/asWraps/geneDef.h"
#include "NuchartII/asWraps/edgeDef.h"
#include "NuchartII/asWraps/fragsDef.h"
#include "NuchartII/asWraps/miniEdgeDef.h"

#include "NuchartII/Finder.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
SEXP runNuChart2(SEXP gene, std::string samfile, std::string genefile, std::string fragsfile,
                std::string collsfile, std::string expfile, std::string uname, int nwrks,
                int lvl, int grain, int p, int colls, int only) {

    int nw=nwrks, level=lvl, chunk=grain, plot=p, intrs=colls, exprs=0;
	float e_lm=0.0;
	uint_64 st=0, sp=0;
	std::string coords, chr, crd, clstr, mld="";
	std::string::size_type found;

	std::string genes = genefile;
	std::string frags = fragsfile;
	std::string sam = samfile;
	std::string locs = collsfile;

	std::unordered_map<std::string, std::string> clusters;
	clustersMap(clusters);
	CharacterVector gn(gene);
	std::vector<std::string> gvec(gn.begin(), gn.end());

	if( gvec.size() == 1 ) {
	    if (gvec[0].compare(0, 3, "chr") == 0) {
	        coords = gvec[0];
	        clstr = gvec[0];
	    } else if(clusters.count(gvec[0]) > 0) {
	        coords = clusters[gvec[0]];
	        clstr = gvec[0];
	    }
	}

	if(intrs && locs.empty()) {
		Rcpp::stop("\nERROR: Interactions detection requires a collisions file to be specified. Aborting...");
        //exit(EXIT_FAILURE);
	}

	int cpus = getNumCpus();
	Rcpp::Rcout << "\n System configuration:\n"
             << "   CPUs: " << cpus << " (physical)\n"
             << "   RAM:  " << (getTotRAM() >> 30) << " GB (total)" << std::endl;

	bool active=false;
	if(cpus > 8) active=true;

	Finder *f = new Finder();
	f->parseFiles(genes, frags, sam);
	if(intrs) {
	  f->detectIntersections(locs);
      found = locs.find("/mld");
	  mld = locs.substr(found+1, 5);
	}

	if(!coords.empty()) {
	    found = coords.rfind(":");
	    chr = coords.substr(0, found);
	    crd = coords.substr(found+1, coords.length());
	    eraseTrailingWhiteSpaces(crd);
	    found = crd.rfind(",");
	    st = atol( crd.substr(0, found).c_str() );
	    sp = atol( crd.substr(found+1, crd.length()).c_str() );
	    gvec.clear();

		std::vector<Gene> gc = f->getGenesByCoordinates(chr, st, sp);
		for(unsigned g=0; g<gc.size(); ++g)
			gvec.push_back( gc[g].getSymbol() );
	}

	Rcpp::Rcout << "\n Using " << nw << " worker(s) | Search distance: " << level << "\n"
           << " Grain of parallel execution: " << chunk << "\n"
           << " SAM file:\t\t\'" << sam << "\'\n"
           << " Genes file:\t\t\'" << genes << "\'\n"
           << " Fragments file:\t\'" << frags << "\'\n"
           << " Intersections:\t\t\'" << locs << "\'\n";

	if(!expfile.empty()) {
	    if(exist_f(expfile.c_str())) {
	        f->parseExpressions(expfile);
	        exprs = 1;
	    } else {
	        Rcpp::Rcout << " Genes Expression file not found!\n";
	    }
	}

	Rcpp::Rcout << " Genes Expression:\t'" << expfile << "\'\n\n"
                << " Starting Gene(s): ";

	for(unsigned i=0; i<gvec.size(); ++i)
		Rcpp::Rcout << gvec[i] << " ";
	Rcpp::Rcout << "\n" << std::endl;

	// fix SAM filename for plotting utility
	found = sam.find("ord_");
	if(found == std::string::npos) {
	  found = sam.rfind("/");
	  sam = sam.substr(found+1, sam.length());
	} else sam = sam.substr(found+4, sam.length());

	found = sam.find(".sam");
	sam = sam.substr(0, found);

    std::vector<Edge> r_ed;

	int ret = f->parFindConnections(gvec, r_ed, sam, level, nw, clstr, chunk, e_lm,
                                    plot, intrs, exprs, uname, mld, active, only);

	//f->~Finder();
	delete f;
	return wrap( r_ed );
}

// [[Rcpp::export]]
SEXP runNormalisation(Rcpp::DataFrame mEd, std::string samfile, std::string featfile,
                      int nwrks, float e_lm, int chunk, int numnodes) {
    Finder *f = new Finder();
    f->reducedParseFiles(samfile, featfile);

    Rcpp::Rcout << " [ContactMaps] Building Contact Maps for the given SAM file..." << std::endl;
    double cmt = tmn::elapsedTime(0);
    size_t nCMaps = f->buildCMaps();
    cmt = tmn::elapsedTime(1);
    Rcpp::Rcout << " [ContactMaps] Done. " << nCMaps << " ContactMaps built (" << cmt << "ms) | ";
    long rsz = ( (48+sizeof(umCM))+(nCMaps*(sizeof(shPtrCM)+sizeof(uint_64)+32)) ) >> 10;
    Rcpp::Rcout << rsz << " KB" << std::endl;

    std::vector<miniEdge> miniEd;
    size_t numRows = mEd.nrows();
    miniEd.reserve(numRows);
    size_t m_f = f->getMaxFeatQtValue();

    int cpus = getNumCpus();
    bool active=false;
    if(cpus > 8) active=true;

    std::vector<long> st1   = Rcpp::as< std::vector<long> >(mEd[0]);
    std::vector<long> st2   = Rcpp::as< std::vector<long> >(mEd[1]);
    std::vector<double> wgt = Rcpp::as< std::vector<double> >(mEd[2]);
    std::vector<long> hs1   = Rcpp::as< std::vector<long> >(mEd[4]);
    std::vector<long> hs2   = Rcpp::as< std::vector<long> >(mEd[5]);

    for(int i=0; i<numRows; ++i) {
        miniEdge m_ed(hs1[i], hs2[i], st1[i], st2[i], false, wgt[i]);
        miniEd.push_back( m_ed );
    }

    double **gcc  = new double*[nwrks];
    double **map  = new double*[nwrks];
    uint_64 **len = new uint_64*[nwrks];

    for(unsigned f=0; f<nwrks; ++f) {
        gcc[f] = new double[2*m_f]();
        map[f] = new double[2*m_f]();
        len[f] = new uint_64[2*m_f]();
    }

    // duplicate SamData into reduced, memory efficient, SamDataT*
    std::vector<SamDataT*> *samt;
    void *pt;

#ifdef USE_NUMA
    if(numa_available >= 0) {
        pt = numa_alloc_interleaved(sizeof(SamDataT*)+24);
        samt = new (pt) std::vector<SamDataT*>();
    } else {
        pt = ::malloc(sizeof(SamDataT*)+24);
        samt = new (pt) std::vector<SamDataT*>();
    }
#else
    pt = ::malloc(sizeof(SamDataT*)+24);
    samt = new (pt) std::vector<SamDataT*>();
#endif

    size_t sam_size = f->samVector()->size();
    samt->reserve(sam_size);
    for(unsigned i=0; i<sam_size; ++i)
        samt->push_back(f->samVector()->at(i).toSamDataT());

    ParallelFor ffpf(nwrks, true /*spin-wait*/, true /*spin-barrier*/);
    ffpf.disableScheduler(active);
    //ffpf.parallel_for(0, numRows, 1, chunk, [] (const unsigned j) {j;}); // warm up

    Rcpp::Rcout << "\n -- START SCORE COMPUTATION FOR " << numRows << "EDGES..." << std::endl;
    double et2 = tmn::elapsedTime(0);
    ffpf.parallel_for_thid(0, numRows, 1, chunk, [&] (const unsigned j, const int thid) {
        f->calculateEdgeScore(miniEd[j], /*samt,*/ gcc[thid], map[thid], len[thid], e_lm);
    }, nwrks);

    et2 = tmn::elapsedTime(1);
    Rcpp::Rcout << " -- STOP NORMALISATION: " << et2 << " ms." << std::endl;

    double max_sc=0, min_sc=0;
    for(unsigned i=0; i<numRows; ++i) {
        if(miniEd[i].score>max_sc)
            max_sc = miniEd[i].score;
        if(min_sc==0 || miniEd[i].score<min_sc)
            min_sc = miniEd[i].score;
    }
    Rcpp::Rcout << " \nMax Score: " << max_sc << "; Min Score: " << min_sc << std::endl;

    double ead=0.0, psum=0.0;
    for(unsigned i=0; i<numRows; ++i) {
        miniEd[i].prob = (miniEd[i].score - min_sc) / (max_sc - min_sc);

#ifdef TEST_EDGES
        miniEd[i].prob = miniEd[i].prob > e_lm ? miniEd[i].prob : -0.01;
#endif

        psum += miniEd[i].prob;
    }

    long probSum = std::round(psum);
    ead = ((double) 2 / (double) numnodes) * psum; // expected average degree
    Rcpp::Rcout << " Cumulative Edges Probabilty: " << probSum
                << ";\n Expected Average Degree: " << ead << std::endl;

    for(unsigned i=0; i<sam_size; ++i)
        ::operator delete(samt->at(i));
    samt->~vector();
    samt->clear();

    delete f;
    return wrap( miniEd );
}

// [[Rcpp::export]]
SEXP extendGenesCoords(std::string genesfile, std::string fragsfile, int nwrks) {
    std::vector<Fragment> frags;
    std::vector<Gene> gn;
    std::unordered_map<size_t, size_t> frag_qt;
    std::unordered_map<size_t, size_t> frag_pos;

    Parser pr;
    pr.miniParseGenesFile(genesfile, gn);
    pr.miniParseFragmentsFile(fragsfile, frags);
    size_t genes_size = gn.size();
    size_t frags_size = frags.size();

    Rcpp::IntegerVector starts(genes_size);
    Rcpp::IntegerVector stops(genes_size);

    size_t qt=0, i=0, chr=0, chunk=0;
    frag_qt.reserve(24);  // chr 1 - 22 + x, y
    frag_pos.reserve(24); // chr 1 - 22 + x, y

    ParallelFor ffpf(nwrks, false /*spin-wait*/, false /*spin-barrier*/);

    while(i < frags.size()) {
        qt = 0;
        chr = frags[i].getHS();
        frag_pos.emplace( chr, i );
        while( i < frags.size() && frags[i].getHS() == chr ) {
            ++i;
            ++qt;
        }
        frag_qt.emplace( chr, qt );
    }

    ffpf.parallel_for(0, genes_size, 1, chunk, [&] (const unsigned s) {
        uint_64 start=0, stop=0;
        uint_64 p = frag_pos[gn[s].getHS()];
        uint_64 q = frag_qt[gn[s].getHS()];
        uint_64 stp = p + q;
        bool found = false;
        while(p < stp && !found) {
            if( ((frags[p].getStart() <= gn[s].getStart()) &&
                (frags[p].getStop() >= gn[s].getStart()))
            ) {
                start = frags[p].getStart();
                while(p < stp) {
                    if( ( (frags[p].getStop() >= gn[s].getStop()) &&
                        (frags[p].getStart() <= gn[s].getStop()))
                    ) {
                        stop = frags[p].getStop();
                        found = true;
                        break;
                    }
                    ++p;
                }
            }
            ++p;
        }
        if (start != 0) {
            starts[s] = start; // gn[s].setExtendedStart(start);
            stops[s]  = stop;  // gn[s].setExtendedStop(stop);
        }
    }, nwrks);

    freeVector(gn);
    freeVector(frags);
    malloc_trim(0);

    return Rcpp::DataFrame::create(Rcpp::Named("starts") = starts,
                                   Rcpp::Named("stops") = stops );
}


// [[Rcpp::export]]
SEXP readFragmentFile(std::string file) {
    Parser p;
    std::vector<Fragment> fn;
    p.parseFragmentsFile(file, fn);

    return wrap(fn);
}

// [[Rcpp::export]]
SEXP readGenesFile(std::string file) {
    Parser p;
    std::vector<Gene> fn;
    p.parseGenesFile(file, fn);

    return wrap(fn);
}

// [[Rcpp::export]]
SEXP readSamFile(std::string file) {
    Parser p;
    std::vector<SamData> fn;
    p.scanSamFile(file, fn);

    return wrap(fn);
}
