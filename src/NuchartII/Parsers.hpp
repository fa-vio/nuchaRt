
#ifndef PARSERS_HPP_
#define PARSERS_HPP_

#include "common.hpp"
#include "MemoryMapped.h"
#include "SamData.hpp"
#include "Gene.hpp"
#include "Fragment.hpp"
#include "PrintOut.hpp"

#include "FileIO.hpp"

#include <Rcpp.h>

struct Tfeat {
	std::string chr, st, end, len, gcc, map;
};

struct Tgene {
	std::string chr, strand, start, stop, symbol, id;
};

struct Tfrag {
	std::string chr, start, end, frgNum; //, r1, r5, r3;
};

struct Intr {
	uint_64 pos;
	size_t hs;
};


class Parser {
private:
	std::string genefile;
	std::string fragfile;
	std::string samfile;

public:

	Parser() { }
	Parser(std::string genes, std::string frags, std::string sam) :
		genefile(genes), fragfile(frags), samfile(sam) { }

	void parseFiles(std::vector<Gene> &gn, std::vector<Fragment> &fg, std::vector<SamData> &sm) {
		scanSamFile(samfile, sm);
		parseFragmentsFile(fragfile, fg);
		parseGenesFile(genefile, gn);
	}

	void parseIntersecFile(std::vector<Intr>& locations, std::string infile) {
		std::ifstream fs;
		std::string line, tmp, chr, pt;
		std::stringstream iss;
		int lns = 0;
		Intr point;

		if( exist_f(infile.c_str()) ) {
		    fs.open(infile.c_str(), std::ifstream::in);
		    if (!fs.good()) 
		        Rcpp::stop("Error opening intersections file\n");
		     
		    lns = countLines(infile);
		    locations.reserve(lns);
		} else {
		    Rcpp::stop(" [Parse Integration file] Error: File does not exists.\nMake sure the path to the file is correct.\n");
		}
		
		fs.ignore(10000, '\n');

		errno=0;
		char *end;
		while(fs) {
		    if ( std::getline(fs, chr, '\t') &&
                    std::getline(fs, pt, '\t') &&
                    std::getline(fs, tmp)
		    ) {
		        point.pos = strtol(pt.c_str(), &end, 10);
		        if(errno != 0) {
		            std::cerr << " [Parse Integration file] Position - Conversion error:\n "
                        << strerror(errno) << std::endl;
		        } else if (*end) {
		            std::cerr << " [Parse Integration file] Position - converted partially.\n"
                        << "Non-convertible part: "	<< end << std::endl;
		        }
		        point.hs = DJBHash(chr);
		        locations.push_back(point);
		    }
		}
		fs.close();
	}

	void buildCisTransMatrix(std::vector<CisTrans>& csts, std::string infile) {
	    std::ifstream fs;
	    std::string deg, tnid, cs, ts;
	    CisTrans cst;
	    
	    if( exist_f(infile.c_str()) ) {
	        fs.open(infile.c_str(), std::ifstream::in);
	        if (!fs.good())
	            Rcpp::stop("Error opening Cis-Trans file\n");

	        int lns = countLines(infile);
	        csts.reserve(lns);
	    } else {
	        Rcpp::stop("Error: Cis-Trans File does not exists.\nMake sure the path to the file is correct.\n");
	    }
	    
	    fs.ignore( 10000, '\n');
	    
	    errno=0;
	    while(fs) {
	        if ( std::getline(fs, tnid, '\t') &&
              std::getline(fs, deg, '\t') &&
              std::getline(fs, cs, '\t') &&
              std::getline(fs, ts)
	        ) {
	            try {
	                cst.nid = std::stol(tnid); //	, (char**)NULL, 10);
	            } catch (const std::invalid_argument& ia) {
	                std::cerr << "NodeId - Invalid argument: " << ia.what() << '\n';
	                continue;
	            }
	            
	            try {
	                cst.deg   = std::stol(deg); //, (char**)NULL, 10);
	            } catch (const std::invalid_argument& ia) {
	                std::cerr << "Degree - Invalid argument: " << ia.what() << '\n';
	                continue;
	            }
	            
	            try {
	                cst.cis   = std::stol(cs); //, (char**)NULL, 10);
	            } catch (const std::invalid_argument& ia) {
	                std::cerr << "Cis - Invalid argument: " << ia.what() << '\n';
	                continue;
	            }
	            
	            try {
	                cst.trans   = std::stol(ts); //, (char**)NULL, 10);
	            } catch (const std::invalid_argument& ia) {
	                std::cerr << "Trans - Invalid argument: " << ia.what() << '\n';
	                continue;
	            }
	            
	            csts.push_back(cst);
	        }
	    }
	    fs.close();
	}

	// parse 'gf_full.txt' file, if no other specified
	void parseFeaturesFile(std::vector<Features*>& uniques, std::string pt="./extdata/human/GF/gf_full.txt") {
	    std::ifstream fs;
	    std::string s_path(pt), line;
	    Features gf;
	    Tfeat tg;
	    
	    if( exist_f(s_path.c_str()) ) {
	        fs.open(s_path.c_str(), std::ifstream::in);
	        if (!fs.good())
	            Rcpp::stop("Error opening features file\n");

	        uniques.reserve( countLines(s_path) );
	    } /*else {
        int r = getGFfile("full");
	    if(r>0) {
	    fs.open(s_path.c_str(), std::ifstream::in);
	    if (!fs.good()) {
	    std::cerr << "Error: could not open \'"<< s_path << "\'" << std::endl;
	    exit(EXIT_FAILURE);
	    }
	    uniques.reserve( countLines(s_path) );
	    }*/
	    else {
	        Rcpp::stop(" [Parse Features file] ERROR - Features file does not exist\n");
	        //exit(EXIT_FAILURE);
	    }
	    
	    errno=0;
	    while(fs) {
	        if ( std::getline(fs, tg.chr, '\t') 	&&
                    std::getline(fs, tg.st, '\t') 	&&
                    std::getline(fs, tg.end, '\t')  &&
                    std::getline(fs, tg.len, '\t')  &&
                    std::getline(fs, tg.gcc, '\t')	&&
                    std::getline(fs, tg.map)
	        ) {
	            try {
	                gf.start = std::stod(tg.st); //.c_str(), (char**)NULL);
	            } catch (const std::invalid_argument& ia) {
	                std::cerr << " [Parse Features file] start - Invalid argument: " << ia.what() << '\n';
	                continue;
	            }
	            
	            try {
	                gf.stop = std::stod(tg.end); //.c_str(), (char**)NULL);
	            } catch (const std::invalid_argument& ia) {
	                std::cerr << " [Parse Features file] stop - Invalid argument: " << ia.what() << '\n';
	                continue;
	            }
	            
	            try {
	                gf.len = std::stol(tg.len); //.c_str(), (char**)NULL);
	            } catch (const std::invalid_argument& ia) {
	                std::cerr << " [Parse Features file] len - Invalid argument: " << ia.what() << '\n';
	                continue;
	            }
	            
	            if(tg.gcc.compare("NaN") == 0)
	                gf.gcc = 0.0;
	            else {
	                try {
	                    gf.gcc = std::stod(tg.gcc); //.c_str(), (char**)NULL);
	                } catch (const std::invalid_argument& ia) {
	                    std::cerr << " [Parse Features file] gcc - Invalid argument: " << ia.what() << '\n';
	                    continue;
	                }
	            }
	            
	            if(tg.map.compare("NaN") == 0)
	                gf.map = 0.0;
	            else {
	                try {
	                    gf.map = std::stod(tg.map); //.c_str(), (char**)NULL);
	                } catch (const std::invalid_argument& ia) {
	                    std::cerr << " [Parse Features file] map - Invalid argument: " << ia.what() << '\n';
	                    continue;
	                }
	            }
	            
	            gf.hs = DJBHash(tg.chr);
	            
	            uniques.push_back( new Features(gf) );
	        }
	    }
	    fs.close();
	}
    
    // ----------------------- Expressions ------------------------------------
    void parseExpressionFile(std::vector<Expression*>& ev, std::string infile="src/PnuChart/extdata/human/expression.txt") {
        if( !exist_f(infile.c_str()) )
            Rcpp::stop(" [Parse Expression file] Error: expressions File does not exists.\n");
        
        std::ifstream fs(infile.c_str(), std::ifstream::in);
        if (!fs.good())
            Rcpp::stop("Error: could not open expressions file");
        
        ev.reserve( countLines(infile) );
        fs.ignore(10000, '\n');
        
        Expression fm;
        std::string c1, c2, c3, c4;
        
        errno=0;
        while(fs) {
            if ( std::getline(fs, c1, ' ') &&
                 std::getline(fs, c2, ' ') &&
                 std::getline(fs, c3, ' ') &&
                 std::getline(fs, c4)
            ) {
                long ic1 = strtol(c1.c_str(), (char**)NULL, 10);
                if( ic1 == 0) {
                    continue;
                }
                else fm.eid = strtol(c1.c_str(), (char**)NULL, 10);
                
                if(c2.compare("NA") == 0) {
                    continue;
                }
                
                if(c3.compare("NA") == 0) {
                    continue;
                }
                else fm.logfc = strtod(c3.c_str(), (char**)NULL);
                
                if(c4.compare("NA") == 0 ) {
                    continue;
                }
                else fm.pval = strtod(c4.c_str(), (char**)NULL);
                
                if(errno != 0)
                    printf(" [Parse Expression file] - Conversion error: %s\n", strerror(errno));
                
                fm.hs = DJBHash(c2);
                
                ev.push_back( new Expression(fm.hs, fm.eid, fm.logfc, fm.pval) );
            }
        }
        fs.close();
    }
    
    // TODO: sistemare, magari generalizzare per differenti BED files; confrontare con versione R
    // ----------------------- BindingSite ------------------------------------
    void parseBedFile(std::vector<BindingSite*>& beds) {
        std::ifstream fs;
        
        char curPath[512];
        
        if( !getcwd(curPath, sizeof(curPath)) ) {
            std::cout << "Features file not found." << std::endl;
            exit(EXIT_FAILURE);
        }
        
        std::string s_path(curPath);
        s_path.append("/NgrapH/src/PnuChart/extdata/BED/full.bed");
        
        std::string line, c1, c2, chrom;
        std::stringstream iss;
        
        if( exist_f(s_path.c_str()) ) {
            beds.reserve( countLines(s_path) );
            fs.open(s_path.c_str(), std::ifstream::in);
            if (!fs.good()) {
                std::cerr << "Error opening \'"<< s_path.c_str() << "\'" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        errno=0;
        char *end;
        while(std::getline(fs, line)) {
            if(line.at(0) != '#') { // header
                BindingSite gf;
                iss << line;
                
                if ( std::getline(iss, chrom, '\t') 	 &&
                     std::getline(iss, c1, '\t') 	 &&
                     std::getline(iss, c2, '\t') 	 &&
                     std::getline(iss, gf.seq, '\t')
                ) {
                    gf.st = strtol(c1.c_str(), &end, 10);
                    if(errno != 0) {
                        printf("BindingSite_start - Conversion error: %s\n", strerror(errno));
                    } else if (*end) {
                        printf("Converted partially: non-convertible part: %s\n", end);
                    }
                    gf.sp  = strtol(c2.c_str(), &end, 10);
                    if(errno != 0) {
                        printf("BindingSite_stop - Conversion error: %s\n", strerror(errno));
                    } else if (*end) {
                        printf("Converted partially: non-convertible part: %s\n", end);
                    }
                    
                    gf.hs = DJBHash(chrom);
                    beds.push_back( new BindingSite(gf.hs,gf.seq,gf.st, gf.sp) );
                }
                iss.clear(); iss.str(std::string()); line.clear(); chrom.clear();
                c1.clear(); c2.clear();
            }
        }
        fs.close();
    }

	// ----------------------- Sam ------------------------------------
	void parseOrderedSamFile(std::string infile, std::vector<SamData>& sm, bool skip_head=true) {
	    if( !exist_f(infile.c_str()) )
	        Rcpp::stop(" [Parse SAM file] Error: SAM File does not exists.\nMake sure the path to your SAM file is correct.\n");
	    
	    MemoryMapped data(infile, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
	    if (!data.isValid())
	        Rcpp::stop(" [Parse SAM file] Error opening SAM file.\n");
	    
	    //double elpsd = tmn::elapsedTime(0);
	    int lns = countLines(&data);
	    //elpsd = tmn::elapsedTime(1);
	    if(sm.capacity() == 0) sm.reserve( lns );
	    
	    std::cout << " [SAM file] parsing..." << std::endl;
	    
	    SamData sd;
	    std::istringstream iss;
	    std::string st1, st2, rn, chr1, chr2, seq1;
	    
	    char *buffer = (char*) data.getData();
	    size_t data_sz = data.size();
	    iss.rdbuf()->pubsetbuf(buffer, data_sz);
	    
	    if(skip_head)
	        iss.ignore( 10000, '\n' );
	    
	    errno=0;
	    char *end;
	    double elpsd = tmn::elapsedTime(0);
	    while(iss) {
	        if ( std::getline(iss, rn, '\t') 	 &&
              std::getline(iss, chr1, '\t') &&
              std::getline(iss, st1, '\t')  &&
              std::getline(iss, chr2, '\t') &&
              std::getline(iss, st2, '\t')  &&
              std::getline(iss, seq1)
	        ) {
	            sd.setChr1(chr1);
	            sd.setChr2(chr2);
	            sd.setSeq(seq1);
	            
	            sd.setStart1( strtol(st1.c_str(), &end, 10) );
	            if(errno != 0) {
	                std::cerr << " [Parse SAM file] start1 - Conversion error:\n" << strerror(errno) << std::endl;
	                errno=0;
	            } else if (*end) {
	                std::cerr << " [Parse SAM file] Converted partially: non-convertible part:\n" << end << std::endl;
	            }
	            
	            sd.setStart2( strtol(st2.c_str(), &end, 10) );
	            if(errno != 0) {
	                std::cerr << " [Parse SAM file] start2 - Conversion error:\n" << strerror(errno) << std::endl;
	                errno=0;
	            } else if (*end) {
	                std::cerr << " [Parse SAM file] Converted partially: non-convertible part:\n" << end << std::endl;
	            }
	            sd.setHS1( DJBHash(chr1) );
	            sd.setHS2( DJBHash(chr2) );
	            
	            sm.push_back( sd );
	        }
	    }
	    elpsd = tmn::elapsedTime(1);
	    Rcpp::Rcout << " [SAM file] parsing completed (" << elpsd << " ms)" << std::endl;
	    iss.clear(); iss.str(std::string());
	    data.close();
	    
	    for(unsigned s=0; s<sm.size(); ++s) {
	        sm[s].setId(s);
	        if( sm[s].getChr2().compare(0,1,"=") == 0 )
	            sm[s].setChr2( sm[s].getChr1() );
	    }
	    int rsz = ((sizeof(SamData)*sm.capacity())+24) >> 30;
	    
	    Rcpp::Rcout << " [SAM file] " << sm.size() << " paired-end reads | ";
	    if(rsz==0) {
	        rsz = ((sizeof(SamData)*sm.capacity())+24) >> 20;
	        std::cout << "Reads dataset: " << rsz << " MB" << std::endl;
	    } else
	        std::cout << "Reads dataset: " << rsz << " GB" << std::endl;
	}

    // parse FULL sam file (not ordered, all columns present)
	void parseSamFile(std::string infile, std::vector<SamData>& sm, bool skip_head=true) {
	    if( !exist_f(infile.c_str()) )
	        Rcpp::stop(" [Parse SAM file] Error: SAM File does not exists.\nMake sure the path to your SAM file is correct.\n");
	    
	    MemoryMapped data(infile, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
	    if (!data.isValid())
	        Rcpp::stop(" [Parse SAM file] Error opening SAM file.\n");
	    
	    //double elpsd = tmn::elapsedTime(0);
	    int lns = countLines(&data);
	    //elpsd = tmn::elapsedTime(1);
	    if(sm.capacity() == 0) sm.reserve( lns );
	    
	    SamData sd;
	    std::istringstream iss;
	    std::string line, st1, st2, chr1, chr2, seq1, entry, tmp;
	    tmp.reserve(128);
	    
	    char *buffer = (char*) data.getData();
	    size_t data_sz = data.size();
	    iss.rdbuf()->pubsetbuf(buffer, data_sz);
	    line.reserve(256);
	    
	    if(skip_head)
	        iss.ignore( 10000, '\n' );
	    
	    errno=0;
	    char *end;
	    double elpsd = tmn::elapsedTime(0);
	    while(iss) {
	        if ( std::getline(iss, entry, '\t') 	 &&
              std::getline(iss, tmp, '\t') 	 &&
              std::getline(iss, chr1, '\t') &&
              std::getline(iss, st1, '\t')  &&
              std::getline(iss, tmp, '\t') 	 &&
              std::getline(iss, tmp, '\t') 	 &&
              std::getline(iss, chr2, '\t') &&
              std::getline(iss, st2, '\t')  &&
              std::getline(iss, tmp, '\t') 	 &&
              std::getline(iss, seq1, '\t') &&
              std::getline(iss, tmp)
	        ) {
	            sd.setEntry(entry);
	            sd.setChr1(chr1);
	            sd.setChr2(chr2);
	            sd.setSeq(seq1);
	            
	            sd.setStart1( strtol(st1.c_str(), &end, 10) );
	            if(errno != 0) {
	                printf(" [Parse SAM file] start1 - Conversion error: %s\n", strerror(errno));
	                errno=0;
	            } else if (*end) {
	                printf(" [Parse SAM file] Converted partially: non-convertible part: %s\n", end);
	            }
	            sd.setStart2( strtol(st2.c_str(), &end, 10) );
	            if(errno != 0) {
	                printf(" [Parse SAM file] start2 - Conversion error: %s\n", strerror(errno));
	                errno=0;
	            } else if (*end) {
	                printf(" [Parse SAM file] Converted partially: non-convertible part: %s\n", end);
	            }
	            
	            sd.setHS1( DJBHash(chr1) );
	            sd.setHS2( DJBHash(chr2) );
	            
	            sm.push_back( sd );
	        }
	    }
	    elpsd = tmn::elapsedTime(1);
	    Rcpp::Rcout << " [SAM file] parsed in " << elpsd << " ms.\n";
	    data.close();
	    
	    elpsd = tmn::elapsedTime(0);
	    for(unsigned s=0; s<sm.size(); ++s) {
	        if( sm[s].getChr2().compare(0, 1, "=") == 0 )
	            sm[s].setChr2( sm[s].getChr1() );
	    }
	    elpsd = tmn::elapsedTime(1);
	    Rcpp::Rcout << " [SAM file] fixed cis contacts in " << elpsd << " ms.\n";
	    
	    elpsd = tmn::elapsedTime(0);
	    std::sort( sm.begin(), sm.end() );
	    elpsd = tmn::elapsedTime(1);
	    Rcpp::Rcout << " [SAM file] ordered in " << elpsd << " ms.\n";
	    
	    std::string ninf = printOrderedSam(infile, sm);
	    Rcpp::Rcout << " [SAM file] ordered SAM file written to " << ninf;
	    
	    Rcpp::Rcout << "\nRestart using ordered SAM file...\n" << std::endl;
	    sm.clear();
	    parseOrderedSamFile(ninf, sm);
	}

	// check if it is and ordered sam file
	void scanSamFile(std::string f, std::vector<SamData>& sm) {
		if( f.find("ord_", 0) == std::string::npos )    //f.compare(0, ord.length(), ord) == 0)
			parseSamFile(f, sm);
		else
			parseOrderedSamFile(f, sm);
	}

	// ----------------------- Fragments ------------------------------------

	void parseFragmentsFile(std::string infile, std::vector<Fragment>& fg, bool skip_head=true) {
	    if( !exist_f(infile.c_str()) )
	        Rcpp::stop(" [Parse Fragments file] Error: fragments File does not exists.\n");
	    
	    std::ifstream fs(infile, std::ifstream::in);
	    if (!fs.good()) {
	        Rcpp::stop(" [Parse Fragments file] Error: could not open fragments file.\n");
	        exit(EXIT_FAILURE);
	    }
	    
	    int lns = countLines(infile);
	    fg.reserve( lns );
	    
	    Fragment fm;
	    std::string cn, tmp, line;
	    Tfrag tf;
	    line.reserve(256);
	    
	    // skip two lines
	    if(skip_head)
	        for(unsigned i=0; i<2; ++i)
	            fs.ignore( 10000, '\n' );
	    
	    errno=0;
	    char *end;
	    while(fs) {
	        if ( std::getline(fs, cn, '\t') &&
              std::getline(fs, tf.start, '\t') &&
              std::getline(fs, tf.end, '\t') &&
              std::getline(fs, tf.frgNum, '\t') &&
              std::getline(fs, tmp)
	        ) {
	            fm.setChromosome(cn);
	            
	            fm.setStart( strtol(tf.start.c_str(), &end, 10) );
	            if(errno != 0) {
	                printf(" [Parse Fragments file] start - Conversion error: %s\n", strerror(errno));
	                errno=0;
	            } else if (*end)
	                printf(" [Parse Fragments file] start - Converted partially: non-convertible part: %s\n", end);
	            
	            fm.setStop( strtol(tf.end.c_str(), &end, 10) );
	            if(errno != 0) {
	                printf(" [Parse Fragments file] stop - Conversion error: %s\n", strerror(errno));
	                errno=0;
	            } else if (*end)
	                printf(" [Parse Fragments file] stop - Converted partially: non-convertible part: %s\n", end);
	            
	            fm.setFragNum( strtol(tf.frgNum.c_str(), &end, 10) );
	            if(errno != 0) {
	                printf(" [Parse Fragments file] fragNum - Conversion error: %s\n", strerror(errno));
	                errno=0;
	            } else if (*end)
	                printf(" [Parse Fragments file] fragNum - Converted partially: non-convertible part: %s\n", end);
	            
	            fm.setHS( DJBHash(cn) );
	            
	            fg.push_back( fm );
	        }
	    }
	    fs.close();
	    
	    fg.shrink_to_fit();
	    Rcpp::Rcout << " [Fragments] " << fg.size() << " Fragments | ";
	    Rcpp::Rcout << "Fragments dataset: " << (((sizeof(Fragment)*fg.capacity())+24) >> 20) << " MB" << std::endl;
	}
    
    void miniParseFragmentsFile(std::string infile, std::vector<Fragment>& fg, bool skip_head=true) {
        if( !exist_f(infile.c_str()) )
            Rcpp::stop(" [Parse Fragments file] Error: fragments File does not exists.\n");
        
        std::ifstream fs(infile, std::ifstream::in);
        if (!fs.good()) {
            Rcpp::stop(" [Parse Fragments file] Error: could not open fragments file.\n");
            exit(EXIT_FAILURE);
        }
        
        int lns = countLines(infile);
        fg.reserve( lns );
        
        Fragment fm;
        std::string cn, tmp, line;
        Tfrag tf;
        line.reserve(256);
        
        // skip two lines
        if(skip_head)
            for(unsigned i=0; i<2; ++i)
                fs.ignore( 10000, '\n' );
        
        errno=0;
        char *end;
        while(fs) {
            if ( std::getline(fs, cn, '\t') &&
                 std::getline(fs, tf.start, '\t') &&
                 std::getline(fs, tf.end, '\t') &&
                 std::getline(fs, tmp)
            ) {
                fm.setChromosome(cn);
                
                fm.setStart( strtol(tf.start.c_str(), &end, 10) );
                if(errno != 0) {
                    printf(" [Parse Fragments file] start - Conversion error: %s\n", strerror(errno));
                    errno=0;
                } else if (*end)
                    printf(" [Parse Fragments file] start - Converted partially: non-convertible part: %s\n", end);
                
                fm.setStop( strtol(tf.end.c_str(), &end, 10) );
                if(errno != 0) {
                    printf(" [Parse Fragments file] stop - Conversion error: %s\n", strerror(errno));
                    errno=0;
                } else if (*end)
                    printf(" [Parse Fragments file] stop - Converted partially: non-convertible part: %s\n", end);
                
                fm.setHS( DJBHash(cn) );
                
                fg.push_back( fm );
            }
        }
        fs.close();
        
        fg.shrink_to_fit();
    }

	// ----------------------- Genes ------------------------------------

	void parseGenesFile(std::string infile, std::vector<Gene>& gv, bool skip_head=true) {
	    if( !exist_f(infile.c_str()) )
	        Rcpp::stop(" [Parse Genes file] Error: genes File does not exists.\n");
	    
	    std::ifstream fs(infile, std::ifstream::in);
	    if (!fs.good())
	        Rcpp::stop("Error: could not open genes file.\n");
	    
	    int lns = countLines(infile);
	    gv.reserve( lns );
	    
	    Gene gn;
	    std::string line;
	    size_t ensg, erbs;
	    Tgene tg;
	    line.reserve(256);
	    
	    if(skip_head)
	        fs.ignore( 10000, '\n' );
	    
	    errno=0;
	    char *end;
	    while(fs) {
	        if ( std::getline(fs, tg.chr, '\t') &&
              std::getline(fs, line, '\t') &&
              std::getline(fs, tg.start, '\t') &&
              std::getline(fs, tg.stop, '\t') &&
              std::getline(fs, tg.symbol, '\t') &&
              std::getline(fs, tg.id)
	        ) {
	            if( (ensg = tg.id.find_first_of("G")) != std::string::npos )
	                tg.id = tg.id.substr(ensg+1);
	            
	            if( (erbs = tg.id.find_first_of("_")) != std::string::npos )
	                tg.id = tg.id.substr(erbs+1);
	            
	            if( (tg.symbol.compare("NA") == 0 ) || (tg.symbol.empty()) ||
                 (tg.symbol.compare("") == 0) || (tg.id.compare("NA") == 0) ||
                 (tg.id.compare("0") == 0) || (tg.chr.size() > 5) ) {
	                continue;
	            } else {
	                gn.setChromosomeName(tg.chr);
	                
	                gn.setStartCoordinate( strtol(tg.start.c_str(), &end, 10) );
	                if(errno != 0) {
	                    printf(" [Parse Genes file] start - Conversion error: %s\n", strerror(errno));
	                    errno=0;
	                } else if (*end) {
	                    printf(" [Parse Genes file] start - Converted partially: non-convertible part: %s\n", end);
	                }
	                
	                gn.setStopCoordinate( strtol(tg.stop.c_str(), &end, 10) );
	                if(errno != 0) {
	                    printf(" [Parse Genes file] stop - Conversion error: %s\n", strerror(errno));
	                    errno=0;
	                } else if (*end) {
	                    printf(" [Parse Genes file] stop - Converted partially: non-convertible part: %s\n", end);
	                }
	                
	                gn.setSymbol(tg.symbol);
	                
	                long d = strtol(tg.id.c_str(), &end, 10);
	                if(errno != 0) {
	                    printf(" [Parse Genes file] id - Conversion error: %s\n", strerror(errno));
	                    errno=0;
	                } else if (*end) {
	                    printf(" [Parse Genes file] id - Converted partially: non-convertible part: %s\n", end);
	                }
	                
	                gn.setEID( d );
	                gn.setHS( DJBHash(tg.chr) );
	                
	                gv.push_back(gn);
	            }
	        }
	    }
	    fs.close();
	    
	    std::sort(gv.begin(), gv.end());
	    gv.erase( std::unique(gv.begin(), gv.end()), gv.end() );
	    
	    Rcpp::Rcout << " [Genes] " << gv.size() << " Genes | ";
	    Rcpp::Rcout << "Genes dataset: " << ( ((sizeof(Gene)*gv.capacity())+24) >> 20 ) << " MB" << std::endl;
	    
	    for(size_t i=0; i<gv.size(); ++i)
	        gv[i].setPosition(i);
	    
	    gv.shrink_to_fit();
	}
    
    void miniParseGenesFile(std::string infile, std::vector<Gene>& gv, bool skip_head=true) {
        if( !exist_f(infile.c_str()) )
            Rcpp::stop(" [Parse Genes file] Error: genes File does not exists.\n");
        
        std::ifstream fs(infile, std::ifstream::in);
        if (!fs.good())
            Rcpp::stop("Error: could not open genes file.\n");
        
        int lns = countLines(infile);
        gv.reserve( lns );
        
        Gene gn;
        std::string line;
        size_t ensg, erbs;
        Tgene tg;
        line.reserve(256);
        
        if(skip_head)
            fs.ignore( 10000, '\n' );
        
        errno=0;
        char *end;
        while(fs) {
            if ( std::getline(fs, tg.chr, '\t') &&
                 std::getline(fs, line, '\t') &&
                 std::getline(fs, tg.start, '\t') &&
                 std::getline(fs, tg.stop, '\t') &&
                 std::getline(fs, tg.symbol, '\t') &&
                 std::getline(fs, tg.id)
            ) {
                if( (ensg = tg.id.find_first_of("G")) != std::string::npos )
                    tg.id = tg.id.substr(ensg+1);

                
                if( (erbs = tg.id.find_first_of("_")) != std::string::npos )
                    tg.id = tg.id.substr(erbs+1);
                
                if( (tg.symbol.compare("NA") == 0 ) || (tg.symbol.empty()) ||
                    (tg.symbol.compare("") == 0) || (tg.id.compare("NA") == 0) ||
                    (tg.id.compare("0") == 0) || (tg.chr.size() > 5) ) {
                    continue;
                } else {
                    gn.setChromosomeName(tg.chr);
                    
                    gn.setStartCoordinate( strtol(tg.start.c_str(), &end, 10) );
                    if(errno != 0) {
                        printf(" [Parse Genes file] start - Conversion error: %s\n", strerror(errno));
                        errno=0;
                    } else if (*end) {
                        printf(" [Parse Genes file] start - Converted partially: non-convertible part: %s\n", end);
                    }
                    
                    gn.setStopCoordinate( strtol(tg.stop.c_str(), &end, 10) );
                    if(errno != 0) {
                        printf(" [Parse Genes file] stop - Conversion error: %s\n", strerror(errno));
                        errno=0;
                    } else if (*end) {
                        printf(" [Parse Genes file] stop - Converted partially: non-convertible part: %s\n", end);
                    }
                    
                    gn.setSymbol(tg.symbol);
                    
                    long d = strtol(tg.id.c_str(), &end, 10);
                    if(errno != 0) {
                        printf(" [Parse Genes file] id - Conversion error: %s\n", strerror(errno));
                        errno=0;
                    } else if (*end) {
                        printf(" [Parse Genes file] id - Converted partially: non-convertible part: %s\n", end);
                    }
                    
                    gn.setEID( d );
                    gn.setHS( DJBHash(tg.chr) );
                    
                    gv.push_back(gn);
                }
            }
        }
        fs.close();
        
        std::sort(gv.begin(), gv.end());
        gv.erase( std::unique(gv.begin(), gv.end()), gv.end() );
        
        for(size_t i=0; i<gv.size(); ++i)
            gv[i].setPosition(i);
        
        gv.shrink_to_fit();
    }
};




#endif /* PARSERS_HPP_ */
