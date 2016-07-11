
#ifndef FRAGMENT_HPP_
#define FRAGMENT_HPP_

#include "common.hpp"
#include "FragmentT.hpp"


/** 
 * This class describes a chromosome fragment
 * chrom:       chromosome name or number
 * start_p:     fragment start position
 * stop_p:      fragment stop position
 * frag_num:    fragment number (sort of Unique ID)
 * RE1_rag_num: restriction site number - MAYBE NOT NEEDED!
 * rs5, rs3:    restriction site 5 and 3 - MAYBE NOT NEEDED!
 * 
 */
class Fragment {
	friend class Gene;
	friend class Edge;

public:
	Fragment() :
		start_p(0), stop_p(0), frag_num(0), hs(0) { }
        
    Fragment(std::string chr, uint_64 st, uint_64 end, uint_64 fragNum, size_t h) :
    	start_p(st), stop_p(end), frag_num(fragNum), hs(h) {
        chrom_name = chr;        
    }

	/* Copy from another Fragment object */
	Fragment(const Fragment& fg) {
		*this = fg;
	}

	// move constructor
	Fragment(Fragment&& fg): chrom_name(std::move(fg.chrom_name)), start_p(std::move(fg.start_p)),
			stop_p(std::move(fg.stop_p)), frag_num(std::move(fg.frag_num)), hs(std::move(fg.hs)) { }

	Fragment& operator=(const Fragment& fg) = default;

	Fragment& operator=(Fragment&& fg) {
		if(this != &fg) {
			chrom_name = std::move(fg.chrom_name);
			start_p = std::move(fg.start_p);
			stop_p = std::move(fg.stop_p);
			frag_num = std::move(fg.frag_num);
			hs = std::move(fg.hs);
		}
		return *this;
	}


	/* Destructor */
	~Fragment() {}

	// getters
	inline std::string getChromName() 	const { return chrom_name; }
	inline uint_64 getStart() 			const { return start_p; }
	inline uint_64 getStop() 			const { return stop_p; }
	inline size_t getHS() 				const { return hs; }
	inline uint_64 getFragNum()         const { return frag_num; }

	// setters
	inline void setChromosome(std::string chr) 	{ chrom_name = chr; }
	inline void setStart(uint_64 st)			{ start_p = st; }
	inline void setStop(uint_64 sp)				{ stop_p = sp; }
	inline void setFragNum(uint_64 fn)			{ frag_num = fn; }
	inline void setHS(size_t h)					{ hs = h; }

	FragmentT* toFragmentT() {
		FragmentT *ft = new FragmentT(start_p, stop_p);
		return ft;
	}

	/*
	 * Friendly print a Fragment entry
	 */
	friend std::ostream& operator<<(std::ostream &out, const Fragment &fg) {
		out << "Chrom: "         << fg.chrom_name 	<< " | " <<
				"Start: "        << fg.start_p 	  	<< " | " <<
				"End: "          << fg.stop_p 	  	<< " | " <<
				"Frag.Num: "     << fg.frag_num   	<< " | ";

		return out;
	}

	bool operator<(const Fragment &fg) const {
		return (hs < fg.hs) || (hs == fg.hs && frag_num < fg.frag_num);
	}

	bool operator==(const Fragment& fg) const {
		return (hs == fg.getHS() && frag_num == fg.frag_num);
	}


private:
	std::string chrom_name;
	uint_64 start_p;
	uint_64 stop_p;
	uint_64 frag_num;
	size_t hs;
};

// to be used as a functor
struct CompareFragments : public Fragment {
	bool operator()(const Fragment &fg1, const Fragment &fg2) {
		if( fg1.getStart() < fg2.getStart() && (fg1.getHS() == fg2.getHS()) )
			return true;
		else if( fg1.getStart() <= fg2.getStart() && (fg1.getChromName().compare(fg2.getChromName()) < 0) )
			return true;
		else return false;
	}
};

#endif /* FRAGMENT_HPP_ */
