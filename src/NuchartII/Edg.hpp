
#ifndef EDG_HPP_
#define EDG_HPP_

#include "common.hpp"
#include "Vertex.hpp"

/**
 * Class representig a generic edge object within two Vertices.
 * An edge can have a weight and a probability.
 */

class Edg {
public:
	Edg() : v1(NULL), v2(NULL), weight(0.1), prob(0.0) { }
	Edg(Vertex *v1_, Vertex *v2_, double w, double p) : v1(v1_), v2(v2_), weight(w), prob(p) { }
	Edg(const Edg& e) : v1(e.getVertex1()), v2(e.getVertex2()), 
                        weight(e.getWeight()), prob(e.getProb()) { }
	Edg(Edg&& e) : v1(NULL), v2(NULL), weight(e.getWeight()), prob(e.getProb()) {
		v1 = e.getVertex1();
		v2 = e.getVertex2();
	}

	Edg& operator=(const Edg& e) { // shallow-copy is fine!
		if(this != &e) {
			if(v1) v1=NULL;
			if(v2) v2=NULL;

			v1 = e.getVertex1();
			v2 = e.getVertex2();
			weight = e.getWeight();
			prob = e.getProb();
		}
		return *this;
	}
	Edg& operator=(Edg&& e) {
		if(this != &e) {
			if(v1) v1 = NULL;
			if(v2) v2 = NULL;

			v1 = e.getVertex1();
			v2 = e.getVertex2();
			weight = e.getWeight();
			prob = e.getProb();
		}
		return *this;
	}

	virtual ~Edg() { };

	inline Vertex* getVertex1() const { return v1; }
	inline Vertex* getVertex2() const { return v2; }
	inline double getWeight()   const { return weight; }
	inline double getProb()   	const { return prob; }

	void setVertex1(Vertex *v)  { v1 = v; }
	void setVertex2(Vertex *v)  { v2 = v; }
	void setWeight(double w)    { weight = w; }
	void setProb(double p)		{ prob = p; }

protected:
	Vertex *v1;
	Vertex *v2;
	double weight, prob;

};


#endif /* EDG_HPP_ */
