#ifndef CRPROPA_PARTICLECOLLECTOR_H
#define CRPROPA_PARTICLECOLLECTOR_H
#include <vector>
#include <string>

#include "crpropa/Module.h"
#include "crpropa/ModuleList.h"

namespace crpropa {

/**
 @class ParticleCollector
 @brief A helper ouput mechanism to directly transfer candidates to Python
 */
class ParticleCollector: public Module {
protected:
        typedef std::vector<ref_ptr<Candidate> > tContainer;
        mutable tContainer container;
        std::size_t nBuffer;
	bool clone;
	bool recursive;

public:
        ParticleCollector();
        ParticleCollector(const std::size_t nBuffer);
        ParticleCollector(const std::size_t nBuffer, const bool clone);
        ParticleCollector(const std::size_t nBuffer, const bool clone, const bool recursive);
        ~ParticleCollector();

        void process(Candidate *candidate) const;
	void reprocess(Module *action) const;
	void dump(const std::string &filename) const;
	void load(const std::string &filename);

        std::size_t getCount() const;
	ref_ptr<Candidate> operator[](const std::size_t i) const;
        void clearContainer();
        
	std::string getDescription() const;
	std::vector<ref_ptr<Candidate> > getAll() const;
	void setClone(bool b);

	/** iterator goodies */
        typedef tContainer::iterator iterator;
        typedef tContainer::const_iterator const_iterator;
        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;

	/**
	 Retrieves the trajectory of a detected particle
	 Procedure: take the initial state of the particle, re-run the ModuleList for that particle, save trajectory
	*/
	ref_ptr<ParticleCollector> getTrajectory(ModuleList* mlist, std::size_t i) const;
	ref_ptr<ParticleCollector> getTrajectory(ref_ptr<ModuleList> mlist, std::size_t i) const;
};

} // namespace crpropa

#endif // CRPROPA_PARTICLECOLLECTOR_H
