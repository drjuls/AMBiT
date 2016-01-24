#ifndef ONE_ELECTRON_OPERATOR_H
#define ONE_ELECTRON_OPERATOR_H

#include "HartreeFock/SpinorOperator.h"
#include "Basis/OrbitalManager.h"
#include "Configuration/ElectronInfo.h"

/** OneElectronIntegrals takes a SpinorMatrixElement and maps it to an electron operator for
    ManyBodyOperator, i.e. it implements GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2),
    adding the angular dependence on m1 and m2.
    Reduced matrix elements are stored using CalculateOneElectronIntegrals.
    The template parameter IsHermitianZeroOperator allows for two optimisations when the operator is
    Hermitian and K == 0:
    - If K == 0 it stores the matrix elements rather than reduced matrix elements,
      which saves multiplication by sqrt(2J+1)
    - It stores fewer matrix elements by assuming the Hermitian symmetry.
 */
template <bool IsHermitianZeroOperator = false>
class OneElectronIntegrals
{
public:
    OneElectronIntegrals(pOrbitalManagerConst pOrbitals, pSpinorMatrixElementConst pOperator):
        orbitals(pOrbitals), op(pOperator)
    {
        num_orbitals = orbitals->size();
    }

    /** Calculate one-electron integrals and return number of integrals stored.
        If check_size_only is true, the integrals are not calculated, but the storage size is returned.
     */
    virtual unsigned int CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only = false);

    /** Clear all integrals. */
    void clear() { integrals.clear(); }

    /** Number of stored integrals. */
    unsigned int size() const { return integrals.size(); }

    /** Return matrix element of stretched states orb1 and orb2. */
    double GetMatrixElement(const OrbitalInfo& orb1, const OrbitalInfo& orb2) const
    {
        double matrix_element = 0.0;

        if(IsHermitianZeroOperator)
        {
            if(orb1.Kappa() == orb2.Kappa())
            {
                unsigned int i1 = orbitals->state_index.at(orb1);
                unsigned int i2 = orbitals->state_index.at(orb2);

                matrix_element = integrals.at(GetKey(i1, i2));
            }
        }
        else if(op->IsNonZero(orb1, orb2))
        {
            MathConstant* math = MathConstant::Instance();
            matrix_element = math->Electron3j(orb2.TwoJ(), orb1.TwoJ(), op->GetK(), orb2.TwoJ(), -orb1.TwoJ());

            unsigned int i1 = orbitals->state_index.at(orb1);
            unsigned int i2 = orbitals->state_index.at(orb2);

            matrix_element *= integrals.at(GetKey(i1, i2));
        }

        return matrix_element;
    }

    /** Returns < e1 | O | e2 > including angular factors. */
    double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const
    {
        double matrix_element = 0.0;

        if(op->IsNonZero(e1, e2))
        {
            if(IsHermitianZeroOperator) // Optimization for K == 0
            {
                if(e1.TwoM() == e2.TwoM())
                {
                    unsigned int i1 = orbitals->state_index.at(e1);
                    unsigned int i2 = orbitals->state_index.at(e2);

                    matrix_element = integrals.at(GetKey(i1, i2));
                }
            }
            else
            {   MathConstant* math = MathConstant::Instance();
                matrix_element = math->Electron3j(e2.TwoJ(), e1.TwoJ(), op->GetK(), e2.TwoM(), -e1.TwoM())
                                 * math->minus_one_to_the_power((e1.TwoJ()-e1.TwoM())/2);

                if(matrix_element)
                {
                    unsigned int i1 = orbitals->state_index.at(e1);
                    unsigned int i2 = orbitals->state_index.at(e2);

                    matrix_element *= integrals.at(GetKey(i1, i2));
                }
            }
        }

        return matrix_element;
    }

    /** Return reduced matrix element. */
    double GetReducedMatrixElement(const OrbitalInfo& orb1, const OrbitalInfo& orb2) const
    {
        double matrix_element = 0.0;

        if(op->IsNonZero(orb1, orb2))
        {
            unsigned int i1 = orbitals->state_index.at(orb1);
            unsigned int i2 = orbitals->state_index.at(orb2);

            matrix_element = integrals.at(GetKey(i1, i2));

            if(IsHermitianZeroOperator)
            {
                matrix_element *= sqrt(orb1.TwoJ() + 1);
            }
        }

        return matrix_element;
    }

    /** Read integrals, adding to existing keys or creating new ones. */
    void Read(const std::string& filename);
    void Write(const std::string& filename) const;

    /** Get the spinor operator. */
    pSpinorMatrixElementConst GetOperator() const { return op; }

protected:
    unsigned int GetKey(unsigned int i1, unsigned int i2) const
    {
        if(!IsHermitianZeroOperator || (i1 <= i2))
            return i1 * num_orbitals + i2;
        else
            return i2 * num_orbitals + i1;
    }

    pSpinorMatrixElementConst op;
    pOrbitalManagerConst orbitals;
    unsigned int num_orbitals;

    /** Integrals stores the reduced matrix elements or, if(IsHermitianZeroOperator), the normal matrix element. */
    std::map<unsigned int, double> integrals;
};

typedef OneElectronIntegrals<true> HFIntegrals;
typedef std::shared_ptr<HFIntegrals> pHFIntegrals;
typedef std::shared_ptr<const HFIntegrals> pHFIntegralsConst;

typedef OneElectronIntegrals<false> TransitionIntegrals;
typedef std::shared_ptr<TransitionIntegrals> pTransitionIntegrals;
typedef std::shared_ptr<const TransitionIntegrals> pTransitionIntegralsConst;

// Template functions

template <bool IsHermitianZeroOperator>
unsigned int OneElectronIntegrals<IsHermitianZeroOperator>::CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only)
{
    unsigned int i1, i2;
    pOrbitalConst s1, s2;

    pSpinorOperatorConst spinor_op(std::dynamic_pointer_cast<const SpinorOperator>(op));

    if(spinor_op && !check_size_only)
    {
        int TwoK = 2 * op->GetK();
        MathConstant* math = MathConstant::Instance();
        pIntegrator integrator = op->GetIntegrator();

        auto it_1 = orbital_map_1->begin();
        while(it_1 != orbital_map_1->end())
        {
            i1 = orbitals->state_index.at(it_1->first);
            s1 = it_1->second;

            // Apply for all non-zero kappa_2
            int twoj_min = mmax(1, s1->TwoJ() - TwoK);
            int twoj_max = s1->TwoJ() + TwoK;

            std::vector<const SpinorFunction> op_applied;
            op_applied.reserve((twoj_max - twoj_min)/2 + 1);
            for(int twoj2 = twoj_min; twoj2 <= twoj_max; twoj2 += 2)
            {
                int kappa2 = math->convert_to_kappa(twoj2, s1->GetParity() * op->GetParity());

                if(IsHermitianZeroOperator)
                    op_applied.push_back(spinor_op->ApplyTo(*s1, kappa2));
                else
                    op_applied.push_back(spinor_op->ReducedApplyTo(*s1, kappa2));
            }

            // Take overlaps to get matrix elements
            auto it_2 = orbital_map_2->begin();
            while(it_2 != orbital_map_2->end())
            {
                s2 = it_2->second;

                if(op->IsNonZero(*s1, *s2))
                {
                    i2 = orbitals->state_index.at(it_2->first);

                    unsigned int key = GetKey(i1, i2);

                    if(!integrals.count(key))
                    {
                        const SpinorFunction& rhs = op_applied[(s2->TwoJ() - twoj_min)/2];
                        integrals[key] = integrator->GetInnerProduct(*s2, rhs);
                    }
                }
                it_2++;
            }

            it_1++;
        }
    }
    else
    {   std::set<unsigned int> found_keys;   // For check_size_only

        auto it_1 = orbital_map_1->begin();
        while(it_1 != orbital_map_1->end())
        {
            i1 = orbitals->state_index.at(it_1->first);
            s1 = it_1->second;

            auto it_2 = orbital_map_2->begin();
            while(it_2 != orbital_map_2->end())
            {
                s2 = it_2->second;

                if(op->IsNonZero(*s1, *s2))
                {
                    i2 = orbitals->state_index.at(it_2->first);

                    unsigned int key = GetKey(i1, i2);
                    if(check_size_only)
                        found_keys.insert(key);
                    else if(!integrals.count(key))
                    {
                        if(IsHermitianZeroOperator)
                            integrals[key] = op->GetMatrixElement(*s1, *s2);
                        else
                            integrals[key] = op->GetReducedMatrixElement(*s1, *s2);
                    }
                }
                it_2++;
            }
            it_1++;
        }

        if(check_size_only)
            return found_keys.size();
    }

    return integrals.size();
}

template <bool IsHermitianZeroOperator>
void OneElectronIntegrals<IsHermitianZeroOperator>::Read(const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
    {   *errstream << "OneElectronIntegrals::Read: file " << filename << " not found." << std::endl;
        exit(1);
    }

    OrbitalIndex old_state_index;
    ReadOrbitalIndexes(old_state_index, fp);
    unsigned int old_num_states = old_state_index.size();

    unsigned int num_integrals;
    unsigned int old_key;
    double value;

    fread(&num_integrals, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < num_integrals; i++)
    {
        fread(&old_key, sizeof(unsigned int), 1, fp);
        fread(&value, sizeof(double), 1, fp);

        unsigned int i1 = old_key/old_num_states;
        unsigned int i2 = old_key - i1 * old_num_states;
        unsigned int new_key = GetKey(i1, i2);

        auto it = integrals.find(new_key);
        if(it == integrals.end())
            integrals[new_key] = value;
        else
            it->second += value;
    }

    fclose(fp);
}

template <bool IsHermitianZeroOperator>
void OneElectronIntegrals<IsHermitianZeroOperator>::Write(const std::string& filename) const
{
    FILE* fp = fopen(filename.c_str(), "wb");

    // Write state index
    WriteOrbitalIndexes(orbitals->state_index, fp);

    unsigned int num_integrals = size();
    fwrite(&num_integrals, sizeof(unsigned int), 1, fp);

    for(auto& pair: integrals)
    {
        fwrite(&pair.first, sizeof(unsigned int), 1, fp);
        fwrite(&pair.second, sizeof(double), 1, fp);
    }

    fclose(fp);
}

#endif
