#include "Include.h"
#include "NonRelInfo.h"
#include "Universal/MathConstant.h"

std::string NonRelInfo::Name() const
{
    char buffer[20];
    sprintf(buffer, "%d", pqn);
    std::string ret(buffer);

    ret.append(1, MathConstant::Instance()->GetSpectroscopicNotation(L()));
    return ret;
}

void NonRelInfoSet::AddConfigs(const char* basis_def) 
{

    unsigned int pqn = 0;
    
    unsigned int p = 0;
    bool ReadSpectroscopicInput = true;
    while(basis_def[p]) {
        if(isdigit(basis_def[p])) 
        {
            if(ReadSpectroscopicInput) 
            {
                pqn = atoi(&basis_def[p]);
            }
            ReadSpectroscopicInput = false;
        } 
        else
        {
            for(int i = 0; i <= pqn; i++) 
            {
                if(i > MathConstant::Instance()->GetL(basis_def[p]))
                {
                    insert(NonRelInfo(i, MathConstant::Instance()->GetL(basis_def[p])));
                }
            }
            ReadSpectroscopicInput = true;
        }

        p++;
    }
}

void NonRelInfoSet::EraseConfigs(const char* basis_def) 
{
    unsigned int pqn = 0;

    unsigned int p = 0;
    bool ReadSpectroscopicInput = true;
    while(basis_def[p]) {
        if(isdigit(basis_def[p])) 
        {
            if(ReadSpectroscopicInput)
            {
                pqn = atoi(&basis_def[p]);
            }
            ReadSpectroscopicInput = false;
        } else 
        {
            for(int i = 0; i <= pqn; i++) 
            {
                NonRelInfoSet::iterator it = find(NonRelInfo(i, MathConstant::Instance()->GetSpectroscopicNotation(basis_def[p])));
                if(it != end()){
                    erase(it);
                }
                it++;
            }
            ReadSpectroscopicInput = true;
        }

        p++;
    }
}
