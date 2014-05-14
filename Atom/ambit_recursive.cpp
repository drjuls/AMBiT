#include "ambit.h"

void Ambit::Recursive()
{
/*
    int Z = fileInput("Z", 0);
    std::ofstream outfile;
    outfile.open("temp.out", std::ios_base::trunc);
    if(!outfile.is_open())
    {
        *errstream << "Error: Could not open file 'temp.out' for use in recursive building." << std::endl;
        exit(1);
    }
    std::string configstring = fileInput("HF/Configuration", "1s1:");
    outfile << "ID = " << fileInput("ID", "") << "N" << fileInput("N", 1) << std::endl << "Z = " << Z << std::endl << "N = " << fileInput("N", 1) << std::endl << "NumValenceElectrons = " << 1 << std::endl;
    outfile << "-c" << std::endl << "-d" << std::endl;
    outfile << "--recursive-build" << std::endl << std::endl;
    outfile << "NuclearRadius = " << fileInput("NuclearRadius", 0) << std::endl << "NuclearThickness = " << fileInput("NuclearThickness", 0) << std::endl << std::endl;
    outfile << "[Lattice]" << std::endl << "NumPoints = " << fileInput("Lattice/NumPoints", 1000) << std::endl;
    outfile << "StartPoint = " << fileInput("Lattice/StartPoint", 1.e-6) << std::endl << "EndPoint = " << fileInput("Lattice/EndPoint", 50.) << std::endl << std::endl;
    outfile << "[HF]" << std::endl << "Configuration = '" << configstring << "'" << std::endl << std::endl;
    outfile << "[Basis]" << std::endl << "--bspline-basis" << std::endl << "ValenceBasis = " << fileInput("Basis/ValenceBasis", "4spdf") << std::endl << std::endl;
    outfile << "BSpline/N = " << fileInput("Basis/BSpline/N", 40) << std::endl << "BSpline/K = " << fileInput("Basis/BSpline/K", 7) << std::endl;
    outfile << "BSpline/Rmax = " << fileInput("Basis/BSpline/Rmax", 50.) << std::endl << std::endl;
    outfile.close();
    for(int i = fileInput("N", 1) - 1; i < Z; i++)
    {
        GetPot tempInput("temp.out", "//", "\n", ",");
        Atom* A;
        A = new Atom(tempInput, Z, tempInput("N", 1), tempInput("ID", ""));
        A->Run();

        if(i + 2 <= Z)
        {
            configstring = A->GetNextConfigString();
            *outstream << "+ " <<A->GetIteratorToNextOrbitalToFill()->first.Name() << " = " <<  std::endl;
            outfile.open("temp.out", std::ios_base::trunc);
            outfile << "ID = " << fileInput("ID", "") << "N" << i + 2 << std::endl << "Z = " << Z << std::endl << "N = " << i + 2 << std::endl << "NumValenceElectrons = " << 1 << std::endl;
            outfile << "-c" << std::endl << "-d" << std::endl;
            outfile << "--recursive-build" << std::endl << std::endl;
            outfile << "NuclearRadius = " << fileInput("NuclearRadius", 0) << std::endl << "NuclearThickness = " << fileInput("NuclearThickness", 0) << std::endl << std::endl;
            outfile << "[Lattice]" << std::endl << "NumPoints = " << fileInput("Lattice/NumPoints", 1000) << std::endl;
            outfile << "StartPoint = " << fileInput("Lattice/StartPoint", 1.e-6) << std::endl << "EndPoint = " << fileInput("Lattice/EndPoint", 50.) << std::endl << std::endl;
            outfile << "[HF]" << std::endl << "Configuration = '" << configstring << "'" << std::endl << std::endl;
            outfile << "[Basis]" << std::endl << "--bspline-basis" << std::endl << "ValenceBasis = " << fileInput("Basis/ValenceBasis", "4spdf") << std::endl << std::endl;
            outfile << "BSpline/N = " << fileInput("Basis/BSpline/N", 40) << std::endl << "BSpline/K = " << fileInput("Basis/BSpline/K", 7) << std::endl;
            outfile << "BSpline/Rmax = " << fileInput("Basis/BSpline/Rmax", 50.) << std::endl << std::endl;
            outfile.close();
        }
        delete A;
    }

    *outstream << std::endl << std::endl << "Final configuration: " << configstring;
 */
}