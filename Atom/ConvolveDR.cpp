#include "ConvolveDR.h"

int main(int argc, char* argv[])
{
    Convolver C;
    C.ReadCrossSections("dr.txt");
    C.BinCrossSections(0.001, 1.0);
    C.Convolve(0.01, 0.00015);
    
    return 0;
}

void Convolver::ReadCrossSections(const std::string filename)
{
    FILE* fp = fopen(filename.c_str(), "rt");
    if(fp)
        IntegratedDRCrossSections.clear();
    else
    {   std::cerr << "Unable to open file: " << filename << std::endl;
        exit(1);
    }
    
    char buffer[100];

    double energy, cross_section;
    while(fgets(buffer, 100, fp))
    {
        sscanf(buffer, "%le%le", &energy, &cross_section);

        // Just in case of degeneracy
        if(IntegratedDRCrossSections.find(energy) != IntegratedDRCrossSections.end())
            cross_section += IntegratedDRCrossSections[energy];
        // Cross sections in Mb (1Mb = 1.e-18 cm^2)
        IntegratedDRCrossSections[energy] = cross_section/100.;
    }

    fclose(fp);
}

void Convolver::BinCrossSections(double energy_width, double energy_max)
{
    num_bins = energy_max/energy_width;
    if(num_bins * energy_width < energy_max)
        num_bins++;

    bin_width = energy_width;
    if(bin_cross_sections)
        delete[] bin_cross_sections;

    bin_cross_sections = new double[num_bins];
    memset(bin_cross_sections, 0, num_bins);
    
    unsigned int i = 0;
    std::map<double, double>::const_iterator it = IntegratedDRCrossSections.begin();
    while((i < num_bins) && (it != IntegratedDRCrossSections.end()))
    {
        while((it != IntegratedDRCrossSections.end()) && (it->first < (i+1) * bin_width))
        {   bin_cross_sections[i] += it->second;
            it++;
        }

        i++;
    }
}

void Convolver::Convolve(double temp_perp, double temp_parallel)
{
    FILE* fp = fopen("drplot.txt", "wt");

    unsigned int num_points = num_bins * 5;
    double E_spacing = bin_width / 5.;

    double E = 0.0;
    fprintf(fp, "%12.6e\t%12.6e\n", E, 0.0);

    for(unsigned int i = 0; i < num_points; i++)
    {
        E = E + E_spacing;

        double total = 0.0;
        for(unsigned int m = 0; m < num_bins; m++)
            if(bin_cross_sections[m])
            {
                double Em = bin_width * m;
                double vm_over_2 = 2.966e-11 * sqrt(Em);
                double A = bin_cross_sections[m] * vm_over_2 / temp_perp
                           / sqrt(1.0 - temp_parallel/temp_perp);

                double B = exp(E/(temp_perp - temp_parallel) - Em/temp_perp);

                double z1 = sqrt(Em * (temp_perp - temp_parallel)/(temp_perp * temp_parallel));
                double z2 = sqrt(E * temp_perp/(temp_parallel * (temp_perp - temp_parallel)));                
                double C = ErrorFunction(z1 - z2) + ErrorFunction(z1 + z2);
                
                total = total + A*B*C;
                //std::cout << A << ", " << B << ", " << C << std::endl;
            }
        
        fprintf(fp, "%12.6e\t%12.6e\n", E, total);
    }
    
    fclose(fp);
}

double Convolver::ErrorFunction(double x)
{
    double a[5];
    a[0] =  0.254829592;
    a[1] = -0.284496736;
    a[2] =  1.421413741;
    a[3] = -1.453152027;
    a[4] =  1.061405429;

    double t = 1.0/(1.0 +  0.3275911*fabs(x));

    double sum = 0.0;
    double t_power = 1.;
    for(unsigned int i = 0; i < 5; i++)
    {   t_power *= t;
        sum += a[i] * t_power;
    }

    sum = 1.0 - sum * exp(-x*x);
    
    if(x < 0.0)
        sum = -sum;
    
    return sum;
}
