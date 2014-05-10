#ifndef MULTIRUN_OPTIONS_H
#define MULTIRUN_OPTIONS_H

#include "Atom/GetPot"
#include "Include.h"

/** A subclass of GetPot, MultirunOptions masks input vectors that are designed for
    multiple runs. Indicate which keys have multiple runs with the "Multirun" variable.
    For example, to run the code with multiple values of nuclear inverse mass
    the user can include a vector in the input file
        ...
        Multirun = 'NuclearInverseMass, MBPT/Delta'
        NuclearInverseMass = '-0.001, 0.0, 0.001'
        MBPT/Delta = '0.62, 0.65, 0.68'
        ...
    which should run the code three times with the three values.
    MultirunOptions will mask this, so that, e.g., on the first run a call
        MultirunOptions("NuclearInverseMass", 0.0);
    will simply return -0.001.

    GetPot is inherited privately because polymorphism cannot be supported since GetPot functions are non-virtual.
 */
class MultirunOptions : private GetPot
{
public:
    inline MultirunOptions();
    inline MultirunOptions(const MultirunOptions&);
    inline MultirunOptions(const int argc_, char** argv_,
                           const char* FieldSeparator=0x0);
    inline MultirunOptions(const char* FileName,
                           const char* CommentStart=0x0, const char* CommentEnd=0x0,
                           const char* FieldSeparator=0x0);
    inline MultirunOptions(std::istream& InputStream,
                           const char* CommentStart=0x0, const char* CommentEnd=0x0,
                           const char* FieldSeparator=0x0);
    inline ~MultirunOptions() {}

    inline MultirunOptions& operator=(const MultirunOptions& other);

    // (*) absorbing contents of another MultirunOptions object
    inline void absorb(const MultirunOptions& That);

    inline int GetNumRuns() const;
    inline void SetRun(int run_index);

    /** Attempt to find a run where varying parameter is zero.
        If successful, return parameter name and run where it is zero.
        If unsuccessful, return null string and -1.
     */
    inline std::pair<std::string, int> FindZeroParameterRun() const;

    /** Print special conditions of current run separated by separator. */
    inline void PrintCurrentRunCondition(std::ostream& outstream, const std::string& separator) const;

public:
    using GetPot::next;
    using GetPot::next_nominus;
    using GetPot::operator[];
    using GetPot::search;
    using GetPot::size;
    using GetPot::vector_variable_size;

    // Overwrite scalar double variables to hide multirun vectors.
    using GetPot::operator();
    inline double operator()(const char* VarName, const double& Default) const;

protected:
    inline void ParseMultirun();

protected:
    // Multiple run parameters
    STRING_VECTOR multirun_keys;
    std::vector< std::vector<double> > multirun_values;

    int num_runs;
    int current_run_index;
};

inline MultirunOptions::MultirunOptions(): GetPot(), num_runs(1), current_run_index(0)
{}

inline MultirunOptions::MultirunOptions(const MultirunOptions& other): GetPot()
{
    MultirunOptions::operator=(other);
}

inline MultirunOptions::MultirunOptions(const int argc_, char** argv_, const char* FieldSeparator):
    GetPot(argc_, argv_, FieldSeparator), num_runs(1), current_run_index(0)
{   ParseMultirun();
}

inline MultirunOptions::MultirunOptions(const char* FileName, const char* CommentStart, const char* CommentEnd, const char* FieldSeparator):
    GetPot(FileName, CommentStart, CommentEnd, FieldSeparator), num_runs(1), current_run_index(0)
{   ParseMultirun();
}

inline
MultirunOptions::MultirunOptions(std::istream& InputStream,
               const char* CommentStart  /* = 0x0 */, const char* CommentEnd /* = 0x0 */,
               const char* FieldSeparator/* = 0x0 */):
    GetPot()
{
    // if specified -> overwrite default strings
    if( CommentStart )   _comment_start = std::string(CommentStart);
    if( CommentEnd )     _comment_end   = std::string(CommentEnd);
    if( FieldSeparator ) _field_separator = FieldSeparator;

    STRING_VECTOR _apriori_argv;
    // First argument is not parsed for
    //    variable assignments or nominuses.
    _apriori_argv.push_back(std::string("StreamInput"));
    
    STRING_VECTOR args = __read_in_stream(InputStream);
    _apriori_argv.insert(_apriori_argv.begin()+1, args.begin(), args.end());
    __parse_argument_vector(_apriori_argv);

    ParseMultirun();
}

inline MultirunOptions& MultirunOptions::operator=(const MultirunOptions& other)
{
    GetPot::operator=(other);
    multirun_keys = other.multirun_keys;
    multirun_values = other.multirun_values;
    num_runs = other.num_runs;
    current_run_index = other.current_run_index;

    return *this;
}

inline void MultirunOptions::absorb(const MultirunOptions& other)
{
    GetPot::absorb(other);
    if(other.num_runs)
    {
        if(multirun_keys.size() && num_runs != other.num_runs)
        {   *errstream << "MultirunOptions::absorb() cannot absorb another object with a different run size." << std::endl;
            exit(1);
        }

        num_runs = other.num_runs;
        // Should probably check that there are no key overlaps, but can't be bothered.
        multirun_keys.insert(multirun_keys.end(), other.multirun_keys.begin(), other.multirun_keys.end());
        multirun_values.insert(multirun_values.end(), other.multirun_values.begin(), other.multirun_values.end());
    }
}

inline double MultirunOptions::operator()(const char* VarName, const double& Default) const
{
    double ret = 0.0;
    int multirun_index = -1;

    for(int i = 0; i < multirun_keys.size(); i++)
    {   if(!strcmp(VarName, multirun_keys[i].c_str()))
        {   multirun_index = i;
            break;
        }
    }

    if(multirun_index >= 0)
        ret = multirun_values[multirun_index][current_run_index];
    else
        ret = GetPot::operator()(VarName, Default);

    return ret;
}

inline int MultirunOptions::GetNumRuns() const
{
    return num_runs;
}

inline void MultirunOptions::SetRun(int run_index)
{
    if(0 <= run_index && run_index < num_runs)
        current_run_index = run_index;
    else
    {   *errstream << "MultirunOptions::SetRun(): run_index out of bounds." << std::endl;
        exit(1);
    }
}

inline std::pair<std::string, int> MultirunOptions::FindZeroParameterRun() const
{
    for(int i = 0; i < multirun_keys.size(); i++)
        for(int run_index = 0; run_index < num_runs; run_index++)
            if(multirun_values[i][run_index] == 0.0)
                return std::make_pair(multirun_keys[i], run_index);

    return std::make_pair("", -1);
}

inline void MultirunOptions::PrintCurrentRunCondition(std::ostream& outstream, const std::string& separator) const
{
    for(int i = 0; i < multirun_keys.size(); i++)
    {
        outstream << multirun_keys[i] << " = " << multirun_values[i][current_run_index] << separator;
    }
}

inline void MultirunOptions::ParseMultirun()
{
    // Get number of multirun variables
    int length = GetPot::vector_variable_size("Multirun");

    // Parse each multirun variable in turn
    for(int i = 0; i < length; i++)
    {
        // Get key and number of runs and check for consistency
        std::string key = GetPot::operator()("Multirun", "", i);
        int num_vals = GetPot::vector_variable_size(key.c_str());

        if(num_vals == 0)
        {   *errstream << "MultirunOptions: " << key << " not found (ignoring)." << std::endl;
        }
        else if(num_vals == 1)
        {   *errstream << "MultirunOptions: " << key << " is a variable of length one (ignoring)." << std::endl;
        }
        else if(multirun_keys.size() && (num_runs != num_vals))
        {   *errstream << "MultirunOptions: " << key << " has wrong length." << std::endl;
            exit(1);
        }
        else
        {   multirun_keys.push_back(key);
            num_runs = num_vals;

            // Get individual values for each run
            std::vector<double> values;
            for(int v = 0; v < num_vals; v++)
                values.push_back(GetPot::operator()(key.c_str(), 0.0, v));
            multirun_values.push_back(values);
        }
    }
}

#endif
