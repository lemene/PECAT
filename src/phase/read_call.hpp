#pragma once

#include <string>
#include <unordered_map>

#include "sequence.hpp"
#include "utils/program.hpp"
#include "utils/project_file.hpp"
#include "utils/string_pool.hpp"


namespace fsa {


class ReadStore;

  
class ReadCall : public Program {
public:
    virtual ~ReadCall() {}
    bool ParseArgument(int argc, const char* const argv[]);
    virtual void Running();
    void Usage();

protected:
    ArgumentParser GetArgumentParser();
    void CheckArguments();

protected:
    std::string vcf_fname_;

    StringPool string_pool_;

    SnpStore snp_store_ { string_pool_ };
};

} // namespace fsa {     
