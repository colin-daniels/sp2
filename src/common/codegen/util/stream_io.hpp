#ifndef CODEGEN_OUTPUT_HPP
#define CODEGEN_OUTPUT_HPP

#include <fstream>
#include <vector>
#include <string>

void write_array(std::ostream &out, const std::vector<double> &data,
    const std::string prefix = "\t", const std::string suffix = "\n");

#endif // CODEGEN_OUTPUT_HPP
