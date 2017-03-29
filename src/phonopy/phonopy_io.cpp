#include "phonopy_io.hpp"
#include "yaml-cpp/yaml.h"

#include <map>
#include <iostream>

void print_type(const YAML::Node &node)
{
    std::string type = "";
    switch (node.Type()) {
    case YAML::NodeType::Null:      type = "Null"; break;
    case YAML::NodeType::Scalar:    type = "Scalar"; break;
    case YAML::NodeType::Sequence:  type = "Sequence"; break;
    case YAML::NodeType::Map:       type = "Map"; break;
    case YAML::NodeType::Undefined: type = "Undefined"; break;
    }

    std::cerr << "Node type: " << type << std::endl;
}

std::vector<std::string> sp2::phonopy::read_irreps(std::string filename)
{
    YAML::Node node = YAML::LoadFile(filename);

    if (!node || !node["normal_modes"])
        return {};


    std::map<int, std::string> irrep_labels;
    for (const auto &elem : node["normal_modes"])
    {
        if (!elem["band_indices"])
            return {};

        std::string label = elem["ir_label"] ?
            elem["ir_label"].as<std::string>() : "Unknown";

        for (const auto &idx_node : elem["band_indices"])
        {
            int idx = idx_node.as<int>();
            irrep_labels[idx] = label;
        }
    }

    std::vector<std::string> output;
    for (auto &elem : irrep_labels)
        output.push_back(elem.second);

    return output;
}
