#include "ModelData.hpp"
#include <cstring>
#include "ns3/log.h"
#include "ns3/simulator.h"
#include "ns3/callback.h"

ModelData::ModelData() : parameters(300, 0.0f){}

void serializeModelData(const ModelData& modelData, std::vector<uint8_t>& buffer){
    buffer.resize(modelData.parameters.size() * sizeof(float));
    std::memcpy(buffer.data(), modelData.parameters.data(), buffer.size());
}

bool deserializeModelData(const std::vector<uint8_t>& buffer, ModelData& modelData){
    if (buffer.size() != modelData.parameters.size() * sizeof(float)){
        std::cout << "Error when deserializing data!" << std::endl;
        return false;
    }

    std::memcpy(modelData.parameters.data(), buffer.data(), buffer.size());
    return true;
}

