#include "RestartObject.h"

#include <fstream>

#include "Parameter/Parameter.h"
#include "basics/utilities/UbMath.h"

void RestartObject::deserialize(const std::string &filename, std::shared_ptr<Parameter>& para)
{
    deserialize_internal(filename);

    for (int j = para->getCoarse(); j <= para->getFine(); j++) {
        std::vector<real> vec;
        fs.push_back(vec);

        for (unsigned int i = 0; i < (para->getD3Qxx() * para->getParH(j)->size_Mat_SP); i++) {
            para->getParH(j)->d0SP.f[0][i] = fs[j][i];
        }
    }
}

void RestartObject::serialize(const std::string &filename, const std::shared_ptr<Parameter>& para)
{
    if (fs.size() > 0) {
        clear(para);
    }
    for (int j = para->getCoarse(); j <= para->getFine(); j++) {
        std::vector<real> vec;
        fs.push_back(vec);

        for (unsigned int i = 0; i < (para->getD3Qxx() * para->getParH(j)->size_Mat_SP); i++) {
            if (UbMath::isNaN(para->getParH(j)->d0SP.f[0][i])) {
                fs[j].push_back((real)0.0);
            } else {
                fs[j].push_back(para->getParH(j)->d0SP.f[0][i]);
            }
        }
    }

    serialize_internal(filename);
}

void RestartObject::clear(const std::shared_ptr<Parameter>& para)
{
    for (int j = para->getCoarse(); j <= para->getFine(); j++) {
        fs[j].resize(0);
    }
    fs.resize(0);
}
//////////////////////////////////////////////////////////////////////////

void ASCIIRestartObject::serialize_internal(const std::string &filename)
{
    std::ofstream stream(filename + ".txt");

    if (!stream) {
        stream.clear();
        std::string path = UbSystem::getPathFromString(filename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            stream.open(filename.c_str());
        }

        if (!stream)
            throw UbException(UB_EXARGS, "couldn't open file " + filename);
    }

    for (std::vector<std::vector<real>>::size_type i = 0; i < fs.size(); i++) {
        for (std::vector<real>::size_type j = 0; j < fs[i].size(); j++) {
            stream << fs[i][j] << " ";
        }
        stream << "\n";
    }

    stream.close();
}

void ASCIIRestartObject::deserialize_internal(const std::string &filename)
{
    std::ifstream stream(filename + ".txt", std::ios_base::in);

    if (!stream.is_open())
        throw UbException(UB_EXARGS, "couldn't open check point file " + filename);

    std::string line;
    while (std::getline(stream, line)) {
        std::vector<real> lineData;
        std::stringstream lineStream(line);

        real value;
        while (lineStream >> value) {
            lineData.push_back(value);
        }
        fs.push_back(lineData);
    }

    stream.close();
}

void BinaryRestartObject::serialize_internal(const std::string &filename)
{
    std::ofstream stream(filename + ".bin", std::ios_base::binary);

    if (!stream) {
        stream.clear();
        std::string path = UbSystem::getPathFromString(filename);
        if (path.size() > 0) {
            UbSystem::makeDirectory(path);
            stream.open(filename.c_str());
        }

        if (!stream)
            throw UbException(UB_EXARGS, "couldn't open file " + filename);
    }

    // Store size of the outer vector
    int s1 = fs.size();
    stream.write(reinterpret_cast<const char *>(&s1), sizeof(s1));

    // Now write each vector one by one
    for (auto &v : fs) {
        // Store its size
        int size = v.size();
        stream.write(reinterpret_cast<const char *>(&size), sizeof(size));

        // Store its contents
        stream.write(reinterpret_cast<const char *>(&v[0]), v.size() * sizeof(float));
    }

    stream.close();
}

void BinaryRestartObject::deserialize_internal(const std::string &filename)
{
    std::ifstream stream(filename + ".bin", std::ios_base::in | std::ifstream::binary);

    if (!stream.is_open())
        throw UbException(UB_EXARGS, "couldn't open check point file " + filename);

    int size = 0;
    stream.read(reinterpret_cast<char *>(&size), sizeof(size));
    fs.resize(size);
    for (int n = 0; n < size; ++n) {
        int size2 = 0;
        stream.read(reinterpret_cast<char *>(&size2), sizeof(size2));
        real f;
        for (int k = 0; k < size2; ++k) {
            stream.read(reinterpret_cast<char *>(&f), sizeof(f));
            fs[n].push_back(f);
        }
    }

    stream.close();
}
