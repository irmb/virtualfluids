#ifndef QDEBUG_HPP
#define QDEBUG_HPP

#include <stdio.h>
#include <fstream>
#include <sstream>
// #include <math.h>
#include <cmath>
#include "Calculation/Calculation.h"
#include "lbm/constants/D3Q27.h"
#include "Parameter/Parameter.h"
#include "basics/utilities/UbSystem.h"
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include "StringUtilities/StringUtil.h"

namespace QDebugWriter
{
    void writeQValues(QforBoundaryConditions &Q, int* k, int kq, const std::string &name)
    {

        std::vector< std::vector<real> > qs;
        for (int j = 0; j < kq; j++)
        {
            uint32_t qKey = 0;
            std::vector<real> qNode;

            for (int i = 26; i >= 0; i--)
            {
                real q = Q.q27[i][j];
                if (q > 0) {
                    qKey += (uint32_t)pow(2, 26 - i);
                    qNode.push_back(q);
                }
            }
            if (qKey > 0) {
                float transportKey = *((float*)&qKey);
                qNode.push_back((real)transportKey);
                qNode.push_back((real)k[j]);
                qs.push_back(qNode);
            }
            qNode.clear();

        }

        SPtr<std::ofstream> outQ(new std::ofstream);
        outQ->open(name.c_str(), std::ios::out | std::ios::binary);



        for (std::size_t index = 0; index < qs.size(); index++) {
            std::vector<real> bcs = qs[index];
            uint32_t key = *((uint32_t*)&bcs[bcs.size() - 2]);
            int qIndex = (int)bcs[bcs.size() - 1];

            *outQ << qIndex << " " << key;

            for (std::size_t i = 0; i < bcs.size() - 2; i++) {
                *outQ << " " << std::fixed << std::setprecision(16) << bcs[i];
            }

            *outQ << "\n";
        }

        outQ->close();

    }
}

#endif
