//
// Created by imartin on 03-Apr-22.
//

#ifndef TVB_CPP_CSV_TOOLS_H
#define TVB_CPP_CSV_TOOLS_H

#include "definitions.h"

namespace tvb {
    tvb::TArray2d csv_load(const std::string &path);

    void csv_save(const std::string &path, const tvb::TArray2d& data);

    template<typename Numeric>
    void csv_save(const std::string &filename, const std::vector<std::vector<Numeric>>& data) {
        std::ofstream myFile(filename);

        // Send data to the stream
        for(int i = 0; i < data.size(); ++i) {
            for (int j = 0; j < data[i].size(); ++j) {
                if (j > 0)
                    myFile << ", ";
                myFile << data[i][j];
            }
            myFile << std::endl;
        }
        // Close the file
        myFile.close();

    }
}

#endif //TVB_CPP_CSV_TOOLS_H
