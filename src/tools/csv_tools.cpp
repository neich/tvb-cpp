//
// Created by imartin on 03-Apr-22.
//

#include <tools/csv_tools.h>

using namespace tvb;
using namespace std;

TArray2d tvb::csv_load(const string& path) {
    vector<Float> matrixEntries;
    //vector<double> matrixEntries;

    // in this object we store the data from the matrix
    ifstream matrixDataFile(path);

    if (!matrixDataFile.good())
        throw std::runtime_error(string_format("Cannot open file: %s", path.c_str()));

    // this variable is used to store the row of the matrix that contains commas
    string matrixRowString;

    // this variable is used to store the matrix entry;
    string matrixEntry;

    // this variable is used to track the number of rows
    int matrixRowNumber = 0;


    while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
    {
        stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

        while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    }

    // here we convet the vector variable into the matrix and return the resulting object,
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
    return Eigen::Map<TArray2d>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
}

void tvb::csv_save(const std::string &filename, const tvb::TArray2d& data) {
    std::ofstream myFile(filename);

    if (!myFile.good())
        throw std::runtime_error(string_format("Cannot open file: %s", filename.c_str()));

    // Send data to the stream
    for(int i = 0; i < data.rows(); ++i) {
        for (int j = 0; j < data.row(i).size(); ++j) {
            if (j > 0)
                myFile << ", ";
            myFile << data(i, j);
        }
        myFile << std::endl;
    }
    // Close the file
    myFile.close();
}