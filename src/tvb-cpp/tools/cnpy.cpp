//Copyright (C) 2011  Carl Rogers
//Released under MIT License
//license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#include "cnpy.h"
#include <complex>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdint.h>
#include <stdexcept>
#include <regex>

cnpy::NpyArray load_the_npz_array(std::fstream& fp, uint32_t compr_bytes, uint32_t uncompr_bytes);


char cnpy::BigEndianTest() {
    int x = 1;
    return (((char *)&x)[0]) ? '<' : '>';
}

char cnpy::map_type(const std::type_info& t)
{
    if(t == typeid(float) ) return 'f';
    if(t == typeid(double) ) return 'f';
    if(t == typeid(long double) ) return 'f';

    if(t == typeid(int) ) return 'i';
    if(t == typeid(char) ) return 'i';
    if(t == typeid(short) ) return 'i';
    if(t == typeid(long) ) return 'i';
    if(t == typeid(long long) ) return 'i';

    if(t == typeid(unsigned char) ) return 'u';
    if(t == typeid(unsigned short) ) return 'u';
    if(t == typeid(unsigned long) ) return 'u';
    if(t == typeid(unsigned long long) ) return 'u';
    if(t == typeid(unsigned int) ) return 'u';

    if(t == typeid(bool) ) return 'b';

    if(t == typeid(std::complex<float>) ) return 'c';
    if(t == typeid(std::complex<double>) ) return 'c';
    if(t == typeid(std::complex<long double>) ) return 'c';

    else return '?';
}

template<> std::vector<char>& cnpy::operator+=(std::vector<char>& lhs, const std::string rhs) {
    lhs.insert(lhs.end(),rhs.begin(),rhs.end());
    return lhs;
}

template<> std::vector<char>& cnpy::operator+=(std::vector<char>& lhs, const char* rhs) {
    //write in little endian
    size_t len = strlen(rhs);
    lhs.reserve(len);
    for(size_t byte = 0; byte < len; byte++) {
        lhs.push_back(rhs[byte]);
    }
    return lhs;
}

void cnpy::parse_npy_header(unsigned char* buffer,size_t& word_size, std::vector<size_t>& shape, bool& fortran_order) {
    //std::string magic_string(buffer,6);
    uint8_t major_version = *reinterpret_cast<uint8_t*>(buffer+6);
    uint8_t minor_version = *reinterpret_cast<uint8_t*>(buffer+7);
    uint16_t header_len = *reinterpret_cast<uint16_t*>(buffer+8);
    std::string header(reinterpret_cast<char*>(buffer+9),header_len);

    size_t loc1, loc2;

    //fortran order
    loc1 = header.find("fortran_order")+16;
    fortran_order = (header.substr(loc1,4) == "True" ? true : false);

    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1+1,loc2-loc1-1);
    while(std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stoi(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    //endian, word size, data type
    //byte order code | stands for not applicable.
    //not sure when this applies except for byte array
    loc1 = header.find("descr")+9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    assert(littleEndian);

    //char type = header[loc1+1];
    //assert(type == map_type(T));

    std::string str_ws = header.substr(loc1+2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0,loc2).c_str());
}

void cnpy::parse_npy_header(std::fstream& fp, size_t& word_size, std::vector<size_t>& shape, bool& fortran_order) {
    char buffer[256];

    fp.read(buffer, 11);
    if (!fp.good()) throw std::runtime_error("parse_npy_header: failed fread");

    fp.getline(buffer, 256);

    if (fp.fail())
        throw std::runtime_error("parse_npy_header: failed to read ny header");

    std::string header(buffer);

    size_t loc1, loc2;

    //fortran order
    loc1 = header.find("fortran_order");
    if (loc1 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'fortran_order'");
    loc1 += 16;
    fortran_order = (header.substr(loc1,4) == "True" ? true : false);

    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");
    if (loc1 == std::string::npos || loc2 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: '(' or ')'");

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1+1,loc2-loc1-1);
    while(std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stoi(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    //endian, word size, data type
    //byte order code | stands for not applicable.
    //not sure when this applies except for byte array
    loc1 = header.find("descr");
    if (loc1 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'descr'");
    loc1 += 9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    assert(littleEndian);

    //char type = header[loc1+1];
    //assert(type == map_type(T));

    std::string str_ws = header.substr(loc1+2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0,loc2).c_str());
}

void cnpy::parse_zip_footer(std::fstream& fp, uint16_t& nrecs, size_t& global_header_size, size_t& global_header_offset)
{
    std::vector<char> footer(22);
    fp.seekg(-22, std::ios::end);
    fp.read(&footer[0], 22);
    if (!fp)
        throw std::runtime_error("parse_zip_footer: failed fread");

    uint16_t disk_no, disk_start, nrecs_on_disk, comment_len;
    disk_no = *(uint16_t*) &footer[4];
    disk_start = *(uint16_t*) &footer[6];
    nrecs_on_disk = *(uint16_t*) &footer[8];
    nrecs = *(uint16_t*) &footer[10];
    global_header_size = *(uint32_t*) &footer[12];
    global_header_offset = *(uint32_t*) &footer[16];
    comment_len = *(uint16_t*) &footer[20];

    assert(disk_no == 0);
    assert(disk_start == 0);
    assert(nrecs_on_disk == nrecs);
    assert(comment_len == 0);
}

cnpy::NpyArray load_the_npy_file(std::fstream& fp) {
    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    cnpy::parse_npy_header(fp,word_size,shape,fortran_order);

    cnpy::NpyArray arr(shape, word_size, fortran_order);
    fp.read(arr.data<char>(), arr.num_bytes());
    if (!fp)
        throw std::runtime_error("load_the_npy_file: failed fread");

    return arr;
}

cnpy::NpyArray load_the_npz_array(std::fstream& fp, uint32_t compr_bytes, uint32_t uncompr_bytes) {

    std::vector<char> buffer_compr(compr_bytes);
    std::vector<char> buffer_uncompr(uncompr_bytes);
    try {
        fp.read(&buffer_compr[0], compr_bytes);
    }
    catch (std::runtime_error& e) {
        throw std::runtime_error("load_the_npy_file: failed fread");
    }

    int err;
    z_stream d_stream;

    d_stream.zalloc = Z_NULL;
    d_stream.zfree = Z_NULL;
    d_stream.opaque = Z_NULL;
    d_stream.avail_in = 0;
    d_stream.next_in = Z_NULL;
    err = inflateInit2(&d_stream, -MAX_WBITS);

    d_stream.avail_in = compr_bytes;
    d_stream.next_in = (unsigned char*)&buffer_compr[0];
    d_stream.avail_out = uncompr_bytes;
    d_stream.next_out = (unsigned char*)&buffer_uncompr[0];

    err = inflate(&d_stream, Z_FINISH);
    err = inflateEnd(&d_stream);

    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    cnpy::parse_npy_header((unsigned char*)&buffer_uncompr[0], word_size, shape, fortran_order);

    cnpy::NpyArray array(shape, word_size, fortran_order);

    size_t offset = uncompr_bytes - array.num_bytes();
    memcpy(array.data<unsigned char>(), &buffer_uncompr[0] + offset, array.num_bytes());

    return array;
}


cnpy::npz_t cnpy::npz_load(const std::string& fname) {
    std::fstream fp(fname, std::ios::binary | std::ios::in);

    if(fp.bad()) {
        throw std::runtime_error("npz_load: Error! Unable to open file "+fname+"!");
    }

    cnpy::npz_t arrays;

    while(1) {
        std::vector<char> local_header(30);
        fp.read(&local_header[0], 30);
        if (fp.fail())
            throw std::runtime_error("npz_load: error reading  header");

        //if we've reached the global header, stop reading
        if(local_header[2] != 0x03 || local_header[3] != 0x04) break;

        //read in the variable name
        uint16_t name_len = *(uint16_t*) &local_header[26];
        std::string varname(name_len,' ');
        fp.read(&varname[0], name_len);
        if (!fp)
            throw std::runtime_error("npz_load: error reading variable name");

        //erase the lagging .npy
        varname.erase(varname.end()-4,varname.end());

        //read in the extra field
        uint16_t extra_field_len = *(uint16_t*) &local_header[28];
        if(extra_field_len > 0) {
            std::vector<char> buff(extra_field_len);
            fp.read(&buff[0], extra_field_len);
            if (!fp)
                throw std::runtime_error("npz_load: error reading extra field");
        }

        uint16_t compr_method = *reinterpret_cast<uint16_t*>(&local_header[0]+8);
        uint32_t compr_bytes = *reinterpret_cast<uint32_t*>(&local_header[0]+18);
        uint32_t uncompr_bytes = *reinterpret_cast<uint32_t*>(&local_header[0]+22);

        if(compr_method == 0) {arrays[varname] = load_the_npy_file(fp);}
        else {arrays[varname] = load_the_npz_array(fp,compr_bytes,uncompr_bytes);}
    }

    fp.close();
    return arrays;
}

cnpy::NpyArray cnpy::npz_load(const std::string& fname, const std::string& varname) {
    std::fstream fp(fname, std::ios::binary | std::ios::in);

    if(fp.bad())
        throw std::runtime_error("npz_load: Unable to open file "+fname);

    while(true) {
        std::vector<char> local_header(30);
        fp.read(&local_header[0], 30);
        if (!fp)
            throw std::runtime_error("npz_load: error reading  header");

        //if we've reached the global header, stop reading
        if(local_header[2] != 0x03 || local_header[3] != 0x04) break;

        //read in the variable name
        uint16_t name_len = *(uint16_t*) &local_header[26];
        std::string vname(name_len,' ');
        fp.read(&vname[0], name_len);
        if (!fp)
            throw std::runtime_error("npz_load: error reading variable name");
        vname.erase(vname.end()-4,vname.end()); //erase the lagging .npy

        //read in the extra field
        uint16_t extra_field_len = *(uint16_t*) &local_header[28];
        fp.seekg(extra_field_len, std::ios_base::cur); //skip past the extra field

        uint16_t compr_method = *reinterpret_cast<uint16_t*>(&local_header[0]+8);
        uint32_t compr_bytes = *reinterpret_cast<uint32_t*>(&local_header[0]+18);
        uint32_t uncompr_bytes = *reinterpret_cast<uint32_t*>(&local_header[0]+22);

        if(vname == varname) {
            NpyArray array  = (compr_method == 0) ? load_the_npy_file(fp) : load_the_npz_array(fp,compr_bytes,uncompr_bytes);
            fp.close();
            return array;
        }
        else {
            //skip past the data
            uint32_t size = *(uint32_t*) &local_header[22];
            fp.seekg(size, std::ios_base::cur);
        }
    }

    fp.close();

    //if we get here, we haven't found the variable in the file
    throw std::runtime_error("npz_load: Variable name "+varname+" not found in "+fname);
}

cnpy::NpyArray cnpy::npy_load(const std::string& fname) {

    std::fstream fp(fname, std::ios::binary | std::ios::in | std::ios::out);

    if(!fp.good()) 
        throw std::runtime_error("npy_load: Unable to open file "+fname);

    NpyArray arr = load_the_npy_file(fp);

    fp.close();
    return arr;
}



