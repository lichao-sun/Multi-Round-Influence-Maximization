//
//  mi_util.cpp
//  max_influence
//
//  Created by Tian on 1/18/16.
//  Copyright © 2016 Tian. All rights reserved.
//

#include "mi_util.h"



std::istream& PyHelper::GetStdin()
{
    return std::cin;
}

std::ostream& PyHelper::GetStdout()
{
    return std::cout;
}

PyInputStream PyHelper::GetInputFileStream(const char* filename)
{
    std::ifstream* fin = new std::ifstream(filename);
    PyInputStream o;
    o.ptr.reset(fin);
    return o;
}

PyOutputStream PyHelper::GetOutputFileStream(const char* filename)
{
    std::ofstream* fout = new std::ofstream(filename);
    PyOutputStream o;
    o.ptr.reset(fout);
    return o;
}
