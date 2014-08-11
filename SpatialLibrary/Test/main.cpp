//
//  main.cpp
//  Test
//
//  Created by Aaron Meyer on 7/1/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include <iostream>
#include "BlasHeader.h"



int main () {
    
    double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
        0.0085047, 1, 0.058122, 359.46, 0, 0};
    
    double data[1];
    double retVal;

    for (size_t tt = 0; tt < 1000; tt++) {
        retVal = matlabEntryA549(data, params, 1);
    }
    
    
    
    return 0;
}