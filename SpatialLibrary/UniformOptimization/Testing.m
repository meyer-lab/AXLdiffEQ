//
//  Testing.m
//  UniformOptimization
//
//  Created by Aaron Meyer on 6/10/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#import <XCTest/XCTest.h>
#import "BlasHeader.h"


@interface Testing : XCTestCase

@end

@implementation Testing

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

//- (void)testDiffModel {
//    // This is an example of a functional test case.
//    
//    // Just check that diffusion model does not error out.
//    
//    const size_t gridSize = 50;
//    
//    double tps[] = { 30 };
//    double *ttps = tps;
//    double data[1];
//    double *ddata = data;
//    double GasIn[gridSize];
//    
//    double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
//        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
//        0.0085047, 1, 0.0019396, 0.058122, 155.7, 359.46 };
//    double *pparams = params;
//    
//    for (size_t ii = 0; ii < gridSize; ii++) GasIn[ii] = ((double) rand()) / ((double) RAND_MAX);
//    double *GGasIn = GasIn;
//    
//    double dIn[] = { 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//    double *ddIn = dIn;
//    
////    XCTAssert( matlabDiffTPS_pYavg(data, 359.46, GasIn, gridSize, 0.001, params, tps, 1, dIn, 1, 1, 0) == 0, @"Pass");
////    
////    [self measureBlock:^{
////        matlabDiffTPS_pYavg(ddata, 359.46, GGasIn, gridSize, 0.001, pparams, ttps, 1, ddIn, 1, 1, 0);
////    }];
//}


- (void)testMatlabEntry {
    __block double dP = 0;
    
    double params[] = {1.0, 1.0, 1.0, 1.0, 1.03695332e-01,16,7.59242750e-04,0.26,1.56016470e-04,   1.14961764e-02,   6.67336796e-02,   9.32114079e+01, 6.05123140e-02, 1};
    double *pP = &params[0];
    
    [self measureBlock:^{
        
        for (int ii = 0; ii < 10; ii++)
            dP = pyEntry(pP);
    }];
    
    XCTAssertEqualWithAccuracy(dP, 5751.623, 0.01, @"Pass");
}


//- (void)testProfile {
//    double data[1];
//    
//    double params[13] = {1.24e-01, 1.44e+03, 8.52e-03, 4.4147e-01, 6.6524e+00, 3.1110e-04, 1.6978e-01, 4.1833e-01, 1.4011e-02, 6.6132e-03, 4.6966e-02, 3.5584e+02, 1.5983e-02};
//    
//    double tps[] = {10};
//    
//    XCTAssert( calcProfileMatlab(data, params, tps, 1, 10, 1) == 0, @"Pass");
//    
//    //XCTAssertEqualWithAccuracy(data[0], 102.838, 0.01, @"Pass");
//}

//- (void) testDiffModelA {
//    // Just check that diffusion model does not error out.
//
//    double params[13] = {1.24e-01, 1.44e+03, 8.52e-03, 4.4147e-01, 6.6524e+00, 3.1110e-04, 1.6978e-01, 4.1833e-01, 1.4011e-02, 6.6132e-03, 4.6966e-02, 3.5584e+02, 1.5983e-02};
//    double *pPtr = params;
//
//    double dIn[] = { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//    double *dInPtr = dIn;
//    __block double retVal = 1;
//
//    [self measureBlock:^{
//        retVal = pyDiffTPS_Activation(16, pPtr, 30, dInPtr);
//    }];
//    
//    XCTAssertEqualWithAccuracy(retVal, 102.838, 1000, @"Pass");
//}


//- (void) testDiffModelA {
//    // Just check that diffusion model does not error out.
//    
//    const size_t gridSize = 50;
//    
//    double tps[] = { 30 };
//    double *tPtr = tps;
//    double data[1];
//    double GasIn[gridSize];
//    double *GasInPtr = GasIn;
//    double *dataPtr = data;
//    
//    double params[] = {0.00062517, 5.9978, 0.0041106, 3214.4, 1.0002e-05, 0.11529, 7.2654, 8.3606, 0.0095521, 0.0010002, 0.099954, 0.0016248, 0.094319, 0.021543, 0.9997, 0.034269, 548.02 };
//    
//    double *pPtr = params;
//    
//    for (size_t ii = 0; ii < gridSize; ii++) GasIn[ii] = ((double) rand()) / ((double) RAND_MAX);
//    
//    double dIn[] = { 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//    double *dInPtr = dIn;
//    __block int retVal = 1;
//    
////    [self measureBlock:^{
////        retVal = matlabDiffTPS_pYavg(dataPtr, 548.02, GasInPtr, gridSize, 0.034269, pPtr, tPtr, 1, dInPtr, 1, 1, 0);
////    }];
////    
////    XCTAssert( retVal == 0, @"Pass");
//}


@end
