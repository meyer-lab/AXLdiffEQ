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

- (void)testDiffModel {
    // This is an example of a functional test case.
    
    // Just check that diffusion model does not error out.
    
    const size_t gridSize = 50;
    
    double tps[] = { 30 };
    double data[1];
    double GasIn[gridSize];
    
    double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
        0.0085047, 1, 0.0019396, 0.058122, 155.7, 359.46 };
    
    for (size_t ii = 0; ii < gridSize; ii++) GasIn[ii] = ((double) rand()) / ((double) RAND_MAX);
    
    double dIn[] = { 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    XCTAssert( matlabDiffTPS_pYavg(data, 359.46, GasIn, gridSize, 0.001, params, tps, 1, dIn, 1, 1, 0) == 0, @"Pass");
}

- (void)testPerformancematlabEntryA549 {
    double data[1];
    
    double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
        0.0085047, 1, 0.058122, 359.46, 0, 0};
    
    void *paramP = params;
    void *dataP = data;
    
    __block int retVal = 1;
    
    
    [self measureBlock:^{
        retVal = matlabEntryA549(dataP, paramP, 1);
    }];
    
    XCTAssert( retVal == 0, @"Pass");
}

- (void)testPerformancematlabEntryBT549py {
    
    double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
        0.0085047, 1, 0.058122, 359.46, 0.5, 600};
    
    void *paramP = params;
    
    __block double retVal = 1;
    
    
    [self measureBlock:^{
        retVal = matlabEntryBT549VaryEndoPy(paramP);
    }];
    
    //XCTAssert( retVal == 0, @"Pass");
}

- (void)testPerformancematlabEntryA549pyRed {
    double params[] = { 0.054, 24, 0.00081, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
        0.0085047, 0.058122, 359.46, 0.5, 600};
    
    double retVal = matlabEntryA549VaryEndoPyRed(params);

    //XCTAssertEqualWithAccuracy(retVal, 344.04, 0.01, @"Pass");
}


- (void)testPerformanceMultiPy {
    
    double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045,
        0.0085047, 1, 0.5, 600, 359.46, 359.46, 359.46, 359.46, 0.1, 0.1, 0.1, 0.1, 0.058122, 0.058122, 0.058122, 0.058122};
    
    void *paramP = params;
    
    __block double retVal = 1;
    
    
    [self measureBlock:^{
        retVal = multiPyEntry(paramP);
    }];
    
    //XCTAssertEqual( retVal, 0, @"Pass");
}


- (void)testMatlabEntry {
    double data[1];
    
    double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
        0.0085047, 1, 0.0019396, 0.058122, 155.7, 359.46};
    
    XCTAssert( matlabEntry(data, params, 1) == 0, @"Pass");
    
    XCTAssertEqualWithAccuracy(data[0], 102.838, 0.01, @"Pass");
}



- (void) testDiffModelA {
    // Just check that diffusion model does not error out.
    
    const size_t gridSize = 50;
    
    double tps[] = { 30 };
    double *tPtr = tps;
    double data[1];
    double GasIn[gridSize];
    double *GasInPtr = GasIn;
    double *dataPtr = data;
    
    double params[] = {0.00062517, 5.9978, 0.0041106, 3214.4, 1.0002e-05, 0.11529, 7.2654, 8.3606, 0.0095521, 0.0010002, 0.099954, 0.0016248, 0.094319, 0.021543, 0.9997, 0.034269, 548.02 };
    
    double *pPtr = params;
    
    for (size_t ii = 0; ii < gridSize; ii++) GasIn[ii] = ((double) rand()) / ((double) RAND_MAX);
    
    double dIn[] = { 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double *dInPtr = dIn;
    __block int retVal = 1;
    
    //[self measureBlock:^{
        retVal = matlabDiffTPS_pYavg(dataPtr, 548.02, GasInPtr, gridSize, 0.034269, pPtr, tPtr, 1, dInPtr, 1, 1, 0);
    //}];
    
    XCTAssert( retVal == 0, @"Pass");
}


/*- (void)testPerformanceExample {
 // This is an example of a performance test case.
 [self measureBlock:^{
 // Put the code you want to measure the time of here.
 }];
 }*/

@end
