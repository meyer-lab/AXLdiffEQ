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

- (void)testMatlabEntry {
    __block double dP = 0;
    
    double params[] = {27.003, 5.0203e-06, 4.4569, 0.0018128, 4.467, 0.07214, 0.0010026, 0.0017783, 0.19216, 0.13615, 10190, 0.038456, 1};
    double *pP = &params[0];
    
    [self measureBlock:^{
        for (int ii = 0; ii < 2; ii++)
            dP = pyEntry(pP);
    }];
}



- (void)testMatlabEntryVal {
    double dP = 0;
    
    double params[] = {27.003, 5.0203e-06, 4.4569, 0.0018128, 4.467, 0.07214, 0.0010026, 0.0017783, 0.19216, 0.13615, 10190, 0.038456, 1};
    
    dP = pyEntry(params);

    XCTAssertEqualWithAccuracy(dP, 2882.157, 0.01, @"Pass");
}

@end
