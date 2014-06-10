
#include "Testing.h"
#include "BlasHeader.h"
#include <iostream>
#include "CVode/sundials_nvector.h"
#include "ModelRunning.h"
#include <cmath>
#define BOOST_TEST_MODULE MyTest
#include <boost/test/included/unit_test.hpp>

using namespace boost::unit_test;
using namespace std;

BOOST_AUTO_TEST_CASE( testDiffModel ) {
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

    BOOST_CHECK_EQUAL( matlabDiffTPS_pYavg(data, 359.46, GasIn, gridSize, 0.001, params, tps, 1, dIn, 1, 1, 0) , 0);
}




BOOST_AUTO_TEST_CASE( testDiffModelA ) {
    // Just check that diffusion model does not error out.
    
    const size_t gridSize = 50;
    
    double tps[] = { 30 };
    double data[1];
    double GasIn[gridSize];
    
    double params[] = {0.00062517, 5.9978, 0.0041106, 3214.4, 1.0002e-05, 0.11529, 7.2654, 8.3606, 0.0095521, 0.0010002, 0.099954, 0.0016248, 0.094319, 0.021543, 0.9997, 0.034269, 548.02 };
    
    for (size_t ii = 0; ii < gridSize; ii++) GasIn[ii] = ((double) rand()) / ((double) RAND_MAX);
    
    double dIn[] = { 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    BOOST_CHECK_EQUAL( matlabDiffTPS_pYavg(data, 548.02, GasIn, gridSize, 0.034269, params, tps, 1, dIn, 1, 1, 0) , 0);
}


BOOST_AUTO_TEST_CASE( matlabEntry_Test ) {
	double data[1];

	double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
			0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
			0.0085047, 1, 0.0019396, 0.058122, 155.7, 359.46};

	BOOST_TEST_CHECKPOINT( "Calling matlabEntry." );

    BOOST_CHECK_EQUAL( matlabEntry(data, params, 1) , 0);

    BOOST_CHECK_EQUAL( fabs(data[0] - 102.838) < 0.01 , 1);
}


BOOST_AUTO_TEST_CASE( matlabEntryA549_Test ) {
	double data[1];
    
	double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
        0.0085047, 1, 0.058122, 359.46, 0, 0};
    
	BOOST_TEST_CHECKPOINT( "Calling matlabEntryA549." );
    
    BOOST_CHECK_EQUAL( matlabEntryA549(data, params, 1) , 0);
    
    BOOST_TEST_MESSAGE( "A549 result: " << data[0] );
}





BOOST_AUTO_TEST_CASE( matlabEntry_Test_SI ) {
	double data[1];
    
	double params[] = { 1.2, 0.054435, 0.042, 24.392, 0.00081113, 0.34571,
        0.0010493, 0.017322, 1e-06, 3.183, 0.0056061, 0.002045, 0.1,
        0.0085047, 1, 0.0019396, 0.058122, 155.7, 359.46};
    
	BOOST_TEST_CHECKPOINT( "Calling matlabEntry." );
    
    BOOST_CHECK_EQUAL( matlabEntryWithSi(data, params, 1) , 0);
    
    BOOST_TEST_MESSAGE( "Si output: " << data[0] );
}



