#include <UnitTest++.h>
#include "TestReporterStdout.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>
#include "json.hpp"

using json = nlohmann::json;


int main( int argc, char** argv )
{
  if( argc > 1 )
  {
      //if first arg is "suite", we search for suite names instead of test names
    const bool bSuite = strcmp( "suite", argv[ 1 ] ) == 0;

      //walk list of all tests, add those with a name that
      //matches one of the arguments  to a new TestList
    const UnitTest::TestList& allTests( UnitTest::Test::GetTestList() );
    UnitTest::TestList selectedTests;
    UnitTest::Test* p = allTests.GetHead();
    while( p )
    {
      for( int i = 1 ; i < argc ; ++i )
        if( strcmp( bSuite ? p->m_details.suiteName
                           : p->m_details.testName, argv[ i ] ) == 0 )
          selectedTests.Add( p );
      p = p->m_nextTest;
    }

      //run selected test(s) only
    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner( reporter );
    return runner.RunTestsIf( selectedTests, 0, UnitTest::True(), 0 );
  }
  else
  {
    return UnitTest::RunAllTests();
  }
}