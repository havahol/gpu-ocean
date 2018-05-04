import unittest
import sys
import time

import xmlrunner
# install xmlrunner by
# $ sudo easy_install unittest-xml-reporting

#import testUtils
from driftTests.CPUDrifter_test import CPUDrifterTest
from driftTests.GPUDrifter_test import GPUDrifterTest

def printSupportedTests():
    print ("Supported tests:")
    print ("0: All, 1: CPUDrifter, 2: GPUDrifter" )


if (len(sys.argv) < 2):
    print ("Usage:")
    print ("\t %s tests  [jenkins]" % sys.argv[0])
    printSupportedTests()
    exit()
tests = int(sys.argv[1])

# In order to format the test report so that Jenkins can read it:
jenkins = False
if (len(sys.argv) > 2):
    if (sys.argv[1].lower() == "jenkins"):
        jenkins = True

if (jenkins):
    unittest.main(testRunner=xmlrunner.XMLTestRunner(output='test-reports'))


# Define the tests that will be part of our test suite:
test_classes_to_run = None
if tests == 0:
    test_classes_to_run = [CPUDrifterTest, GPUDrifterTest]
elif tests == 1:
    test_classes_to_run = [CPUDrifterTest]
elif tests == 2:
    test_classes_to_run = [GPUDrifterTest]
else:
    print ("Error: " + str(tests) + " is not a supported test number...")
    printSupportedTests()
    exit()

    
loader = unittest.TestLoader()
suite_list = []
for test_class in test_classes_to_run:
    suite = loader.loadTestsFromTestCase(test_class)
    suite_list.append(suite)

big_suite = unittest.TestSuite(suite_list)
results = unittest.TextTestRunner(verbosity=2).run(big_suite)

sys.exit(not results.wasSuccessful())
