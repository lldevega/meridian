/// Common types and functions for the Meridian solver
/*
 *
 */

#include <string>
#include <vector>
#include <map>
#include <cassert>

namespace meridian
{
/// Typename of a vector of strings.
using StringVector = typename std::vector<std::string>;
/// Typename of a vector of doubles.
using DoubleVector = typename std::vector<double>;
/// Typename of a vector of integers.
using IntVector = typename std::vector<int>;

/// Typename of a map holding keys (strings) and values (doubles).
using ValueDict = typename std::map<std::string, double>;
/// Typename of a map holding keys (strings) and values (strings).
using StringDict = typename std::map<std::string, std::string>;

/// Check whether a given value dict contains a expected list of keys.
/*
 * Requires the value dict to be examined and the list of expected keys.
 * Stops the code execution if any of the expected keys is not found.
 * The expected keys get copied as they are usually passed from an on-the-fly
 * instantiated object.
 */
void CheckValueDict(const ValueDict &valueDict, const StringVector expectedKeys,
	const char* msg = "Not all the expected keys are present in the value dict")
{
	int foundKeyNumber = 0;
	int expectedKeyNumber = expectedKeys.size();

	// loop over the expected keys and check if they are present of the value dict
	for (auto it = expectedKeys.begin(); it != expectedKeys.end(); ++it)
	    foundKeyNumber += valueDict.count(*it);

	assert(foundKeyNumber == expectedKeyNumber && msg);
}

/// Check whether a given string dict contains a expected list of keys.
/*
 * Requires the value dict to be examined and the list of expected keys.
 * Stops the code execution if any of the expected keys is not found.
 * The expected keys get copied as they are usually passed from an on-the-fly
 * instantiated object.
 */
void CheckStringDict(const StringDict &stringDict, const StringVector expectedKeys,
	const char* msg = "Not all the expected keys are present in the value dict")
{
	int foundKeyNumber = 0;
	int expectedKeyNumber = expectedKeys.size();

	// loop over the expected keys and check if they are present of the value dict
	for (auto it = expectedKeys.begin(); it != expectedKeys.end(); ++it)
	    foundKeyNumber += stringDict.count(*it);

	assert(foundKeyNumber == expectedKeyNumber && msg);
}

} // end of the meridian namespace
