

#include "../external/catch.hpp"
#include "../filesystem.h"

TEST_CASE("Filesystem join, parent, stem, extension", "[filesystem]") {
    REQUIRE(filesystem::stem("xyz.txt") == "xyz");
    REQUIRE(filesystem::stem("xyz") == "xyz");
    REQUIRE(filesystem::stem("/xyz") == "xyz");
    REQUIRE(filesystem::stem("/xyz.abc") == "xyz");
    REQUIRE(filesystem::stem("/xyz.abc.xyz") == "xyz.abc");
    REQUIRE(filesystem::stem("/dsdsd/xyz.abc.xyz") == "xyz.abc");

    REQUIRE(filesystem::parent("xyz.txt") == "");
    REQUIRE(filesystem::parent("/xyz.txt") == "/");
    REQUIRE(filesystem::parent("dir/xyz.txt") == "dir");
    REQUIRE(filesystem::parent("/dir/xyz.txt") == "/dir");

    // TODO further parent tests on directories (which must exist)...

    REQUIRE(filesystem::extension("xyz.txt") == "txt");
    REQUIRE(filesystem::extension("xyz") == "");

    REQUIRE(filesystem::filename("/xyz.txt") == "xyz.txt");
    REQUIRE(filesystem::filename("ddd/dsds/xyz.txt") == "xyz.txt");
    REQUIRE(filesystem::filename("xyz.txt") == "xyz.txt");
}