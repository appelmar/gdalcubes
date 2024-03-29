/*
    MIT License

    Copyright (c) 2019 Marius Appel <marius.appel@hs-bochum.de>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#include "../external/catch.hpp"
#include "../filesystem.h"

using namespace gdalcubes;

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
