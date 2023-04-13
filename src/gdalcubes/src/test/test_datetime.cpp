/*
    MIT License

    Copyright (c) 2019 Marius Appel <marius.appel@uni-muenster.de>

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
#include <string>

#include "../datetime.h"
#include "../external/catch.hpp"

using namespace gdalcubes;

TEST_CASE("Deriving datetime unit from string", "[datetime]") {
    REQUIRE(datetime::from_string("2002-03-04 12:13:14").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_string("2002-03-04 12:13").unit() == datetime_unit::MINUTE);
    REQUIRE(datetime::from_string("2002-03-04 12").unit() == datetime_unit::HOUR);
    REQUIRE(datetime::from_string("2002-03-04").unit() == datetime_unit::DAY);
    REQUIRE(datetime::from_string("2002-03").unit() == datetime_unit::MONTH);
    REQUIRE(datetime::from_string("2002").unit() == datetime_unit::YEAR);

    REQUIRE(datetime::from_YmdHMS_digits("2002-03-04 12:13:14").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_YmdHMS_digits("2002-03-04 12:13").unit() == datetime_unit::MINUTE);
    REQUIRE(datetime::from_YmdHMS_digits("2002-03-04 12").unit() == datetime_unit::HOUR);
    REQUIRE(datetime::from_YmdHMS_digits("2002-03-04").unit() == datetime_unit::DAY);
    REQUIRE(datetime::from_YmdHMS_digits("2002-03").unit() == datetime_unit::MONTH);
    REQUIRE(datetime::from_YmdHMS_digits("2002").unit() == datetime_unit::YEAR);

    REQUIRE(duration::from_string("PT3S").dt_unit == datetime_unit::SECOND);
    REQUIRE(duration::from_string("PT3M").dt_unit == datetime_unit::MINUTE);
    REQUIRE(duration::from_string("PT3H").dt_unit == datetime_unit::HOUR);
    REQUIRE(duration::from_string("P3D").dt_unit == datetime_unit::DAY);
    REQUIRE(duration::from_string("P3M").dt_unit == datetime_unit::MONTH);
    REQUIRE(duration::from_string("P3Y").dt_unit == datetime_unit::YEAR);

    REQUIRE(datetime::from_YmdHMS_digits("2021-06-29T10:03:38.494534Z").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_YmdHMS_digits("2021-06-29T10:03:38.494534").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_YmdHMS_digits("2021-06-29T10:03:38.494534+01:00").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_YmdHMS_digits("2021-06-29T10:03:38.494534-01:00").unit() == datetime_unit::SECOND);

    REQUIRE(datetime::from_string("2021-06-29T10:03:38.494534Z").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_string("2021-06-29T10:03:38.494534").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_string("2021-06-29T10:03:38.494534+01:00").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_string("2021-06-29T10:03:38.494534-01:00").unit() == datetime_unit::SECOND);
}

TEST_CASE("Datetime arithmetics", "[datetime]") {
    datetime x = datetime::from_string("2002-03-04");
    datetime y = datetime::from_string("2002-03-10");
    duration d = duration::from_string("P5D");

    REQUIRE((x + d).to_string() == "2002-03-09");
    REQUIRE((x + d).unit() == datetime_unit::DAY);
    REQUIRE((d * 5).to_string() == "P25D");
    REQUIRE((y - x).to_string() == "P6D");
    REQUIRE(x.to_double() == 20020304);

    // TODO test month / year / hours / minutes / seconds / also mixed

    x = datetime::from_string("2002-03");
    y = datetime::from_string("2002-08");
    d = duration::from_string("P3M");

    REQUIRE((x + d).unit() == datetime_unit::MONTH);
    REQUIRE((x + d).to_string() == "2002-06-01");
    REQUIRE((x + d).to_double() == 200206);
}

TEST_CASE("Datetime Access Functions", "[datetime]") {
    datetime x = datetime::from_string("2002-03-04T13:40:21");

    REQUIRE(x.seconds() == 21);
    REQUIRE(x.minutes() == 40);
    REQUIRE(x.hours() == 13);
    REQUIRE(x.dayofmonth() == 04);
    REQUIRE(x.dayofyear() == 63);
    REQUIRE(x.dayofweek() == 1);  // MON
    REQUIRE(x.month() == 03);
    REQUIRE(x.year() == 2002);

    datetime xx = datetime::from_YmdHMS_digits("2002-03-04T13:40:21");
    REQUIRE(xx == x);
    // REQUIRE(x.epoch_time() == 1015249221); // might depend on OS

    datetime y = datetime::from_string("2005-06");

    REQUIRE(y.seconds() == 0);
    REQUIRE(y.minutes() == 0);
    REQUIRE(y.hours() == 0);
    REQUIRE(y.dayofmonth() == 1);
    REQUIRE(y.dayofyear() == 152);
    REQUIRE(y.dayofweek() == 3);  // WED
    REQUIRE(y.month() == 6);
    REQUIRE(y.year() == 2005);

    datetime yy = datetime::from_YmdHMS_digits("2005-06");
    REQUIRE(yy == y);



    REQUIRE(datetime::from_string("2001-03-01").dayofyear() == 60);  // no leap year
    REQUIRE(datetime::from_string("2000-03-01").dayofyear() == 61);  // leap year



    auto xxx = datetime::from_string("2021-01-31") + duration::from_string("P3M");
    REQUIRE(xxx.year() == 2021);
    REQUIRE(xxx.month() == 4);
    REQUIRE(xxx.dayofmonth() == 30);

}