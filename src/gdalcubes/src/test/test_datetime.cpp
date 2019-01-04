/*
   Copyright 2018 Marius Appel <marius.appel@uni-muenster.de>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include <string>
#include "../datetime.h"
#include "../external/catch.hpp"
TEST_CASE("Deriving datetime unit from string", "[datetime]") {
    REQUIRE(datetime::from_string("2002-03-04 12:13:14").unit() == datetime_unit::SECOND);
    REQUIRE(datetime::from_string("2002-03-04 12:13").unit() == datetime_unit::MINUTE);
    REQUIRE(datetime::from_string("2002-03-04 12").unit() == datetime_unit::HOUR);
    REQUIRE(datetime::from_string("2002-03-04").unit() == datetime_unit::DAY);
    REQUIRE(datetime::from_string("2002-03").unit() == datetime_unit::MONTH);
    REQUIRE(datetime::from_string("2002").unit() == datetime_unit::YEAR);

    REQUIRE(duration::from_string("PT3S").dt_unit == datetime_unit::SECOND);
    REQUIRE(duration::from_string("PT3M").dt_unit == datetime_unit::MINUTE);
    REQUIRE(duration::from_string("PT3H").dt_unit == datetime_unit::HOUR);
    REQUIRE(duration::from_string("P3D").dt_unit == datetime_unit::DAY);
    REQUIRE(duration::from_string("P3M").dt_unit == datetime_unit::MONTH);
    REQUIRE(duration::from_string("P3Y").dt_unit == datetime_unit::YEAR);
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
    REQUIRE((x + d).to_string() == "2002-06");
    REQUIRE((x + d).to_double() == 200206);
}
