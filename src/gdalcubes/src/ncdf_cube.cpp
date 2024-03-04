/*
    MIT License

    Copyright (c) 2021 Marius Appel <marius.appel@hs-bochum.de>

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
#include "ncdf_cube.h"

#include <netcdf.h>

#include "filesystem.h"

namespace gdalcubes {

// Helper function to convert ncdf text attributes to std::string
std::string ncdf_attr_to_string(int ncid, int varid, std::string attr_name) {
    std::size_t len = 0;
    char *attr_text = nullptr;

    int retval = nc_inq_attlen(ncid, varid, attr_name.c_str(), &len);
    if (retval != NC_NOERR) {
        GCBS_DEBUG("Failed to find attribute '" + attr_name + "' of variable '" + std::to_string(varid) + "' in netCDF file");
        return "";
    }
    attr_text = (char *)std::malloc(len + 1);
    retval = nc_get_att_text(ncid, varid, attr_name.c_str(), attr_text);
    attr_text[len] = '\0';
    if (retval != NC_NOERR) {
        std::free(attr_text);
        GCBS_DEBUG("Failed to read attribute '" + attr_name + "' of variable '" + std::to_string(varid) + "' in netCDF file");
        return "";
    }
    std::string out = attr_text;
    std::free(attr_text);
    return out;
}

ncdf_cube::ncdf_cube(std::string path, bool auto_unpack) : cube(), _auto_unpack(auto_unpack), _path({path}), _orig_bands(), _band_selection(), _mutex() {
    if (!filesystem::is_regular_file(path)) {
        GCBS_ERROR("NetCDF file '" + path + "' does not exist or is not a file");
        throw std::string("NetCDF file '" + path + "' does not exist or is not a file");
    }

    // Open file
    int ncfile;
    int retval = nc_open(path.c_str(), NC_NOWRITE, &ncfile);
    if (retval != NC_NOERR) {
        GCBS_ERROR("Failed to open netCDF file '" + path + "'; nc_open() returned " + std::to_string(retval));
        throw std::string("Failed to open netCDF file '" + path + "'; nc_open() returned " + std::to_string(retval));
    }

    // find dimensions
    int ndims = 0;
    retval = nc_inq_ndims(ncfile, &ndims);
    if (retval != NC_NOERR) {
        GCBS_ERROR("Failed to read dimension number of netCDF file '" + path + "'; nc_inq_ndims() returned " + std::to_string(retval));
        throw std::string("Failed to read dimension number of netCDF file '" + path + "'; nc_inq_ndims() returned " + std::to_string(retval));
    }

    int dim_id_x = -1;
    int dim_id_y = -1;
    int dim_id_t = -1;

    retval = nc_inq_dimid(ncfile, "x", &dim_id_x);
    if (retval != NC_NOERR) {
        retval = nc_inq_dimid(ncfile, "longitude", &dim_id_x);
    }

    retval = nc_inq_dimid(ncfile, "y", &dim_id_y);
    if (retval != NC_NOERR) {
        retval = nc_inq_dimid(ncfile, "latitude", &dim_id_y);
    }

    retval = nc_inq_dimid(ncfile, "time", &dim_id_t);
    if (retval != NC_NOERR) {
        retval = nc_inq_dimid(ncfile, "t", &dim_id_t);
    }

    if (dim_id_x < 0 || dim_id_x >= ndims ||
        dim_id_y < 0 || dim_id_y >= ndims ||
        dim_id_t < 0 || dim_id_t >= ndims) {
        GCBS_ERROR("Failed to identify x,y,t dimensions in netCDF file '" + path + "'");
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("Failed to identify x,y,t dimensions in netCDF file '" + path + "'");
    }

    int var_id_x = -1;
    int var_id_y = -1;
    int var_id_t = -1;

    retval = nc_inq_varid(ncfile, "x", &var_id_x);
    if (retval != NC_NOERR) {
        retval = nc_inq_varid(ncfile, "longitude", &var_id_x);
    }

    retval = nc_inq_varid(ncfile, "y", &var_id_y);
    if (retval != NC_NOERR) {
        retval = nc_inq_varid(ncfile, "latitude", &var_id_y);
    }

    retval = nc_inq_varid(ncfile, "time", &var_id_t);
    if (retval != NC_NOERR) {
        retval = nc_inq_varid(ncfile, "t", &var_id_t);
    }

    if (var_id_x < 0 || var_id_y < 0 || var_id_t < 0) {
        GCBS_ERROR("Failed to identify x,y,t variables in netCDF file '" + path + "'");
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("Failed to identify x,y,t variables in netCDF file '" + path + "'");
    }

    /* Read SRS and Geotransform metadata */
    int v_crs = -1;
    retval = nc_inq_varid(ncfile, "crs", &v_crs);
    if (retval != NC_NOERR) {
        GCBS_ERROR("Failed to find crs variable in netCDF file '" + path + "'");
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("Failed to find crs variable in netCDF file '" + path + "'");
    }

    std::size_t len = 0;
    char *attr_text = nullptr;

    retval = nc_inq_attlen(ncfile, v_crs, "spatial_ref", &len);
    if (retval != NC_NOERR) {
        GCBS_ERROR("Failed to find spatial_ref attribute of crs variable in netCDF file '" + path + "'");
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("Failed to find spatial_ref attribute of crs variable in netCDF file '" + path + "'");
    }
    attr_text = (char *)std::malloc(len + 1);
    retval = nc_get_att_text(ncfile, v_crs, "spatial_ref", attr_text);
    attr_text[len] = '\0';
    if (retval != NC_NOERR) {
        GCBS_ERROR("Failed to read spatial_ref attribute of crs variable in netCDF file '" + path + "'");
        std::free(attr_text);
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("Failed to read spatial_ref attribute of crs variable in netCDF file '" + path + "'");
    }
    std::string spatial_ref_str = attr_text;
    std::free(attr_text);

    retval = nc_inq_attlen(ncfile, v_crs, "GeoTransform", &len);
    if (retval != NC_NOERR) {
        GCBS_ERROR("Failed to find GeoTransform attribute of crs variable in netCDF file '" + path + "'");
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("Failed to find GeoTransform attribute of crs variable in netCDF file '" + path + "'");
    }
    attr_text = (char *)std::malloc(len + 1);
    retval = nc_get_att_text(ncfile, v_crs, "GeoTransform", attr_text);
    attr_text[len] = '\0';
    if (retval != NC_NOERR) {
        GCBS_ERROR("Failed to read GeoTransform attribute of crs variable in netCDF file '" + path + "'");
        std::free(attr_text);
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("Failed to read GeoTransform attribute of crs variable in netCDF file '" + path + "'");
    }
    std::string geotransform_str = attr_text;
    std::free(attr_text);

    // TODO: parse geotransform
    std::istringstream iss(geotransform_str);
    std::vector<std::string> geotranform_parts(
        std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>());
    if (geotranform_parts.size() != 6) {
        GCBS_ERROR("Failed to parse GeoTransform attribute of crs variable in netCDF file '" + path + "'");
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("Failed to parse GeoTransform attribute of crs variable in netCDF file '" + path + "'");
    }

    double left = std::atof(geotranform_parts[0].c_str());
    double dx = std::abs(std::atof(geotranform_parts[1].c_str()));
    double top = std::atof(geotranform_parts[3].c_str());
    double dy = std::abs(std::atof(geotranform_parts[5].c_str()));

    // find out nx / ny
    std::size_t nx = -1;
    std::size_t ny = -1;
    retval = nc_inq_dimlen(ncfile, dim_id_x, &nx);  // TODO: error handling
    retval = nc_inq_dimlen(ncfile, dim_id_y, &ny);  // TODO: error handling

    double right = left + dx * nx;
    double bottom = top - dy * ny;

    /* Read time metadata */

    // 0. try to read from NC_GLOBAL metadata (gdalcubes version >= 0.3.2)

    std::string str_datetime_type = ncdf_attr_to_string(ncfile, NC_GLOBAL, "gdalcubes_datetime_type");
    std::string str_t0 = ncdf_attr_to_string(ncfile, NC_GLOBAL, "gdalcubes_datetime_t0");
    std::string str_t1 = ncdf_attr_to_string(ncfile, NC_GLOBAL, "gdalcubes_datetime_t1");
    std::string str_dt = ncdf_attr_to_string(ncfile, NC_GLOBAL, "gdalcubes_datetime_dt");

    bool finished_st_reference = false;
    if (!str_datetime_type.empty() &&
        !str_t0.empty() &&
        !str_t1.empty() &&
        !str_dt.empty()) {
        if (str_datetime_type == "regular") {
            cube_stref_regular ref;
            ref.set_x_axis(left, right, (uint32_t)nx);
            ref.set_y_axis(bottom, top, (uint32_t)ny);
            ref.srs(spatial_ref_str);
            ref.set_t_axis(datetime::from_string(str_t0), datetime::from_string(str_t1), duration::from_string(str_dt));
            _st_ref = std::make_shared<cube_stref_regular>(ref);
            finished_st_reference = true;
        }
    }

    if (!finished_st_reference) {  // gdalcubes version < 0.3.2 / labeled / previous error
        // 1. Find out nt using nc_inq_dimlen(), allocate int buffer
        // 2. Read all integer time values using nc_get_var_int()
        // 3. Iterate over integer values and find out whether time is regular or not! and create new st_ref accordingly
        // 4. Get datetime unit and start with nc_inq_attlen() nc_get_att_text()  [var = id of time, ttribute name = "units"]
        // 5. free int buffer
        std::size_t nt = -1;
        retval = nc_inq_dimlen(ncfile, dim_id_t, &nt);  // TODO: error handling
        int *tvalues = (int *)std::malloc(nt * sizeof(int));
        nc_get_var_int(ncfile, var_id_t, tvalues);  // TODO: error handling
        bool time_is_regular = true;
        int delta = tvalues[0];
        if (nt > 1) {
            delta = tvalues[1] - tvalues[0];
            for (uint32_t i = 2; i < nt; ++i) {
                if (tvalues[i] - tvalues[i - 1] != delta) {
                    time_is_regular = false;
                    break;
                }
            }
        }

        retval = nc_inq_attlen(ncfile, var_id_t, "units", &len);  // TODO: error handling
        attr_text = (char *)std::malloc(len + 1);
        retval = nc_get_att_text(ncfile, var_id_t, "units", attr_text);  // TODO: error handling
        attr_text[len] = '\0';
        std::string datetime_str = attr_text;
        std::free(attr_text);

        // find out t0 and dt unit from datetime_str
        std::istringstream isd(datetime_str);
        std::vector<std::string> datetime_parts(
            std::istream_iterator<std::string>{isd},
            std::istream_iterator<std::string>());
        if (datetime_parts.size() != 3) {
            GCBS_ERROR("Failed to parse time units in netCDF file '" + path + "'");
            retval = nc_close(ncfile);
            if (retval != NC_NOERR) {
                GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
            }
            throw std::string("Failed to parse time units in netCDF file '" + path + "'");
        }

        datetime_unit tunit;
        if (datetime_parts[0] == "years") {
            tunit = datetime_unit::YEAR;
        } else if (datetime_parts[0] == "months") {
            tunit = datetime_unit::MONTH;
        } else if (datetime_parts[0] == "days") {
            tunit = datetime_unit::DAY;
        } else if (datetime_parts[0] == "hours") {
            tunit = datetime_unit::HOUR;
        } else if (datetime_parts[0] == "minutes") {
            tunit = datetime_unit::MINUTE;
        } else if (datetime_parts[0] == "seconds") {
            tunit = datetime_unit::SECOND;
        } else {
            GCBS_ERROR("Failed to parse datetime unit '" + datetime_parts[0] + "' in netCDF file '" + path + "'");
            retval = nc_close(ncfile);
            if (retval != NC_NOERR) {
                GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
            }
            throw std::string("Failed to parse datetime unit '" + datetime_parts[0] + "' in netCDF file '" + path + "'");
        }

        duration dt;
        if (time_is_regular) {
            dt.dt_interval = delta;
            dt.dt_unit = tunit;
        } else {
            dt.dt_interval = 1;
            dt.dt_unit = tunit;
        }

        datetime t0 = datetime::from_string(datetime_parts[2]);
        if (time_is_regular) {
            datetime t1 = t0 + duration(tvalues[nt - 1] + delta - 1, tunit);
            // If nt == 1, delta cannot be derived and hence dt interval is unclear
            if (delta == 0) {
                t1 = t0;
                dt.dt_interval = 1;
                GCBS_WARN("Setting dt = " + dt.to_string() + " due to missing metadata in netCDF file");
            }
            cube_stref_regular ref;
            ref.set_x_axis(left, right, (uint32_t)nx);
            ref.set_y_axis(bottom, top,  (uint32_t)ny);
            ref.srs(spatial_ref_str);
            ref.set_t_axis(datetime::from_string(str_t0), datetime::from_string(str_t1), duration::from_string(str_dt));
            assert(ref.nt() == nt);
            //ref.nt(nt);
            _st_ref = std::make_shared<cube_stref_regular>(ref);
        } else {
            cube_stref_labeled_time ref;
            ref.set_x_axis(left, right, (uint32_t)nx);
            ref.set_y_axis(bottom, top, (uint32_t)ny);
            ref.srs(spatial_ref_str);

            //ref.dt(dt);
            std::vector<datetime> labels;
            for (uint32_t i = 0; i < nt; ++i) {
                labels.push_back(t0 + dt * tvalues[i]);
            }
            ref.set_time_labels(labels);
            _st_ref = std::make_shared<cube_stref_labeled_time>(ref);
        }
        std::free(tvalues);
    }

    /* Find band variables and metadata (having x,y,t dimensions) */

    int nvars = -1;
    retval = nc_inq_nvars(ncfile, &nvars);  // TODO: error handling

    if (nvars < 1) {
        GCBS_ERROR("No variables found in netCDF file '" + path + "'");
        retval = nc_close(ncfile);
        if (retval != NC_NOERR) {
            GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
        }
        throw std::string("No variables found in netCDF file '" + path + "'");
    } else {
        for (uint16_t i = 0; i < nvars; ++i) {
            if (i == var_id_t ||
                i == var_id_y ||
                i == var_id_x ||
                i == v_crs) continue;
            char *name = (char *)std::malloc(NC_MAX_NAME * sizeof(char));
            retval = nc_inq_varname(ncfile, i, name);  // TODO: error handling
            std::string varname = name;
            std::free(name);
            int ndims = -1;
            retval = nc_inq_varndims(ncfile, i, &ndims);  // TODO: error handling
            if (ndims != 3) {
                GCBS_DEBUG("Variable '" + varname + "' in netCDF file '" + path + "' has " + std::to_string(ndims) + " dimensions and will be skipped");
                continue;
            }
            int *dims = (int *)std::malloc(sizeof(int) * ndims);
            retval = nc_inq_vardimid(ncfile, i, dims);  // TODO: error handling

            if (dims[0] == dim_id_t &&
                dims[1] == dim_id_y &&
                dims[2] == dim_id_x) {
                band b(varname);
                double v;
                if (nc_get_att_double(ncfile, i, "_FillValue", &v) == NC_NOERR) {
                    b.no_data_value = std::to_string(v);
                }
                if (nc_get_att_double(ncfile, i, "scale_factor", &v) == NC_NOERR) {
                    b.scale = v;
                }
                if (nc_get_att_double(ncfile, i, "add_offset", &v) == NC_NOERR) {
                    b.offset = v;
                }
                if (nc_inq_attlen(ncfile, v_crs, "type", &len) == NC_NOERR) {
                    attr_text = (char *)std::malloc(len + 1);
                    if (nc_get_att_text(ncfile, v_crs, "type", attr_text) == NC_NOERR) {
                        attr_text[len] = '\0';
                        b.type = attr_text;
                    }
                    std::free(attr_text);
                }
                if (nc_inq_attlen(ncfile, v_crs, "units", &len) == NC_NOERR) {
                    attr_text = (char *)std::malloc(len + 1);
                    if (nc_get_att_text(ncfile, v_crs, "units", attr_text) == NC_NOERR) {
                        attr_text[len] = '\0';
                        b.unit = attr_text;
                    }
                    std::free(attr_text);
                }
                _bands.add(b);
                _orig_bands.add(b);

                // if available use chunk size of netCDF file
                // Assumption here is that all band variables have the same chunk sizes
                std::size_t chunksize_in[3];
                int storage_in;
                if (nc_inq_var_chunking(ncfile, i, &storage_in, chunksize_in) == NC_NOERR) {
                    if (storage_in == NC_CHUNKED) {
                        _chunk_size[0] = chunksize_in[0];
                        _chunk_size[1] = chunksize_in[1];
                        _chunk_size[2] = chunksize_in[2];
                    }
                }
            } else {
                GCBS_DEBUG("Variable '" + varname + "' in netCDF file '" + path + "' will be skipped; dimensions do not match (t/y/x)");
            }
            std::free(dims);
        }
    }

    std::string process_graph;  // TODO: store as member variable
    if (nc_inq_attlen(ncfile, NC_GLOBAL, "process_graph", &len) == NC_NOERR) {
        attr_text = (char *)std::malloc(len + 1);
        if (nc_get_att_text(ncfile, v_crs, "GeoTransform", attr_text) == NC_NOERR) {
            attr_text[len] = '\0';
            process_graph = attr_text;
        }
        std::free(attr_text);
    }

    retval = nc_close(ncfile);
    if (retval != NC_NOERR) {
        GCBS_DEBUG("Failed to properly close netCDF file '" + path + "'; nc_close() returned " + std::to_string(retval));
    }
}

std::shared_ptr<chunk_data> ncdf_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("ncdf_cube::read_chunk(" + std::to_string(id) + ")");

    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks()) {
        // chunk is outside of the cube, we don't need to read anything.
        GCBS_DEBUG("Chunk id " + std::to_string(id) + " is out of range");
        return out;
    }

    // Derive how many pixels the chunk has (this varies for chunks at the boundary)
    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {_bands.count(), size_tyx[0], size_tyx[1], size_tyx[2]};
    out->size(size_btyx);

    if (size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3] == 0)
        return out;

    // Fill buffers accordingly
    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
    double *begin = (double *)out->buf();
    double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    std::fill(begin, end, NAN);

    bounds_nd<uint32_t, 3> climits = chunk_limits(id);
    std::size_t startp[] = {climits.low[0], climits.low[1], climits.low[2]};
    std::size_t countp[] = {size_btyx[1], size_btyx[2], size_btyx[3]};

    _mutex.lock();
    // Open file
    int ncfile;
    int retval = nc_open(_path.c_str(), NC_NOWRITE, &ncfile);
    if (retval != NC_NOERR) {
        GCBS_ERROR("Failed to open netCDF file '" + _path + "'; nc_open() returned " + std::to_string(retval));
        throw std::string("Failed to open netCDF file '" + _path + "'; nc_open() returned " + std::to_string(retval));
    }

    int varid = -1;
    if (nc_inq_varid(ncfile, "chunk_status", &varid) == NC_NOERR) {
        int s = 0;
        std::size_t nc_chunk_id = std::size_t(id);
        nc_get_var1_int(ncfile, varid, &nc_chunk_id, &s);
        out->set_status(static_cast<chunk_data::chunk_status>(s));
    }
    else {
        GCBS_DEBUG("NetCDF input file does not contain chunk status data. ");
        out->set_status(chunk_data::chunk_status::UNKNOWN);
    }

    for (uint16_t i = 0; i < size_btyx[0]; ++i) {
        if (nc_inq_varid(ncfile, _bands.get(i).name.c_str(), &varid) == NC_NOERR) {
            nc_get_vara_double(ncfile, varid, startp, countp, (&((double *)(out->buf()))[i * size_btyx[1] * size_btyx[2] * size_btyx[3]]));
        } else {
            GCBS_ERROR("Failed to read band '" + _bands.get(i).name + "' for chunk " + std::to_string(id) + " from netCDF file");
            throw std::string("Failed to read band '" + _bands.get(i).name + "' for chunk " + std::to_string(id) + " from netCDF file");
        }
    }

    retval = nc_close(ncfile);
    if (retval != NC_NOERR) {
        GCBS_DEBUG("Failed to properly close netCDF file '" + _path + "'; nc_close() returned " + std::to_string(retval));
    }

    _mutex.unlock();

    double nodata = NAN;
    if (_auto_unpack) {
        for (uint16_t ib = 0; ib < size_btyx[0]; ++ib) {
            if (!_bands.get(ib).no_data_value.empty()) {
                nodata = std::atof(_bands.get(ib).no_data_value.c_str());  // TODO ignore if this fails
            }
            for (uint32_t i = 0; i < size_btyx[1] * size_btyx[2] * size_btyx[3]; ++i) {
                double &v = ((double *)out->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3] + i];
                if (v == nodata) {
                    v = NAN;
                } else {
                    v = _bands.get(ib).offset + v * _bands.get(ib).scale;
                }
            }
        }
    }

    // check if chunk is completely NAN and if yes, return empty chunk
    if (out->all_nan()) {
        auto s = out->status();
        out = std::make_shared<chunk_data>();
        out->set_status(s);
    }
    return out;
}

}  // namespace gdalcubes
