/*
    MIT License

    Copyright (c) 2026 Marius Appel <marius.appel@hs-bochum.de>

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
#include "zarr_cube.h"

//#include "filesystem.h"

namespace gdalcubes {


// TODO: Add opening options?

zarr_cube::zarr_cube(std::string path, bool auto_unpack) : cube(), _auto_unpack(auto_unpack), _path({path}), _orig_bands(), _band_selection(), _mutex() {
    if (!filesystem::exists(path)) {
        GCBS_ERROR("Zarr dataset '" + path + "' does not exist");
        throw std::string("Zarr dataset '" + path + "' does not exist.");
    }


    auto dataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(path.c_str(), GDAL_OF_MULTIDIM_RASTER ));
    if( !dataset )
    {
        throw std::string("ERROR: Failed to open Zarr file");
    }
    auto root = dataset->GetRootGroup();
    if( !root )
    {
        throw std::string("ERROR: Failed to read Zarr file");
    }

    // std::shared_ptr<GDALDimension> dim_t, dim_y, dim_x;
    //
    // auto dimensions = root->GetDimensions();
    // for (int i=0; i<dimensions.size(); ++i) {
    //     if (dimensions[i]->GetName() == "x") dim_t = dimensions[i];
    //     if (dimensions[i]->GetName() == "y") dim_y = dimensions[i];
    //     if (dimensions[i]->GetName() == "t") dim_x = dimensions[i];
    // }
    // if (dim_x == nullptr  || dim_y == nullptr || dim_t == nullptr) {
    //     GCBS_ERROR("Failed to identify x,y,t dimensions in Zarr dataset '" + path + "'");
    //     throw std::string("Failed to identify x,y,t dimensions in Zarr dataset '" + path + "'");
    // }
    // uint32_t nt = dim_t->GetSize();
    // uint32_t ny = dim_y->GetSize();
    // uint32_t nx = dim_x->GetSize();



    // Find 3-dimensional arrays as well as arrays for dimension values
    std::string spatial_ref_str;
    std::shared_ptr<GDALMDArray> arr_t, arr_y, arr_x;
    int32_t nt = 0;
    uint32_t ny = 0;
    uint32_t nx = 0;
    auto var_names = root->GetMDArrayNames();
    for (auto vname : var_names) {
        auto arr = root->OpenMDArray(vname);
        auto dims = arr->GetDimensions();

        if (dims.size() == 3) {
            char *wkt_out;
            arr->GetSpatialRef()->exportToPrettyWkt(&wkt_out); // TODO: check if GetSPatialRef must be destroyed...
            spatial_ref_str = wkt_out;
            CPLFree(wkt_out);

            // Assuming the same chunk size for all bands
            auto csize = arr->GetBlockSize();
            _chunk_size[0] = csize[0]; // TODO; How to make sure that the order is t,y,x?
            _chunk_size[1] = csize[1];
            _chunk_size[2] = csize[2];

            // TODO: if no spatial reference, skip
            // TODO: support two-dimensional arrays with dropped time dimension if their length is 1?

            band b(vname);
            b.no_data_value = std::to_string(arr->GetNoDataValueAsDouble());
            b.scale = arr->GetScale();
            b.offset = arr->GetOffset();
            b.unit = arr->GetUnit();
            b.type = arr->GetDataType().GetName();
            _bands.add(b);
            _orig_bands.add(b);
        }
        else if (dims.size() == 1) {
            if (dims[0]->GetName() == "t") {
                arr_t = arr;
                nt = dims[0]->GetSize();
            }
            else if (dims[0]->GetName() == "y") {
                arr_y = arr;
                ny = dims[0]->GetSize();
            }
            else if (dims[0]->GetName() == "x") {
                arr_x= arr;
                nx = dims[0]->GetSize();
            }
        }
    }

    if (_bands.count() == 0) {
        std::string msg = "Zarr dataset does not contain any three-dimensional (t,y,x) arrays";
        GCBS_ERROR(msg);
        throw(msg);
    }



    // Extract geotransform
    std::string geotransform_str = root->GetAttribute("gdalcubes.geotransform")->ReadAsString();
    std::istringstream iss(geotransform_str);
    std::vector<std::string> geotranform_parts(
        std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>());
    if (geotranform_parts.size() != 6) {
        GCBS_ERROR("Failed to parse GeoTransform attribute of Zarr dataset '" + path + "'");
        throw std::string("Failed to parse GeoTransform attribute of Zarr dataset '" + path + "'");
    }

    double left = std::atof(geotranform_parts[0].c_str());
    double dx = std::abs(std::atof(geotranform_parts[1].c_str()));
    double top = std::atof(geotranform_parts[3].c_str());
    double dy = std::abs(std::atof(geotranform_parts[5].c_str()));
    double right = left + dx * nx;
    double bottom = top - dy * ny;


    std::string str_datetime_type = root->GetAttribute("gdalcubes.datetime_type")->ReadAsString();
    std::string str_t0 = root->GetAttribute("gdalcubes.datetime_t0")->ReadAsString();
    std::string str_t1 = root->GetAttribute("gdalcubes.datetime_t1")->ReadAsString();
    std::string str_dt = root->GetAttribute("gdalcubes.datetime_dt")->ReadAsString();



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
        }
        else {
            // 1. Extract time dimension values and check if t is regular anyway
            int32_t *tvalues = (int32_t*)std::malloc(nt * sizeof(int32_t));
            GUInt64 startp[] = {0};
            std::size_t countp[] = {static_cast<std::size_t>(nt)};
            arr_t->Read(startp,countp,nullptr,nullptr, GDALExtendedDataType::Create(GDT_Int32),tvalues);
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

            duration dt = duration::from_string(str_dt);
            datetime t0 = datetime::from_string(str_t0);

            if (time_is_regular) {
                cube_stref_regular ref;
                ref.set_x_axis(left, right, (uint32_t)nx);
                ref.set_y_axis(bottom, top,  (uint32_t)ny);
                ref.srs(spatial_ref_str);
                ref.set_t_axis(datetime::from_string(str_t0), datetime::from_string(str_t1), duration::from_string(str_dt));
                assert(ref.nt() == nt);
                _st_ref = std::make_shared<cube_stref_regular>(ref);
            }
            else {
                cube_stref_labeled_time ref;
                ref.set_x_axis(left, right, (uint32_t)nx);
                ref.set_y_axis(bottom, top, (uint32_t)ny);
                ref.srs(spatial_ref_str);
                std::vector<datetime> labels;
                for (uint32_t i = 0; i < nt; ++i) {
                    labels.push_back(t0 + dt * tvalues[i]);
                }
                ref.set_time_labels(labels);
                _st_ref = std::make_shared<cube_stref_labeled_time>(ref);
            }
            std::free(tvalues);
        }
    }
}

std::shared_ptr<chunk_data> zarr_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("zarr_cube::read_chunk(" + std::to_string(id) + ")");

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
    GUInt64 startp[] = {climits.low[0], climits.low[1], climits.low[2]};
    std::size_t countp[] = {size_btyx[1], size_btyx[2], size_btyx[3]};

    _mutex.lock(); // TODO: needed?

    auto dataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(_path.c_str(), GDAL_OF_MULTIDIM_RASTER ));
    if( !dataset) {
        throw std::string("ERROR: Failed to open Zarr file");
    }
    auto root = dataset->GetRootGroup();
    if( !root ) {
        throw std::string("ERROR: Failed to read Zarr file");
    }

    // TODO: chunk status variables as in netCDF?
    //GCBS_DEBUG("Zarr input file does not contain chunk status data. ");
    out->set_status(chunk_data::chunk_status::UNKNOWN);


    auto f64t = GDALExtendedDataType::Create(GDT_Float64);
    for (uint16_t i = 0; i < size_btyx[0]; ++i) {
        auto arr =  root->OpenMDArray(_bands.get(i).name.c_str());
        if (!arr->Read(startp, countp, nullptr, nullptr, f64t, (&((double *)(out->buf()))[i * size_btyx[1] * size_btyx[2] * size_btyx[3]]))) {
            throw std::string("Failed to read band '" + _bands.get(i).name + "' for chunk " + std::to_string(id) + " from Zarr dataset");
        }
    }
    _mutex.unlock(); // TODO: needed?

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
