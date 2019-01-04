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

#include "reduce.h"

std::shared_ptr<chunk_data> reduce_cube::read_chunk(chunkid_t id) {
    GCBS_DEBUG("reduce_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id < 0 || id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    // If input cube is already "reduced", simply return corresponding input chunk
    if (_in_cube->size_t() == 1) {
        return _in_cube->read_chunk(id);
    }

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {_bands.count(), 1, size_tyx[1], size_tyx[2]};
    out->size(size_btyx);

    // Fill buffers accordingly
    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
    double *begin = (double *)out->buf();
    double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    std::fill(begin, end, NAN);

    reducer *r = nullptr;
    if (_reducer == "min") {
        r = new min_reducer();
    } else if (_reducer == "max") {
        r = new max_reducer();
    } else if (_reducer == "mean") {
        r = new mean_reducer();
    } else if (_reducer == "median") {
        r = new median_reducer();
    } else if (_reducer == "sum") {
        r = new sum_reducer();
    } else if (_reducer == "count") {
        r = new count_reducer();
    } else if (_reducer == "prod") {
        r = new prod_reducer();
    } else if (_reducer == "var") {
        r = new var_reducer();
    } else if (_reducer == "sd") {
        r = new sd_reducer();
    } else
        throw std::string("ERROR in reduce_cube::read_chunk(): Unknown reducer given");

    r->init(out);

    // iterate over all chunks that must be read from the input cube to compute this chunk
    for (chunkid_t i = id; i < _in_cube->count_chunks(); i += _in_cube->count_chunks_x() * _in_cube->count_chunks_y()) {
        std::shared_ptr<chunk_data> x = _in_cube->read_chunk(i);
        r->combine(out, x);
    }

    r->finalize(out);
    if (r != nullptr) delete r;

    return out;
}

void reduce_cube::write_gdal_image(std::string path, std::string format, std::vector<std::string> co, std::shared_ptr<chunk_processor> p) {
    std::shared_ptr<progress> prg = config::instance()->get_default_progress_bar()->get();
    prg->set(0);  // explicitly set to zero to show progress bar immediately
    GDALDriver *drv = (GDALDriver *)GDALGetDriverByName(format.c_str());
    if (!drv) {
        throw std::string("ERROR in reduce_cube::write_gdal_image(): Cannot find GDAL driver for given format.");
    }
    // TODO: Check whether driver supports Create()

    CPLStringList out_co;
    for (uint16_t i = 0; i < co.size(); ++i) {
        out_co.AddString(co[i].c_str());
    }

    GDALDataset *gdal_out = drv->Create(path.c_str(), size_x(), size_y(), bands().count(), GDT_Float64, out_co.List());
    if (!gdal_out) {
        throw std::string("ERROR in reduce_cube::write_gdal_image(): cannot create output image");
    }

    OGRSpatialReference proj_out;
    proj_out.SetFromUserInput(_st_ref->proj().c_str());
    char *out_wkt;
    proj_out.exportToWkt(&out_wkt);

    double affine[6];
    affine[0] = _st_ref->win().left;
    affine[3] = _st_ref->win().top;
    affine[1] = _st_ref->dx();
    affine[5] = -_st_ref->dy();
    affine[2] = 0.0;
    affine[4] = 0.0;

    gdal_out->SetProjection(out_wkt);
    gdal_out->SetGeoTransform(affine);
    CPLFree(out_wkt);

    // The following loop seems to be needed for some drivers only
    for (uint16_t b = 0; b < _bands.count(); ++b) {  //            gdal_out->GetRasterBand(b+1)->SetNoDataValue(NAN);
        if (!_bands.get(b).no_data_value.empty()) {
            gdal_out->GetRasterBand(b + 1)->Fill(std::stod(_bands.get(b).no_data_value));
            gdal_out->GetRasterBand(b + 1)->SetNoDataValue(std::stod(_bands.get(b).no_data_value));
        }
    }

    GDALClose(gdal_out);
    std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f = [this, &path, prg](chunkid_t id, std::shared_ptr<chunk_data> dat, std::mutex &m) {
        if (dat->empty()) {
            GCBS_WARN("Output GDAL image contains empty chunk " + std::to_string(id));
            prg->increment((double)1 / (double)this->count_chunks());
            return;
        }
        m.lock();
        GDALDataset *gdal_out = (GDALDataset *)GDALOpen(path.c_str(), GA_Update);
        m.unlock();
        //bounds_nd<uint32_t, 3> cb = chunk_limits(id);
        chunk_coordinate_tyx ct = chunk_coords_from_id(id);
        for (uint16_t b = 0; b < _bands.count(); ++b) {
            uint32_t yoff = std::max(0, ((int)size_y() - ((int)ct[1] + 1) * (int)_chunk_size[1]));
            uint32_t xoff = ct[2] * _chunk_size[2];
            uint32_t xsize = dat->size()[3];
            uint32_t ysize = dat->size()[2];
            m.lock();
            CPLErr res = gdal_out->GetRasterBand(b + 1)->RasterIO(GF_Write, xoff, yoff, xsize,
                                                                  ysize, ((double *)dat->buf()) + b * dat->size()[2] * dat->size()[3], dat->size()[3], dat->size()[2],
                                                                  GDT_Float64, 0, 0, NULL);
            if (res != CE_None) {
                GCBS_WARN("RasterIO (write) failed for " + std::string(gdal_out->GetDescription()));
            }
            m.unlock();
        }
        m.lock();
        GDALClose(gdal_out);
        m.unlock();
        prg->increment((double)1 / (double)this->count_chunks());
    };

    p->apply(shared_from_this(), f);
    prg->finalize();
}