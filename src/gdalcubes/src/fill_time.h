
#ifndef FILL_TIME_H
#define FILL_TIME_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that applies reducer functions over selected bands of a data cube over time
 * @note This is a reimplementation of reduce_cube. The new implementation allows to apply different reducers to different bands instead of just one reducer to all bands of the input data cube
 */
class fill_time_cube : public cube {
   public:
    /**
    * @brief Create a data cube that fills NAN values of a given input data cube by time series interpolation
    * @note This static creation method should preferably be used instead of the constructors as
    * the constructors will not set connections between cubes properly.
    * @param in input data cube
    * @param method interpolation method, currently "near" or "bilinear"
    * @return a shared pointer to the created data cube instance
    */
    static std::shared_ptr<fill_time_cube> create(std::shared_ptr<cube> in, std::string method = "near") {
        std::shared_ptr<fill_time_cube> out = std::make_shared<fill_time_cube>(in, method);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    fill_time_cube(std::shared_ptr<cube> in, std::string method = "near") : cube(in->st_reference()->copy()), _in_cube(in), _method(method) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
            band b = in->bands().get(i);
            _bands.add(b);
        }

        if (!_st_ref->has_regular_time()) {
            GCBS_WARN("Cube has irregular time dimension, interpolation currently uses index-based time distances");
        }

        if (method != "near" && method != "linear" && method != "locf" && method != "nocb") {
            GCBS_WARN("Invalid time-series interpolation method, using default (nearest neighbor)");
            _method = "near";
        }
    }

   public:
    ~fill_time_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "fill_time";
        out["method"] = _method;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::string _method;
};

}  // namespace gdalcubes

#endif  // FILL_TIME_H
