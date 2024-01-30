

#ifndef STREAM_REDUCE_TIME_H
#define STREAM_REDUCE_TIME_H

#include "cube.h"

namespace gdalcubes {

/**
* @brief A data cube that applies a user-defined function over pixel time-series
*/

class stream_reduce_time_cube : public cube {
   public:
    /**
        * @brief Create a data cube that applies a user-defined function on a given input data cube over time
        * @note This static creation method should preferably be used instead of the constructors as
        * the constructors will not set connections between cubes properly.
        * @param in input data cube
        * @param cmd external program call
        * @param nbands number of bands in the output cube
        * @param names string vector of output band names, must have the size nbands or empty (the default). If empty,
        * output bands will be named "band1", "band2", ...
        * @return a shared pointer to the created data cube instance
        */
    static std::shared_ptr<stream_reduce_time_cube>
    create(std::shared_ptr<cube> in, std::string cmd, uint16_t nbands,
           std::vector<std::string> names = std::vector<std::string>()) {
        std::shared_ptr<stream_reduce_time_cube> out = std::make_shared<stream_reduce_time_cube>(in, cmd, nbands,
                                                                                                 names);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    stream_reduce_time_cube(std::shared_ptr<cube> in, std::string cmd, uint16_t nbands,
                            std::vector<std::string> names = std::vector<std::string>()) : cube(in->st_reference()->copy()), _in_cube(in), _cmd(cmd), _nbands(nbands), _names(names) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well

        if (cube_stref::type_string(in->st_reference()) == "cube_stref_regular") {
            std::shared_ptr<cube_stref_regular> stref = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref);
            duration dt = (stref->t1() - stref->t0() + 1);
            stref->set_t_axis(stref->t0(), stref->t1(), dt);
        } else if (cube_stref::type_string(_st_ref) == "cube_stref_labeled_time") {
            std::shared_ptr<cube_stref_labeled_time> stref = std::dynamic_pointer_cast<cube_stref_labeled_time>(_st_ref);
            //stref->dt((stref->t1() - stref->t0()) + 1);
            //stref->t1(stref->t0()) ;  // set nt=1
            stref->set_time_labels({stref->t0()});
        }
        _chunk_size[0] = 1;
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (!names.empty()) {
            if (names.size() != nbands) {
                GCBS_ERROR("size of names is different to nbands");
                throw std::string(
                    "ERROR in stream_reduce_time_cube::reduce_time_stream_cube(): size of names is different to nbands");
            }
        }

        for (uint16_t i = 0; i < nbands; ++i) {
            std::string name;
            if (!_names.empty())
                name = _names[i];
            else
                name = "X" + std::to_string(i + 1);
            
            if (!std::isalnum(name[0])) {
              GCBS_WARN("Variable name '" + name  + "' is not compatible with netCDF format; replacing with 'X" + name + "'");
              name = "X" + name;
            }
            
            _bands.add(band(name));
        }
    }

   public:
    ~stream_reduce_time_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "stream_reduce_time";
        out["cmd"] = _cmd;
        out["nbands"] = _nbands;
        out["names"] = _names;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::string _cmd;
    uint16_t _nbands;
    std::vector<std::string> _names;
};

}  // namespace gdalcubes

#endif  //STREAM_REDUCE_TIME_H
