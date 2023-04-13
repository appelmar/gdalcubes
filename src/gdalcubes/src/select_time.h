
#ifndef SELECT_TIME_H
#define SELECT_TIME_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that selects one or more time (irregular) time slices of a data cube, the resulting data cube will have an irregular (labeled) temporal reference
 */
class select_time_cube : public cube {
   public:
    /**
    * @brief Create a data cube that selects one or more time (irregular) time slices of a data cube
    * @note This static creation method should preferably be used instead of the constructors as
    * the constructors will not set connections between cubes properly.
    * @param in input data cube
    * @param t vector of datetime objects
    * @return a shared pointer to the created data cube instance
    */
    static std::shared_ptr<select_time_cube> create(std::shared_ptr<cube> in, std::vector<datetime> t) {
        std::shared_ptr<select_time_cube> out = std::make_shared<select_time_cube>(in, t);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

    /**
    * @brief Create a data cube that selects one or more time (irregular) time slices of a data cube
    * @note This static creation method should preferably be used instead of the constructors as
    * the constructors will not set connections between cubes properly.
    * @param in input data cube
    * @param t vector of datetime strings
    * @return a shared pointer to the created data cube instance
    */
    static std::shared_ptr<select_time_cube> create(std::shared_ptr<cube> in, std::vector<std::string> t) {
        std::vector<datetime> dt;
        for (uint32_t i = 0; i < t.size(); ++i) {
            dt.push_back(datetime::from_string(t[i]));
        }
        std::shared_ptr<select_time_cube> out = std::make_shared<select_time_cube>(in, dt);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    select_time_cube(std::shared_ptr<cube> in, std::vector<datetime> t) : cube(in->st_reference()->copy()), _in_cube(in), _t() {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (t.empty()) {
            GCBS_ERROR("ERROR in select_time_cube::select_time_cube(): missing time slices");
            throw std::string("ERROR: missing time slices in select_time()");
        }

        for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
            band b = in->bands().get(i);
            _bands.add(b);
        }

        std::shared_ptr<cube_stref_labeled_time> stref = std::make_shared<cube_stref_labeled_time>();

        if (!_st_ref->has_regular_space()) {
            throw std::string("ERROR: Time slice selection currently does not support irregular spatial dimensions");
        }

        // NOTE: the following will only work as long as all cube st reference types with regular spatial dimensions inherit from  cube_stref_regular class
        std::shared_ptr<cube_stref_regular> stref_orig = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref);

        stref->set_x_axis(stref_orig->left(), stref_orig->right(), stref_orig->nx());
        stref->set_y_axis(stref_orig->bottom(), stref_orig->top(), stref_orig->ny());

        //stref->dt(stref_orig->dt());

        std::vector<datetime> dt;
        for (uint32_t i = 0; i < t.size(); ++i) {
            if (t[i] >= stref_orig->t0() && t[i] <= stref_orig->t1()) {
                dt.push_back(t[i]);
            }
        }
        if (dt.size() < t.size()) {
            GCBS_WARN("One or more time slices are outside of the input cube and will be ignored");
        }
        if (dt.empty()) {
            GCBS_ERROR("ERROR in select_time_cube::select_time_cube(): resulting cube does not contain any time slice");
            throw std::string("ERROR: resulting cube does not contain any time slice");
        }
        _t = dt;
        stref->set_time_labels(dt);
        _st_ref = stref;
        // TODO: what if t is not sorted
    }

    ~select_time_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        std::vector<std::string> dt;
        for (uint32_t i = 0; i < _t.size(); ++i) {
            dt.push_back(_t[i].to_string());
        }
        json11::Json::object out;
        out["cube_type"] = "select_time";
        out["t"] = dt;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::vector<datetime> _t;
};

}  // namespace gdalcubes

#endif  // SELECT_TIME_H
