#ifndef REDUCE_TIME_R_H
#define REDUCE_TIME_R_H


#include "gdalcubes/src/cube.h"

using namespace gdalcubes;


/**
* @brief A data cube that applies a user-defined  R function over pixel time-series
*/

class reduce_time_stream_cube : public cube {
public:
  /**
  * @brief Create a data cube that applies a user-defined R function on a given input data cube over time
  * @note This static creation method should preferably be used instead of the constructors as
  * the constructors will not set connections between cubes properly.
  * @param in input data cube
  * @param reducer reducer function
  * @return a shared pointer to the created data cube instance
  */
  static std::shared_ptr<reduce_time_stream_cube> create(std::shared_ptr<cube> in, std::string cmd, uint16_t nbands, std::vector<std::string> names = std::vector<std::string>()) {
    std::shared_ptr<reduce_time_stream_cube> out = std::make_shared<reduce_time_stream_cube>(in, cmd, nbands, names);
    in->add_child_cube(out);
    out->add_parent_cube(in);
    return out;
  }
  
public:
  reduce_time_stream_cube(std::shared_ptr<cube> in, std::string cmd, uint16_t nbands,  std::vector<std::string> names = std::vector<std::string>()) : cube(std::make_shared<cube_st_reference>(*(in->st_reference()))), _in_cube(in), _cmd(cmd), _nbands(nbands), _names(names){  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
    _st_ref->dt((_st_ref->t1() - _st_ref->t0()) + 1);
    _st_ref->t1() = _st_ref->t0();  // set nt=1
    _chunk_size[0] = 1;
    _chunk_size[1] = _in_cube->chunk_size()[1];
    _chunk_size[2] = _in_cube->chunk_size()[2];
    
    if (!names.empty()) {
      if (names.size() != nbands) {
        GCBS_ERROR("size of names is different to nbands");
        throw std::string("ERROR in reduce_time_stream_cube::reduce_time_stream_cube(): size of names is different to nbands");
      }
    }
  
    for (uint16_t i=0; i<nbands; ++i) {
     std::string name;
      if (!_names.empty()) 
        name = _names[i];
      else 
        name = "band" + std::to_string(i+1);
      _bands.add(band(name));
    }
  }
  
public:
  ~reduce_time_stream_cube() {}
  
  std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;
  
 
  nlohmann::json make_constructible_json() override {
    nlohmann::json out;
    out["cube_type"] = "reduce_time_stream_cube";
    out["cmd"] = _cmd ;
    out["nbands"] = _nbands;
    out["names"] = _names;
    out["in_cube"] = _in_cube->make_constructible_json();
    return out;
  }
  
private:
  std::shared_ptr<cube> _in_cube;
  std::string  _cmd;
  uint16_t _nbands;
  std::vector<std::string> _names;
  
  virtual void set_st_reference(std::shared_ptr<cube_st_reference> stref) override {
    // copy fields from st_reference type
    _st_ref->win() = stref->win();
    _st_ref->srs() = stref->srs();
    _st_ref->ny() = stref->ny();
    _st_ref->nx() = stref->nx();
    _st_ref->t0() = stref->t0();
    _st_ref->t1() = stref->t1();
    _st_ref->dt(stref->dt());
    
    _st_ref->dt((_st_ref->t1() - _st_ref->t0()) + 1);
    _st_ref->t1() = _st_ref->t0();  // set nt=1
    //assert(_st_ref->nt() == 1);
  }
};


#endif
