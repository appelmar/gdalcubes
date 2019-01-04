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

#ifndef STREAM_H
#define STREAM_H

#include "cube.h"

/**
 * @brief A data cube that streams data from another data cube to stdin of an external program and captures stdout as result
 *
 * Assumptions for the implementation below:
 * 1. the stream operation will not change the spatial / temporal extent of a chunk
 * 2. except boundary chunks, all result chunks will have the same dimensions
 * 3. the size of the output chunk is a function of the size of an input chunk, the function is linear in each dimension
 *
 * @TODO: implement constructor that takes output size as argument if already known
 */
class stream_cube : public cube {
   public:
    /**
     * @brief Create a data cube that streams chunk of a given input data cube to an external program
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in_cube input data cube
     * @param cmd external program call
     * @param log_output what to to with the output of the external program, either empty, "stdout", "stderr", or a filename
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<stream_cube> create(std::shared_ptr<cube> in_cube, std::string cmd, std::string log_output = "") {
        std::shared_ptr<stream_cube> out = std::make_shared<stream_cube>(in_cube, cmd, log_output);
        in_cube->add_child_cube(out);
        out->add_parent_cube(in_cube);
        return out;
    }

    stream_cube(std::shared_ptr<cube> in_cube, std::string cmd, std::string log_output = "") : cube(std::make_shared<cube_st_reference>(*(in_cube->st_reference()))), _in_cube(in_cube), _cmd(cmd), _log_output(log_output), _keep_input_nt(false), _keep_input_ny(false), _keep_input_nx(false) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        // Test CMD and find out what size comes out.
        cube_size_tyx tmp = _in_cube->chunk_size(0);
        cube_size_btyx csize_in = {_in_cube->bands().count(), tmp[0], tmp[1], tmp[2]};

        // do not read original chunk data (which can be expensive) but simply stream a dummy chunk with proper size here
        std::shared_ptr<chunk_data> dummy_chunk = std::make_shared<chunk_data>();
        dummy_chunk->size(csize_in);
        dummy_chunk->buf(std::calloc(csize_in[0] * csize_in[1] * csize_in[2] * csize_in[3], sizeof(double)));

        std::shared_ptr<chunk_data> c0 = stream_chunk_stdin(dummy_chunk);

        for (uint16_t ib = 0; ib < c0->size()[0]; ++ib) {
            band b("band" + std::to_string(ib + 1));
            b.no_data_value = std::to_string(NAN);
            b.type = "float64";
            _bands.add(b);
        }

        _chunk_size = {c0->size()[1], c0->size()[2], c0->size()[3]};

        // Make an optimistic guess
        if (c0->size()[1] == 1) {
            _st_ref->nt(in_cube->count_chunks_t());
        } else if (c0->size()[1] == csize_in[1]) {
            _keep_input_nt = true;
            _st_ref->nt(in_cube->size()[1]);
        } else
            throw std::string("ERROR in stream_cube::stream_cube(): could not derive size of result cube");

        if (c0->size()[2] == 1) {
            _st_ref->ny() = in_cube->count_chunks_y();
        } else if (c0->size()[2] == csize_in[2]) {
            _keep_input_ny = true;
            _st_ref->ny() = in_cube->size()[2];
        } else
            throw std::string("ERROR in stream_cube::stream_cube(): could not derive size of result cube");

        if (c0->size()[3] == 1) {
            _st_ref->nx() = in_cube->count_chunks_x();
        } else if (c0->size()[3] == csize_in[3]) {
            _keep_input_nx = true;
            _st_ref->nx() = in_cube->size()[3];
        } else
            throw std::string("ERROR in stream_cube::stream_cube(): could not derive size of result cube");
    }

    virtual std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    virtual nlohmann::json make_constructible_json() override {
        nlohmann::json out;
        out["cube_type"] = "stream";
        out["command"] = _cmd;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::string _cmd;
    std::string _log_output;

    // Variables to help deriving the size when view changes without testing with a dummy chunk
    bool _keep_input_nt;
    bool _keep_input_ny;
    bool _keep_input_nx;

   private:
    std::shared_ptr<chunk_data> stream_chunk_stdin(std::shared_ptr<chunk_data> data);
    virtual void set_st_reference(std::shared_ptr<cube_st_reference> stref) override {
        _st_ref->win() = stref->win();
        _st_ref->proj() = stref->proj();
        _st_ref->ny() = stref->ny();
        _st_ref->nx() = stref->nx();
        _st_ref->t0() = stref->t0();
        _st_ref->t1() = stref->t1();
        _st_ref->dt() = stref->dt();

        if (!_keep_input_nt) _st_ref->nt(count_chunks_t());
        if (!_keep_input_ny) _st_ref->ny() = count_chunks_y();
        if (!_keep_input_nx) _st_ref->nx() = count_chunks_x();
    }
};

#endif  //STREAM_H
