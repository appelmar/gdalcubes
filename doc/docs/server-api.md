# REST-like web API for distributed processing

For distributed processing, gdalcubes and gdalcubes_server communicate via a REST-like API (no HATEOAS).
The following documentation describes available endpoints.

By default, `gdalcubes_server` will listen on `http://0.0.0.0:1111/gdalcubes/api/`.



## `GET /version`
Returns version information as text

## `HEAD /file?name=xyz&size=1234`
Checks whether the server instance has file with name and size equal to the given query parameters. Returns
200, if the file is available, 409 if a file with the same name has different size, ...

## `POST /file?name=xyz`
Uploads the file in the application/octet-stream body as file with name as given in the query parameters.

## `POST /cube`
Uploads a JSON description of a cube and returns a cube ID.

## `POST /cube/{cube_id}/{chunk_id}/start`
Queue a chunk read for a given cube and given chunk id as given in the path, returns immediately.

## `GET /cube/{cube_id}/{chunk_id}/status`
Ask for the status of a chunk read request. Possible return values are `notrequested`, `queued`, `running`, `finished`, and `error`.

## `GET /cube/{cube_id}/{chunk_id}/download`
Download a chunk with given id as application/octet-stream. The chunk must have been queued before. If the chunk has not yet been
read, it will block until the data becomes available.