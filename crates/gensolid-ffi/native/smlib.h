/* Gensolid private SMLib ABI. Apache-2.0. Not part of the public GS_ API.
 * Arrays are borrowed for a call; returned arrays belong to owner until free.
 * All entry points are serialized in C++; no C++ exception crosses this ABI.
 */
#pragma once
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
struct GsSmResult {
  void *owner;
  const double *values;
  size_t values_len;
  const size_t *indices;
  size_t indices_len;
  int64_t upstream_status;
};
int gs_sm_create(double linear, double angular, double size, void **session);
void gs_sm_destroy(void *session);
int gs_sm_call(void *session, uint32_t operation, const size_t *ids,
               size_t ids_len, const double *values, size_t values_len,
               const size_t *sizes, size_t sizes_len,
               struct GsSmResult *result);
void gs_sm_result_free(void *owner);
#ifdef __cplusplus
}
#endif
