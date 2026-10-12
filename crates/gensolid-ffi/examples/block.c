/* Creates a block through the C API, checks it and prints its vertices. */
#include <stdio.h>

#include "gensolid.h"

#define CHECK(call)                                                   \
  do {                                                                \
    GS_ERROR_code_t err_ = (call);                                    \
    if (err_ != GS_ERROR_no_errors) {                                 \
      fprintf(stderr, "%s failed with error %d\n", #call, (int)err_); \
      return 1;                                                       \
    }                                                                 \
  } while (0)

int main(void) {
  GS_SESSION_start_o_t options;
  CHECK(GS_SESSION_start_o_init(&options));
  CHECK(GS_SESSION_start(&options));

  GS_AXIS2_sf_t basis = {{10.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}};
  GS_BODY_t body = 0;
  CHECK(GS_BODY_create_solid_block(2.0, 4.0, 6.0, &basis, &body));

  int n_faults = 0;
  GS_check_fault_t *faults = NULL;
  CHECK(GS_BODY_check(body, &n_faults, &faults));
  GS_MEMORY_free(faults);

  int n_vertices = 0;
  GS_VERTEX_t *vertices = NULL;
  CHECK(GS_BODY_ask_vertices(body, &n_vertices, &vertices));
  printf("body %d: %d faults, %d vertices\n", body, n_faults, n_vertices);
  for (int i = 0; i < n_vertices; ++i) {
    double p[3];
    CHECK(GS_VERTEX_ask_point(vertices[i], &p));
    printf("  vertex %d: (%g, %g, %g)\n", vertices[i], p[0], p[1], p[2]);
  }
  GS_MEMORY_free(vertices);

  CHECK(GS_SESSION_stop());
  return n_faults == 0 ? 0 : 1;
}
