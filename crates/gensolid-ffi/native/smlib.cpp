// Gensolid SMLib adapter, Apache-2.0. Original glue, not an algorithm port.
// Upstream: NVIDIA-Omniverse/usd-brep @
// c9c5979f04cea8c38213004d4f04df0dcfb40d05.
#include "smlib.h"
#include <NL_Nurbsdef.h>
#include <SmApiCurves.h>
#include <SmApiGeneral.h>
#include <SmApiSurfaces.h>
#include <SmNurbsSrf.h>
#include <SmSmlibAll.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <vector>

#ifdef __EMSCRIPTEN__
// The JS transport copies this wasm32 layout; fail compilation on ABI drift.
static_assert(sizeof(GsSmResult) == 32);
static_assert(offsetof(GsSmResult, values) == 4);
static_assert(offsetof(GsSmResult, values_len) == 8);
static_assert(offsetof(GsSmResult, indices) == 12);
static_assert(offsetof(GsSmResult, indices_len) == 16);
static_assert(offsetof(GsSmResult, upstream_status) == 24);
static_assert(sizeof(size_t) == 4);
#endif

namespace {
std::mutex gate;
struct Failure {
  int code;
  SmStatus native;
};
void require(bool valid) {
  if (!valid)
    throw Failure{1, 0};
}
void checked(SmStatus status) {
  if (status == SM_SUCCESS)
    return;
  int code = 10;
  if (status == SM_ERR_INVALID_INPUT)
    code = 1;
  else if (status == SM_ERR_DEGENERATE_SURFACE)
    code = 4;
  else if (status == SM_ERR_OUT_OF_MEMORY)
    code = 8;
  else if (status == SM_ERR_NOT_WITHIN_TOLERANCE)
    code = 6;
  throw Failure{code, status};
}
struct FreeNurb {
  void operator()(void *p) const { smos_Free(p); }
};
struct Session {
  double linear, angular;
  SmContext context;
  std::vector<std::unique_ptr<SmObject>> objects;
  Session(double l, double a, double size)
      : linear(l), angular(a), context(size) {}
  template <class T> T *get(size_t id) {
    require(id < objects.size());
    auto p = dynamic_cast<T *>(objects[id].get());
    if (!p)
      throw Failure{3, 0};
    return p;
  }
};
struct Response {
  std::vector<double> values;
  std::vector<size_t> indices;
};
void point(Response &out, const SmVector3d &p) {
  out.values.insert(out.values.end(), {p.x, p.y, p.z});
}
void control(Response &out, const gw_CPOINT &p) {
  bool rational = p.w != NL_NOW;
  const double w = rational ? p.w : 1.0;
  require(std::isfinite(w) && w > 0);
  out.values.insert(out.values.end(),
                    {p.x / w, p.y / w, p.z == NL_NOZ ? 0.0 : p.z / w, w});
}
void add(Session &s, Response &out, std::unique_ptr<SmObject> p) {
  require(bool(p));
  p->SetContext(&s.context);
  out.indices.push_back(s.objects.size());
  s.objects.push_back(std::move(p));
}
template <class T> struct DeleteOnExit {
  T *&pointer;
  explicit DeleteOnExit(T *&p) : pointer(p) {}
  ~DeleteOnExit() { delete pointer; }
};
template <class T>
void result(Session &s, Response &out, SmStatus status, T *&p) {
  std::unique_ptr<SmObject> owner(p);
  p = nullptr;
  checked(status);
  add(s, out, std::move(owner));
}
template <class T> struct ArrayOwner {
  SmTArray<T *> array;
  ~ArrayOwner() {
    for (ULONG i = 0; i < array.GetSize(); ++i)
      delete array[i];
  }
  void append(std::unique_ptr<T> object) {
    array.Add(object.get());
    object.release();
  }
  void take(Session &s, Response &out) {
    for (ULONG i = 0; i < array.GetSize(); ++i) {
      std::unique_ptr<SmObject> p(array[i]);
      array[i] = nullptr;
      add(s, out, std::move(p));
    }
  }
};
std::unique_ptr<SmBSplineCurve> copy_curve(Session &s, SmBSplineCurve &source) {
  if (!source.GetOrCreateGwNurbPointer())
    throw Failure{6, 0};
  return std::unique_ptr<SmBSplineCurve>(new (s.context)
                                             SmBSplineCurve(source));
}
std::unique_ptr<SmBSplineSurface> copy_surface(Session &s,
                                               SmBSplineSurface &source) {
  if (!source.GetOrCreateGwNurbPointer())
    throw Failure{6, 0};
  return std::unique_ptr<SmBSplineSurface>(new (s.context)
                                               SmBSplineSurface(source));
}
void solutions(Response &out, SmSolutionArray &a, size_t dim) {
  for (ULONG i = 0; i < a.GetSize(); ++i) {
    const auto &x = a[i];
    out.indices.push_back(static_cast<size_t>(x.m_eSolutionType));
    for (size_t j = 0; j < dim; ++j)
      out.values.push_back(x.m_vStart[j]);
    for (size_t j = 0; j < dim; ++j)
      out.values.push_back(x.m_eSolutionType == SM_ST_RANGE_OF_VALUES
                               ? x.m_vEnd[j]
                               : x.m_vStart[j]);
  }
}
// Private ABI: the Rust facade supplies validated dimensions and owns all
// spans.
void dispatch(Session &s, uint32_t op, const std::vector<size_t> &ids,
              const std::vector<double> &d, const std::vector<size_t> &n,
              Response &out) {
  auto c = [&](size_t i) { return s.get<SmCurve>(ids.at(i)); };
  auto b = [&](size_t i) { return s.get<SmBSplineCurve>(ids.at(i)); };
  auto f = [&](size_t i) { return s.get<SmSurface>(ids.at(i)); };
  auto v = [&](size_t i) {
    return SmVector3d(d.at(i), d.at(i + 1), d.at(i + 2));
  };
  auto pts = [&]() {
    require(d.size() % 3 == 0);
    SmTArray<SmPoint3d> p;
    for (size_t i = 0; i < d.size(); i += 3)
      p.Add(v(i));
    return p;
  };
  SmBSplineCurve *bc = nullptr;
  SmCurve *curve = nullptr;
  SmSurface *surface = nullptr;
  DeleteOnExit<SmBSplineCurve> bc_guard(bc);
  DeleteOnExit<SmCurve> curve_guard(curve);
  DeleteOnExit<SmSurface> surface_guard(surface);
  switch (op) {
  case 1: {
    size_t degree = n.at(0), count = n.at(1), knots = n.at(2);
    require(degree > 0 && degree <= SHRT_MAX && count > degree &&
            knots == count + degree + 1 && d.size() == knots + 4 * count);
    std::unique_ptr<void, FreeNurb> memory(
        sm_AllocateNurbCurve(count - 1, degree, knots - 1));
    if (!memory)
      throw std::bad_alloc();
    auto p = static_cast<gw_CURVE *>(memory.get());
    std::copy_n(d.data(), knots, p->knt->U);
    for (size_t i = 0; i < count; ++i) {
      size_t j = knots + 4 * i;
      double w = d[j + 3];
      require(w > 0 && std::isfinite(d[j] * w) && std::isfinite(d[j + 1] * w) &&
              std::isfinite(d[j + 2] * w));
      p->pol->Pw[i] = {d[j] * w, d[j + 1] * w, d[j + 2] * w, w};
    }
    auto object = std::unique_ptr<SmBSplineCurve>(
        new (s.context) SmBSplineCurve(p, 3, FALSE, &s.context));
    memory.release();
    if (op == 1 && !std::equal(d.begin(), d.begin() + knots, p->knt->U))
      throw Failure{6, 0};
    add(s, out, std::move(object));
    break;
  }
  case 2: {
    size_t pu = n.at(0), pv = n.at(1), nu = n.at(2), nv = n.at(3), ku = n.at(4),
           kv = n.at(5);
    require(pu > 0 && pv > 0 && pu <= SHRT_MAX && pv <= SHRT_MAX && nu > pu &&
            nv > pv && ku == nu + pu + 1 && kv == nv + pv + 1 &&
            d.size() == ku + kv + 4 * nu * nv);
    std::unique_ptr<void, FreeNurb> memory(
        sm_AllocateNurbSurface(nu - 1, nv - 1, pu, pv, ku - 1, kv - 1));
    if (!memory)
      throw std::bad_alloc();
    auto p = static_cast<gw_SURFACE *>(memory.get());
    std::copy_n(d.data(), ku, p->knu->U);
    std::copy_n(d.data() + ku, kv, p->knv->U);
    for (size_t i = 0; i < nu; ++i)
      for (size_t j = 0; j < nv; ++j) {
        size_t k = ku + kv + 4 * (i * nv + j);
        double w = d[k + 3];
        require(w > 0 && std::isfinite(d[k] * w) &&
                std::isfinite(d[k + 1] * w) && std::isfinite(d[k + 2] * w));
        p->net->Pw[i][j] = {d[k] * w, d[k + 1] * w, d[k + 2] * w, w};
      }
    auto object = std::unique_ptr<SmBSplineSurface>(
        new (s.context) SmBSplineSurface(p, FALSE, &s.context));
    memory.release();
    if (!std::equal(d.begin(), d.begin() + ku, p->knu->U) ||
        !std::equal(d.begin() + ku, d.begin() + ku + kv, p->knv->U))
      throw Failure{6, 0};
    add(s, out, std::move(object));
    break;
  }
  case 3: {
    auto p = b(0)->GetOrCreateGwNurbPointer();
    require(p);
    out.indices = {size_t(p->p), size_t(p->pol->n + 1), size_t(p->knt->m + 1)};
    out.values.assign(p->knt->U, p->knt->U + p->knt->m + 1);
    for (long i = 0; i <= p->pol->n; ++i)
      control(out, p->pol->Pw[i]);
    break;
  }
  case 4: {
    auto p = s.get<SmBSplineSurface>(ids.at(0))->GetOrCreateGwNurbPointer();
    require(p);
    out.indices = {size_t(p->p),          size_t(p->q),
                   size_t(p->net->n + 1), size_t(p->net->m + 1),
                   size_t(p->knu->m + 1), size_t(p->knv->m + 1)};
    out.values.assign(p->knu->U, p->knu->U + p->knu->m + 1);
    out.values.insert(out.values.end(), p->knv->U, p->knv->U + p->knv->m + 1);
    for (long i = 0; i <= p->net->n; ++i)
      for (long j = 0; j <= p->net->m; ++j)
        control(out, p->net->Pw[i][j]);
    break;
  }
  case 5: {
    require(c(0)->GetNaturalInterval().ContainsValue(d.at(0)));
    SmVector3d a[3];
    checked(c(0)->Evaluate(d[0], 2, FALSE, a, FALSE));
    for (auto &p : a)
      point(out, p);
    break;
  }
  case 6: {
    SmPoint2d uv(d.at(0), d.at(1));
    require(f(0)->GetNaturalUVDomain().ContainsPoint2d(uv));
    SmVector3d a[6];
    checked(f(0)->Evaluate2ndDerivatives(uv, FALSE, FALSE, a[0], a[1], a[2],
                                         a[4], a[3], a[5]));
    for (auto &p : a)
      point(out, p);
    break;
  }
  case 7: {
    auto r = c(0)->GetNaturalInterval();
    out.values = {r.GetMin(), r.GetMax()};
    break;
  }
  case 8: {
    auto r = f(0)->GetNaturalUVDomain();
    out.values = {r.GetUMin(), r.GetUMax(), r.GetVMin(), r.GetVMax()};
    break;
  }
  case 9:
  case 10: {
    SmSolutionArray a;
    if (op == 9)
      checked(c(0)->GlobalPointSolve(c(0)->GetNaturalInterval(), SM_SO_MINIMIZE,
                                     v(0), s.linear, nullptr, nullptr,
                                     SM_SR_ALL, a));
    else
      checked(f(0)->GlobalPointSolve(f(0)->GetNaturalUVDomain(), SM_SO_MINIMIZE,
                                     v(0), SmZoneTol3d(s.linear), nullptr,
                                     SM_SR_ALL, a));
    solutions(out, a, op == 9 ? 1 : 2);
    break;
  }
  case 11:
  case 12: {
    // Object identity is an exact coincidence decision, not a numerical
    // self-intersection query (which searches for distinct parameter pairs).
    if (op == 11 && ids.at(0) == ids.at(1)) {
      auto r = c(0)->GetNaturalInterval();
      out.indices = {size_t(SM_ST_RANGE_OF_VALUES)};
      out.values = {r.GetMin(), r.GetMin(), r.GetMax(), r.GetMax()};
      break;
    }
    SmSolutionArray a;
    if (op == 11)
      checked(c(0)->GlobalCurveIntersect(c(0)->GetNaturalInterval(), *c(1),
                                         c(1)->GetNaturalInterval(), s.linear,
                                         a));
    else
      checked(f(1)->GlobalCurveIntersect(f(1)->GetNaturalUVDomain(), *c(0),
                                         c(0)->GetNaturalInterval(), s.linear,
                                         a));
    solutions(out, a, op == 11 ? 2 : 3);
    break;
  }
  case 13: {
    if (ids.at(0) == ids.at(1))
      throw Failure{6, 0}; // an overlapping area is not a curve branch
    ArrayOwner<SmCurve> xyz, uv1, uv2;
    SmTArray<SmTsectCurveType> types;
    SmTArray<double> deviations;
    SmBoolean edges[2] = {TRUE, TRUE};
    SmApproxTol3d tol(s.linear);
    checked(f(0)->GlobalSurfaceIntersect(
        s.context, f(0)->GetNaturalUVDomain(), *f(1),
        f(1)->GetNaturalUVDomain(), edges, &tol, &s.angular, &xyz.array,
        &uv1.array, &uv2.array, &types, &deviations));
    size_t count = xyz.array.GetSize();
    require(uv1.array.GetSize() == count && uv2.array.GetSize() == count &&
            types.GetSize() == count && deviations.GetSize() == count);
    xyz.take(s, out);
    uv1.take(s, out);
    uv2.take(s, out);
    for (size_t i = 0; i < count; ++i) {
      out.indices.push_back(size_t(types[i]));
      out.values.push_back(deviations[i]);
    }
    break;
  }
  case 20: {
    SmLine *p = nullptr;
    DeleteOnExit<SmLine> guard(p);
    auto st = SmApiCreateLineSegment(v(0), v(3), p);
    result(s, out, st, p);
    break;
  }
  case 21: {
    auto st = SmApiCreateCircle(v(0), d.at(3), bc);
    result(s, out, st, bc);
    break;
  }
  case 22: {
    SmEllipse *p = nullptr;
    DeleteOnExit<SmEllipse> guard(p);
    auto st = SmApiCreateEllipse(v(0), d.at(3), d.at(4), p);
    result(s, out, st, p);
    break;
  }
  case 23: {
    auto st = SmApiCreateArc(v(0), d.at(3), d.at(4), d.at(5), bc);
    result(s, out, st, bc);
    break;
  }
  case 24: {
    auto st = SmApiCreateHelix(v(0), d.at(3), d.at(4), d.at(5), d.at(6),
                               n.at(0) != 0, s.linear, bc);
    result(s, out, st, bc);
    break;
  }
  case 25: {
    auto p = pts();
    auto st = SmApiCreateCurveInterpPoints(
        p, curve, static_cast<SmCurveParameterizationType>(n.at(0)));
    result(s, out, st, curve);
    break;
  }
  case 26: {
    auto p = pts();
    auto st = SmApiCreateCurveApproxPoints(p, curve);
    result(s, out, st, curve);
    break;
  }
  case 27: {
    auto p = pts();
    auto st = SmApiCreateCurve(p, bc);
    result(s, out, st, bc);
    break;
  }
  case 28:
  case 29: {
    ArrayOwner<SmBSplineCurve> a;
    auto st = op == 28
                  ? SmApiCreateRectangle(v(0), d.at(3), d.at(4), a.array)
                  : SmApiCreateRegularPolygon(v(0), n.at(0), d.at(3), a.array);
    checked(st);
    a.take(s, out);
    break;
  }
  case 30: {
    auto st =
        SmApiCreateSurfaceFromCornerPoints(v(0), v(3), v(6), v(9), surface);
    result(s, out, st, surface);
    break;
  }
  case 31:
  case 32:
  case 33: {
    auto p = pts();
    SmStatus st;
    if (op == 31)
      st = SmApiCreateSurface(p, n.at(0), n.at(1), surface);
    else if (op == 32)
      st = SmApiCreateSurfaceFromOrderedPoints(p, n.at(0), n.at(1), surface);
    else
      st = SmApiCreateSurfaceFromRandomPoints(p, surface);
    result(s, out, st, surface);
    break;
  }
  case 34: {
    auto st = SmApiCreateExtrude(b(0), v(0), surface);
    result(s, out, st, surface);
    break;
  }
  case 35: {
    SmBSplineSurface *p = nullptr;
    DeleteOnExit<SmBSplineSurface> guard(p);
    auto st = SmBSplineSurface::CreateSweepSurface(s.context, b(0), b(1),
                                                   nullptr, FALSE, s.linear, p);
    result(s, out, st, p);
    break;
  }
  case 36: {
    auto st = SmApiCreateRuledSurface(b(0), b(1), surface);
    result(s, out, st, surface);
    break;
  }
  case 37: {
    auto st = SmApiCreateSurfaceRevolution(b(0), v(0), v(3), d.at(6), surface);
    result(s, out, st, surface);
    break;
  }
  case 38: {
    SmTArray<SmBSplineCurve *> a;
    for (size_t i = 0; i < ids.size(); ++i)
      a.Add(b(i));
    auto st = SmApiCreateSkin(a, surface);
    result(s, out, st, surface);
    break;
  }
  case 39: {
    auto st = f(0)->CreateOffsetSurface(s.context, d.at(0),
                                        SmApproxTol3d(s.linear), surface);
    result(s, out, st, surface);
    break;
  }
  case 40: {
    SmVector3d normal;
    if (d.size() == 4)
      normal = v(1);
    else if (!b(0)->IsPlanar(s.linear, nullptr, &normal))
      throw Failure{6, 0};
    double deviation = 0;
    auto st = b(0)->CreateSimpleOffset(s.context, s.linear, normal, d.at(0), bc,
                                       deviation, FALSE);
    result(s, out, st, bc);
    break;
  }
  case 41:
  case 42:
  case 43: {
    ArrayOwner<SmBSplineCurve> uv, xyz;
    double drop = 0, gap = 0;
    SmApproxTol3d tol(s.linear);
    if (op != 41) {
      SmStatus st =
          op == 42
              ? f(1)->DropCurve(s.context, f(1)->GetNaturalUVDomain(), *b(0),
                                b(0)->GetNaturalInterval(), tol, drop, gap,
                                uv.array, TRUE)
              : f(1)->DropAndTrimCurve(s.context, f(1)->GetNaturalUVDomain(),
                                       *b(0), b(0)->GetNaturalInterval(), tol,
                                       drop, gap, uv.array, TRUE, FALSE);
      checked(st);
    }
    const size_t count = op == 41 ? 1 : uv.array.GetSize();
    for (size_t i = 0; i < count; ++i) {
      auto input = op == 41 ? b(0) : uv.array[i];
      SmBSplineCurve *lifted = nullptr;
      DeleteOnExit<SmBSplineCurve> guard(lifted);
      auto st = f(1)->LiftCurve(s.context, f(1)->GetNaturalUVDomain(), *input,
                                input->GetNaturalInterval(), tol, gap, lifted);
      std::unique_ptr<SmBSplineCurve> owner(lifted);
      lifted = nullptr;
      checked(st);
      require(bool(owner));
      xyz.array.Add(owner.get());
      owner.release();
    }
    out.values.push_back(double(uv.array.GetSize()));
    uv.take(s, out);
    xyz.take(s, out);
    break;
  }
  case 44:
  case 45: {
    ArrayOwner<SmBSplineCurve> a;
    for (size_t i = 0; i < ids.size(); ++i)
      a.append(copy_curve(s, *b(i)));
    checked(op == 44 ? SmApiMakeCurvesCompatible(a.array)
                     : SmApiOrderCurves(a.array));
    a.take(s, out);
    break;
  }
  case 46:
  case 49:
  case 51:
  case 53:
  case 54: {
    auto p = copy_curve(s, *b(0));
    SmStatus st;
    if (op == 46)
      st = p->RemoveExtraKnots(s.linear);
    else if (op == 49) {
      SmExtent1d range(d.at(0), d.at(1));
      st = p->Trim(range);
    } else if (op == 51) {
      SmExtent1d r;
      st = p->ReverseParameterization(p->GetNaturalInterval(), r);
    } else if (op == 53)
      st = p->DegreeElevate(n.at(0));
    else {
      SmTArray<double> k;
      for (double x : d)
        k.Add(x);
      st = p->InsertKnots(k);
    }
    checked(st);
    add(s, out, std::move(p));
    break;
  }
  case 47: {
    ArrayOwner<SmBSplineCurve> a;
    checked(b(0)->SubdivideToLinesAndArcs(s.context, s.linear, a.array));
    if (a.array.GetSize() == 0)
      a.append(copy_curve(s, *b(0)));
    a.take(s, out);
    break;
  }
  case 48: {
    SmContinuityType t;
    checked(SmApiEvaluateContinuity(c(0), d.at(0), t));
    out.indices.push_back(size_t(t));
    break;
  }
  case 50:
  case 52: {
    auto p = copy_surface(s, *s.get<SmBSplineSurface>(ids.at(0)));
    if (op == 50) {
      SmExtent2d r(d.at(0), d.at(2), d.at(1), d.at(3));
      checked(p->TrimWithDomain(r));
    } else
      checked(p->Reverse(static_cast<SmSurfParamType>(n.at(0))));
    add(s, out, std::move(p));
    break;
  }
  case 55: {
    auto st = f(0)->CreateIsoParametricCurve(
        s.context, static_cast<SmSurfParamType>(n.at(0)), d.at(0),
        SmApproxTol3d(s.linear), bc);
    result(s, out, st, bc);
    break;
  }
  case 56:
  case 57: {
    auto p = copy_surface(s, *s.get<SmBSplineSurface>(ids.at(0)));
    auto direction = static_cast<SmSurfParamType>(n.at(0));
    if (op == 56)
      checked(p->DegreeElevate(direction, n.at(1)));
    else {
      SmTArray<double> k;
      for (double x : d)
        k.Add(x);
      checked(p->InsertKnots(direction, k));
    }
    add(s, out, std::move(p));
    break;
  }
  case 58: {
    SmSurface *left = nullptr, *right = nullptr;
    DeleteOnExit<SmSurface> lg(left), rg(right);
    auto st = s.get<SmBSplineSurface>(ids.at(0))->SplitAt(
        s.context, d.at(0), static_cast<SmSurfParamType>(n.at(0)), left, right);
    std::unique_ptr<SmObject> l(left), r(right);
    left = nullptr;
    right = nullptr;
    checked(st);
    add(s, out, std::move(l));
    add(s, out, std::move(r));
    break;
  }
  case 60: {
    require(n.at(0) > 0 && n[0] <= 8);
    SmBSplineSurface *approximation = nullptr;
    DeleteOnExit<SmBSplineSurface> guard(approximation);
    double achieved = 0;
    auto status = f(0)->ApproximateSurface(s.context, SmApproxTol3d(s.linear),
                                          achieved, static_cast<ULONG>(n[0]),
                                          approximation, TRUE);
    result(s, out, status, approximation);
    out.values.push_back(achieved);
    break;
  }
  default:
    throw Failure{7, 0};
  }
}
} // namespace
extern "C" int gs_sm_create(double l, double a, double size, void **result) {
  if (!result)
    return 1;
  *result = nullptr;
  try {
    std::lock_guard<std::mutex> lock(gate);
    require(std::isfinite(l) && std::isfinite(a) && std::isfinite(size) &&
            l > 0 && a > 0 && size > l);
    SmApiCreateContext();
    *result = new Session(l, a, size);
    return 0;
  } catch (const Failure &e) {
    return e.code;
  } catch (const std::bad_alloc &) {
    return 8;
  } catch (...) {
    return 10;
  }
}
extern "C" void gs_sm_destroy(void *p) {
  try {
    std::lock_guard<std::mutex> lock(gate);
    delete static_cast<Session *>(p);
  } catch (...) {
  }
}
extern "C" void gs_sm_result_free(void *p) {
  delete static_cast<Response *>(p);
}
extern "C" int gs_sm_call(void *ptr, uint32_t op, const size_t *ids, size_t ni,
                          const double *values, size_t nv, const size_t *sizes,
                          size_t ns, GsSmResult *output) {
  if (!output)
    return 1;
  *output = {};
  try {
    std::lock_guard<std::mutex> lock(gate);
    require(ptr && (!ni || ids) && (!nv || values) && (!ns || sizes));
    std::vector<size_t> ii, nn;
    std::vector<double> dd;
    if (ni)
      ii.assign(ids, ids + ni);
    if (ns)
      nn.assign(sizes, sizes + ns);
    if (nv)
      dd.assign(values, values + nv);
    for (double x : dd)
      require(std::isfinite(x));
    auto r = std::make_unique<Response>();
    auto &session = *static_cast<Session *>(ptr);
    const auto before = session.objects.size();
    try {
      dispatch(session, op, ii, dd, nn, *r);
    } catch (...) {
      session.objects.resize(before);
      throw;
    }
    output->values = r->values.data();
    output->values_len = r->values.size();
    output->indices = r->indices.data();
    output->indices_len = r->indices.size();
    output->owner = r.release();
    return 0;
  } catch (const Failure &e) {
    output->upstream_status = e.native;
    return e.code;
  } catch (const std::out_of_range &) {
    return 1;
  } catch (const std::bad_alloc &) {
    return 8;
  } catch (...) {
    return 10;
  }
}
