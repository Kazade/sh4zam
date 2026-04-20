#include "shz_test.h"
#include "shz_test.hpp"
#include "sh4zam/shz_quat.hpp"

#include <print>

#define GBL_SELF_TYPE   shz_quat_test_suite

GBL_TEST_FIXTURE_NONE
GBL_TEST_INIT_NONE
GBL_TEST_FINAL_NONE

namespace {
    shz_vec3_t raylibVector3RotateByQuaternion(shz_quat_t q, shz_vec3_t v) {
        shz_vec3_t result = { 0 };

        result.x = v.x*(q.x*q.x + q.w*q.w - q.y*q.y - q.z*q.z) + v.y*(2*q.x*q.y - 2*q.w*q.z) + v.z*(2*q.x*q.z + 2*q.w*q.y);
        result.y = v.x*(2*q.w*q.z + 2*q.x*q.y) + v.y*(q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z) + v.z*(-2*q.w*q.x + 2*q.y*q.z);
        result.z = v.x*(-2*q.w*q.y + 2*q.x*q.z) + v.y*(2*q.w*q.x + 2*q.y*q.z)+ v.z*(q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z);

        return result;
    }
}

GBL_TEST_CASE(identity)
    const auto q = shz::quat::identity();
    GBL_TEST_VERIFY(q.w == 1.0f);
    GBL_TEST_VERIFY(q.x == 0.0f);
    GBL_TEST_VERIFY(q.y == 0.0f);
    GBL_TEST_VERIFY(q.z == 0.0f);
    GBL_TEST_VERIFY(q == shz::quat(1.0f, 0.0f, 0.0f, 0.0f));
    GBL_TEST_VERIFY(shz_equalf(q.magnitude(), 1.0f));
GBL_TEST_CASE_END

GBL_TEST_CASE(from_axis_angle)
    // zero angle around any axis yields identity
    GBL_TEST_VERIFY(shz::quat::from_axis_angle({1.0f, 0.0f, 0.0f}, 0.0f) == shz::quat::identity());
    GBL_TEST_VERIFY(shz::quat::from_axis_angle({0.0f, 1.0f, 0.0f}, 0.0f) == shz::quat::identity());
    GBL_TEST_VERIFY(shz::quat::from_axis_angle({0.0f, 0.0f, 1.0f}, 0.0f) == shz::quat::identity());

    // 180° around each basis axis: w≈0, corresponding component≈1, others=0
    GBL_TEST_VERIFY(shz::quat::from_axis_angle({1.0f, 0.0f, 0.0f}, SHZ_F_PI) == shz::quat(0.0f, 1.0f, 0.0f, 0.0f));
    GBL_TEST_VERIFY(shz::quat::from_axis_angle({0.0f, 1.0f, 0.0f}, SHZ_F_PI) == shz::quat(0.0f, 0.0f, 1.0f, 0.0f));
    GBL_TEST_VERIFY(shz::quat::from_axis_angle({0.0f, 0.0f, 1.0f}, SHZ_F_PI) == shz::quat(0.0f, 0.0f, 0.0f, 1.0f));

    // 90° around each basis axis: verify axis()/angle() round-trip
    {
        const auto q = shz::quat::from_axis_angle({1.0f, 0.0f, 0.0f}, SHZ_F_PI_2);
        GBL_TEST_VERIFY(q.axis() == shz::vec3(1.0f, 0.0f, 0.0f));
        GBL_TEST_VERIFY(shz_equalf(q.angle(), SHZ_F_PI_2));
        GBL_TEST_VERIFY(shz_equalf(q.magnitude(), 1.0f));
    }
    {
        const auto q = shz::quat::from_axis_angle({0.0f, 1.0f, 0.0f}, SHZ_F_PI_2);
        GBL_TEST_VERIFY(q.axis() == shz::vec3(0.0f, 1.0f, 0.0f));
        GBL_TEST_VERIFY(shz_equalf(q.angle(), SHZ_F_PI_2));
        GBL_TEST_VERIFY(shz_equalf(q.magnitude(), 1.0f));
    }
    {
        const auto q = shz::quat::from_axis_angle({0.0f, 0.0f, 1.0f}, SHZ_F_PI_2);
        GBL_TEST_VERIFY(q.axis() == shz::vec3(0.0f, 0.0f, 1.0f));
        GBL_TEST_VERIFY(shz_equalf(q.angle(), SHZ_F_PI_2));
        GBL_TEST_VERIFY(shz_equalf(q.magnitude(), 1.0f));
    }

    // arbitrary axis: axis()/angle() round-trip
    {
        const shz::vec3 axis  = shz::vec3(1.0f, 1.0f, 0.0f).direction();
        const float     angle = SHZ_F_PI / 3.0f;
        const auto q = shz::quat::from_axis_angle(axis, angle);
        GBL_TEST_VERIFY(q.axis() == axis);
        GBL_TEST_VERIFY(shz_equalf(q.angle(), angle));
        GBL_TEST_VERIFY(shz_equalf(q.magnitude(), 1.0f));
    }
GBL_TEST_CASE_END

GBL_TEST_CASE(from_angles_xyz)
    // all-zero angles yield identity
    GBL_TEST_VERIFY(shz::quat::from_angles_xyz(0.0f, 0.0f, 0.0f) == shz::quat::identity());

    // single-axis 180° must agree with from_axis_angle
    GBL_TEST_VERIFY(shz::quat::from_angles_xyz(SHZ_F_PI,   0.0f,       0.0f      ) ==
                    shz::quat::from_axis_angle({1.0f, 0.0f, 0.0f}, SHZ_F_PI));
    GBL_TEST_VERIFY(shz::quat::from_angles_xyz(0.0f,       SHZ_F_PI,   0.0f      ) ==
                    shz::quat::from_axis_angle({0.0f, 1.0f, 0.0f}, SHZ_F_PI));
    GBL_TEST_VERIFY(shz::quat::from_angles_xyz(0.0f,       0.0f,       SHZ_F_PI  ) ==
                    shz::quat::from_axis_angle({0.0f, 0.0f, 1.0f}, SHZ_F_PI));

    // single-axis 90° must agree with from_axis_angle
    GBL_TEST_VERIFY(shz::quat::from_angles_xyz(SHZ_F_PI_2, 0.0f,       0.0f      ) ==
                    shz::quat::from_axis_angle({1.0f, 0.0f, 0.0f}, SHZ_F_PI_2));
    GBL_TEST_VERIFY(shz::quat::from_angles_xyz(0.0f,       SHZ_F_PI_2, 0.0f      ) ==
                    shz::quat::from_axis_angle({0.0f, 1.0f, 0.0f}, SHZ_F_PI_2));
    GBL_TEST_VERIFY(shz::quat::from_angles_xyz(0.0f,       0.0f,       SHZ_F_PI_2) ==
                    shz::quat::from_axis_angle({0.0f, 0.0f, 1.0f}, SHZ_F_PI_2));

    // produced quaternions must be unit quaternions
    GBL_TEST_VERIFY(shz_equalf(shz::quat::from_angles_xyz(0.5f,       1.0f,       1.5f      ).magnitude(), 1.0f));
    GBL_TEST_VERIFY(shz_equalf(shz::quat::from_angles_xyz(SHZ_F_PI_4, SHZ_F_PI_2, SHZ_F_PI  ).magnitude(), 1.0f));
GBL_TEST_CASE_END

GBL_TEST_CASE(transform_vec3)
    auto check = [](shz::quat q, volatile shz::vec3 v, shz::vec3 expected) {
        return q.transform(v) == expected;
    };

    auto cross_check = [](shz::quat q, volatile shz::vec3 v) {
        return q.transform(v) == shz::vec3(raylibVector3RotateByQuaternion(q, const_cast<shz::vec3&>(v)));
    };

    // identity leaves any vector unchanged
    GBL_TEST_VERIFY(check(shz::quat::identity(), { 1.0f,  2.0f,  3.0f}, { 1.0f,  2.0f,  3.0f}));
    GBL_TEST_VERIFY(check(shz::quat::identity(), {-4.0f,  0.0f,  7.5f}, {-4.0f,  0.0f,  7.5f}));
    GBL_TEST_VERIFY(check(shz::quat::identity(), { 0.0f,  0.0f,  0.0f}, { 0.0f,  0.0f,  0.0f}));

    {
        const auto rot90x  = shz::quat::from_axis_angle({1.0f, 0.0f, 0.0f}, SHZ_F_PI_2);
        const auto rot90y  = shz::quat::from_axis_angle({0.0f, 1.0f, 0.0f}, SHZ_F_PI_2);
        const auto rot90z  = shz::quat::from_axis_angle({0.0f, 0.0f, 1.0f}, SHZ_F_PI_2);
        const auto rot180x = shz::quat::from_axis_angle({1.0f, 0.0f, 0.0f}, SHZ_F_PI);
        const auto rot180y = shz::quat::from_axis_angle({0.0f, 1.0f, 0.0f}, SHZ_F_PI);
        const auto rot180z = shz::quat::from_axis_angle({0.0f, 0.0f, 1.0f}, SHZ_F_PI);

        // 90° around X: X unchanged, Y→Z, Z→-Y
        GBL_TEST_VERIFY(check(rot90x, {1.0f, 0.0f, 0.0f}, { 1.0f,  0.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot90x, {0.0f, 1.0f, 0.0f}, { 0.0f,  0.0f,  1.0f}));
        GBL_TEST_VERIFY(check(rot90x, {0.0f, 0.0f, 1.0f}, { 0.0f, -1.0f,  0.0f}));

        // 90° around Y: Y unchanged, Z→X, X→-Z
        GBL_TEST_VERIFY(check(rot90y, {0.0f, 1.0f, 0.0f}, { 0.0f,  1.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot90y, {0.0f, 0.0f, 1.0f}, { 1.0f,  0.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot90y, {1.0f, 0.0f, 0.0f}, { 0.0f,  0.0f, -1.0f}));

        // 90° around Z: Z unchanged, X→Y, Y→-X
        GBL_TEST_VERIFY(check(rot90z, {0.0f, 0.0f, 1.0f}, { 0.0f,  0.0f,  1.0f}));
        GBL_TEST_VERIFY(check(rot90z, {1.0f, 0.0f, 0.0f}, { 0.0f,  1.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot90z, {0.0f, 1.0f, 0.0f}, {-1.0f,  0.0f,  0.0f}));

        // 180° around X: X unchanged, Y→-Y, Z→-Z
        GBL_TEST_VERIFY(check(rot180x, {1.0f, 0.0f, 0.0f}, { 1.0f,  0.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot180x, {0.0f, 1.0f, 0.0f}, { 0.0f, -1.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot180x, {0.0f, 0.0f, 1.0f}, { 0.0f,  0.0f, -1.0f}));

        // 180° around Y: Y unchanged, X→-X, Z→-Z
        GBL_TEST_VERIFY(check(rot180y, {0.0f, 1.0f, 0.0f}, { 0.0f,  1.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot180y, {1.0f, 0.0f, 0.0f}, {-1.0f,  0.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot180y, {0.0f, 0.0f, 1.0f}, { 0.0f,  0.0f, -1.0f}));

        // 180° around Z: Z unchanged, X→-X, Y→-Y
        GBL_TEST_VERIFY(check(rot180z, {0.0f, 0.0f, 1.0f}, { 0.0f,  0.0f,  1.0f}));
        GBL_TEST_VERIFY(check(rot180z, {1.0f, 0.0f, 0.0f}, {-1.0f,  0.0f,  0.0f}));
        GBL_TEST_VERIFY(check(rot180z, {0.0f, 1.0f, 0.0f}, { 0.0f, -1.0f,  0.0f}));

        // zero vector is invariant under any rotation
        GBL_TEST_VERIFY(check(rot90x,  {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}));
        GBL_TEST_VERIFY(check(rot180y, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}));
    }

    // cross-validate arbitrary rotations against reference implementation
    GBL_TEST_VERIFY(cross_check(shz::quat::from_axis_angle(shz::vec3{1.0f, 1.0f, 0.0f}.direction(), SHZ_F_PI / 3.0f), {1.0f, 2.0f, 3.0f}));
    GBL_TEST_VERIFY(cross_check(shz::quat::from_axis_angle(shz::vec3{1.0f, 1.0f, 1.0f}.direction(), SHZ_F_PI / 4.0f), {-1.0f, 0.5f, 2.0f}));
    GBL_TEST_VERIFY(cross_check(shz::quat::from_angles_xyz(0.5f, 1.0f, 1.5f), {3.0f, -1.0f, 0.0f}));
    GBL_TEST_VERIFY(cross_check(shz::quat::from_angles_xyz(0.3f, 0.7f, 0.2f), {0.0f, 1.0f, -1.0f}));

    GBL_TEST_VERIFY(
        (benchmark_cmp<shz::vec3>)("shz::quat::transform(shz::vec3)",
                                   [](const shz::quat& q, const shz::vec3& v) {
                                        return q.transform(v);
                                   },
                                   "Vector3RotateByQuaternion",
                                   raylibVector3RotateByQuaternion,
                                   shz::quat::from_angles_xyz(0.5f, 1.0f, 1.5f), shz::vec3{ 43.0f, -232.0f, 344.0f })
    );
GBL_TEST_CASE_END

GBL_TEST_REGISTER(identity,
                  from_axis_angle,
                  from_angles_xyz,
                  transform_vec3)