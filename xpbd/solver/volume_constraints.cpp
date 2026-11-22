#include "volume_constraints.hpp"
#include "xpbd/core/particle.hpp"
#include "xpbd/core/constraints.hpp"
#include "../utils.hpp"

#include <vector>


namespace xpbd
{
    float compute_polygon_area(const std::vector<glm::vec2> &positions)
    {
        float area = 0.0f;
        size_t n = positions.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &p0 = positions[j];
            const glm::vec2 &p1 = positions[i];
            area += Utils::cross(p0, p1);
        }
        return 0.5f * area;
    }
    float compute_polygon_area(const std::vector<Particle> &p, const std::vector<size_t> &indices)
    {
        float area = 0.0f;
        size_t n = indices.size();
        for (size_t i = 0, j = n - 1; i < n; j = i++)
        {
            const glm::vec2 &p0 = p[indices[j]].pos;
            const glm::vec2 &p1 = p[indices[i]].pos;
            area += Utils::cross(p0, p1);
        }
        return 0.5f * area;
    }

    void solve_volume_constraints(std::vector<Particle> &p, std::vector<VolumeConstraint> &vc, float dt)
    {
        for (auto &c : vc)
            solve_volume_constraint(p, c, dt);
    }
    void solve_volume_constraint(std::vector<Particle> &p, VolumeConstraint &vc, float dt)
    {
        std::vector<size_t> &indices = vc.indices;
        size_t N = indices.size();

        float volume = compute_polygon_area(p, indices);
        float C = volume - vc.restVolume;

        float denom = 0.0f;
        std::vector<glm::vec2> grads(N);

        for (size_t i = 0; i < N; ++i)
        {
            const glm::vec2 &pi_prev = p[indices[(i + N - 1) % N]].pos;
            const glm::vec2 &pi_next = p[indices[(i + 1) % N]].pos;

            glm::vec2 grad = 0.5f * Utils::perpCW2d(pi_next - pi_prev);
            grads[i] = grad;

            float wi = p[indices[i]].w;
            denom += wi * glm::dot(grad, grad);
        }

        float alphaTilde = vc.compliance / (dt * dt);
        denom += alphaTilde;
        if (denom < 1e-6f)
            return;

        float deltaLambda = (-C - alphaTilde * vc.lambda) / denom;
        vc.lambda += deltaLambda;

        for (size_t i = 0; i < N; ++i)
        {
            uint32_t idx = indices[i];
            float wi = p[idx].w;

            p[idx].pos += wi * deltaLambda * grads[i];
        }
    } 
    void reset_volume_constraints_lambdas(std::vector<VolumeConstraint> &vc)
    {
        for (auto &c : vc)
            c.lambda = 0;
    }
} // namespace xpbd