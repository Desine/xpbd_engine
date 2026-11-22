#include "cluster.hpp"
#include "entity.hpp"
#include "particle.hpp"
#include "constraints.hpp"
#include "colliders.hpp"
#include "actuators.hpp"

#include <vector>
#include <string>

namespace xpbd
{
    void merge_clusters(Cluster &into, Cluster &from)
    {
        if (&into == &from)
            return;

        const size_t offset = into.particles.size();

        // Particles
        into.particles.insert(into.particles.end(),
                                  from.particles.begin(), from.particles.end());

        // Distance constraints
        into.distanceConstraints.reserve(into.distanceConstraints.size() + from.distanceConstraints.size());
        for (auto &dc : from.distanceConstraints)
        {
            dc.i1 += offset;
            dc.i2 += offset;
            into.distanceConstraints.push_back(dc);
        }

        // Angle constraints
        into.angleConstraints.reserve(into.angleConstraints.size() + from.angleConstraints.size());
        for (auto &ac : from.angleConstraints)
        {
            ac.i1 += offset;
            ac.i2 += offset;
            ac.i3 += offset;
            ac.i4 += offset;
            into.angleConstraints.push_back(ac);
        }

        // Volume constraints
        into.volumeConstraints.reserve(into.volumeConstraints.size() + from.volumeConstraints.size());
        for (auto &vc : from.volumeConstraints)
        {
            for (auto &idx : vc.indices)
                idx += offset;
            into.volumeConstraints.push_back(vc);
        }

        // Polygon colliders
        into.polygonColliders.reserve(into.polygonColliders.size() + from.polygonColliders.size());
        for (auto &pc : from.polygonColliders)
        {
            for (auto &idx : pc.indices)
                idx += offset;
            into.polygonColliders.push_back(pc);
        }

        // Points colliders
        into.pointsColliders.reserve(into.pointsColliders.size() + from.pointsColliders.size());
        for (auto &pc : from.pointsColliders)
        {
            for (auto &idx : pc.indices)
                idx += offset;
            into.pointsColliders.push_back(pc);
        }

        // Angular motors
        into.angularMotors.reserve(into.angularMotors.size() + from.angularMotors.size());
        for (auto &am : from.angularMotors)
        {
            am.i1 += offset;
            am.i2 += offset;
            am.i3 += offset;
            am.i4 += offset;
            into.angularMotors.push_back(am);
        }

        // Angular servos
        into.angularServos.reserve(into.angularServos.size() + from.angularServos.size());
        for (auto &as : from.angularServos)
        {
            as.i1 += offset;
            as.i2 += offset;
            as.i3 += offset;
            as.i4 += offset;
            into.angularServos.push_back(as);
        }

        from = Cluster();
    }
} // namespace xpbd