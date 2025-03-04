// Original header license: SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Whitney Armstrong, Sylvester Joosten, Wouter Deconinck, Dmitry Romanov
//

#pragma once


#include <spdlog/logger.h>
#include "JugTrack/TrackingResultTrajectory.hpp"
#include <algorithms/tracking/ParticlesFromTrackFitResult.h>


namespace eicrecon::Reco {

    /** Extract the particles form fit trajectories.
     *
     * \ingroup tracking
     */
    class ParticlesFromTrackFit {
    private:
        std::shared_ptr<spdlog::logger> m_log;


    public:
        void init(std::shared_ptr<spdlog::logger> log);

        ParticlesFromTrackFitResult *execute(const std::vector<const eicrecon::TrackingResultTrajectory *> &trajectories);

    };
} // namespace Jug::Reco
