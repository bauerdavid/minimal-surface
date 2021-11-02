#pragma once
#include <utility>
#include "AreaEikonal.h"

std::pair<IPoi3<double>, IPoi3<double>> best_plane_from_points(const std::unordered_set<IPoi3<double>, IPoi3Hash<double>>& c);