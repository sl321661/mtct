#include <cstdlib>
#define TEST_FRIENDS true

#include "CustomExceptions.hpp"
#include "EOMHelper.hpp"
#include "datastructure/Timetable.hpp"
#include "probleminstances/GeneralPerformanceOptimizationInstance.hpp"
#include "simulator/EventSimulator.hpp"

#include "gtest/gtest.h"
#include <cmath>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Init.h>
#include <plog/Log.h>
#include <plog/Severity.h>

using namespace cda_rail;

#define EXPECT_APPROX_EQ_6(a, b) EXPECT_APPROX_EQ(a, b, 1e-6)

#define EXPECT_APPROX_EQ(a, b, c)                                              \
  EXPECT_TRUE(std::abs((a) - (b)) < (c)) << (a) << " !=(approx.) " << (b)

// NOLINTBEGIN
// (clang-analyzer-deadcode.DeadStores,misc-const-correctness,clang-diagnostic-unused-result)

TEST(EventSimulator, CalculateOptimalTrainMovements) {
  instances::GeneralPerformanceOptimizationInstance instance{};

  simulator::EventSimulator simulator{instance, {}};

  // TODO: limit_speed_by_leaving_edges true oder false?
  std::vector<cda_rail::simulator::EventSimulator::EdgeSegment> free_track = {
      {1000, 50}, {2000, 100}, {200, 90}, {1000, 50}};
  auto initial_speed = 10.0;
  auto train         = Train("Deutsche Bahn", 300, 100, 10, 10);

  auto movements = simulator.calculate_optimal_train_movements(
      free_track, initial_speed, train);

  // EXPECT_EQ(movements[0])
}

// NOLINTEND
// (clang-analyzer-deadcode.DeadStores,misc-const-correctness,clang-diagnostic-unused-result)
