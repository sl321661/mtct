#pragma once

#include "CustomExceptions.hpp"
#include "Definitions.hpp"
#include "datastructure/Train.hpp"
#include "probleminstances/GeneralPerformanceOptimizationInstance.hpp"
#include "simulator/GeneralSimulator.hpp"

// NOLINTNEXTLINE(misc-include-cleaner)
#include "EOMHelper.hpp"

#include "gtest/gtest_prod.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <queue>
#include <span>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef TEST_FRIENDS
#define TEST_FRIENDS false
#endif
#if TEST_FRIENDS
class EventSimulator_CalculateOptimalTrainMovements_Test;
#endif

namespace cda_rail::simulator {

class EventSimulator : public GeneralSimulator {
public:
  // reuse default values
  using GeneralSimulator::simulate;

  // ----------------
  // CONSTRUCTOR
  // ----------------
  explicit EventSimulator(
      cda_rail::instances::GeneralPerformanceOptimizationInstance const&
                                       instance,
      std::vector<cda_rail::index_set> ttd_sections);
  explicit EventSimulator(
      cda_rail::instances::GeneralPerformanceOptimizationInstance const&
                                          instance,
      std::vector<cda_rail::index_set>    ttd_sections,
      std::vector<cda_rail::index_vector> train_edges,
      std::vector<cda_rail::index_vector> ttd_orders,
      std::vector<cda_rail::index_vector> vertex_orders,
      std::vector<std::vector<double>>    stop_positions);

  // Rule of 5
  EventSimulator(EventSimulator const&)            = default;
  EventSimulator(EventSimulator&&)                 = default;
  EventSimulator& operator=(EventSimulator const&) = default;
  EventSimulator& operator=(EventSimulator&&)      = default;
  ~EventSimulator() override                       = default;

  /**
   * @brief Retrieves the typed performance optimization instance.
   *
   * @return Const pointer to the `GeneralPerformanceOptimizationInstance` for
   * this simulator.
   */
  [[nodiscard]] instances::GeneralPerformanceOptimizationInstance const*
  get_instance() const override {
    return dynamic_cast<
        instances::GeneralPerformanceOptimizationInstance const*>(
        GeneralSimulator::get_instance());
  }

  // ------------------
  // SIMULATION
  // ------------------

  /**
   * @brief Runs the simulation.
   *
   * @return SimulatorResults containing the simulation outcome.
   */
  [[nodiscard]] SimulatorResults
  simulate(bool late_entry_possible, bool limit_speed_by_leaving_edges,
           bool save_trajectories,
           bool disappear_at_partial_route_end) const override;

private:
#if TEST_FRIENDS
  FRIEND_TEST(::EventSimulator, CalculateOptimalTrainMovements);
#endif

  /**
   * Describes linear train movement with initial velocity u, final velocity v,
   * acceleration a, distance traveled s.
   */
  struct TrainMovement {
    TrainMovement(double u, double v, double a, double s) {
      this->u = u;
      this->v = v;
      this->a = a;
      this->s = s;
    }

    TrainMovement(double u, double v, double a)
        : TrainMovement(u, v, a, distance_travelled(u, v, a)) {}

    double u;
    double v;
    double a;
    double s;
  };

  /**
   * Describes (part of) an edge, storing only length and max_speed
   */
  struct EdgeSegment {
    double length;
    double max_speed;
  };

  /**
   * Holds a train's front and read positions
   */
  struct TrainPosition {
    double rear;
    double front;
  };

  struct TrainTimestep {
    size_t tr;
    double timestep;
  };

  struct MaxSpeedChangePoint {
    double position;
    double new_max_speed;
  };

  struct Compare {
    bool operator()(const TrainTimestep& a, const TrainTimestep& b) const {
      return a.timestep > b.timestep;
    }
  };

  using train_timestep_queue =
      std::priority_queue<TrainTimestep, std::vector<TrainTimestep>, Compare>;

  /**
   * Get the position in the current_train's path where the train no longer has
   * a free path, either due to another train on a shared edge or a train which
   * is yet to cross a vertex first.
   * @param current_train ID of current train to find the free track end for
   * @param train_positions Train positions
   * @param trains_in_network Trains in network, only these will be considered
   * for collisions on edges
   * @param trains_on_edges
   */
  void find_free_track_end(
      size_t current_train,
      const std::vector<::cda_rail::simulator::EventSimulator::TrainPosition>&
                                                     train_positions,
      const std::unordered_set<size_t>&              trains_in_network,
      const std::vector<std::unordered_set<size_t>>& trains_on_edges) const;

  /**
   *
   * @param free_track Vector of EdgeSegments describing the free track ahead
   * @param initial_speed initial speed at which train is traveling
   * @param train The train to calculate movements for
   * @return Series of train movements to traverse the free track as fast as
   * possible
   */
  std::vector<TrainMovement>
  calculate_optimal_train_movements(const std::vector<EdgeSegment>& free_track,
                                    double       initial_speed,
                                    const Train& train) const;

  /**
   * Calculates train movement which traverses a span of tracks as quickly as
   * possible and appends these to the movements vector passed by reference.
   * @param movements Vector of train movements which the newly generated
   * movements will be appended to
   * @param free_track Span of free tracks which should be traversed as fast as
   * possible
   * @param initial_speed Initial train speed
   * @param train Train object which is traversing the track
   */
  static void
  append_max_permitted_speed_movements(std::vector<TrainMovement>&  movements,
                                       std::span<const EdgeSegment> free_track,
                                       double       initial_speed,
                                       const Train& train);

  static double squared(double x);
};
} // namespace cda_rail::simulator
