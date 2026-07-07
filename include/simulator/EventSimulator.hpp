#pragma once

#include "CustomExceptions.hpp"
#include "Definitions.hpp"
#include "datastructure/Train.hpp"
#include "probleminstances/GeneralPerformanceOptimizationInstance.hpp"
#include "simulator/GeneralSimulator.hpp"

// NOLINTNEXTLINE(misc-include-cleaner)
#include "../../extern/cli11/include/CLI/App.hpp"
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
class EventSimulator_MoveTrain_Test;
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
  FRIEND_TEST(::EventSimulator, MoveTrain);
#endif

  /**
   * Describes linear train movement with initial speed u, final speed v,
   * acceleration a, distance traveled s over time t.
   */
  struct TrainMovement {
    TrainMovement(double u, double v, double a, double s, double t) {
      this->u = u;
      this->v = v;
      this->a = a;
      this->s = s;
      this->t = t;
    }

    TrainMovement(double u, double v, double a, double s)
        : TrainMovement(u, v, a, s, time_taken(u, v, s)) {}

    TrainMovement(double u, double v, double a)
        : TrainMovement(u, v, a, distance_travelled(u, v, a)) {}

    double u;
    double v;
    double a;
    double s;
    double t;
  };

  /**
   * Describes the stop at the end of some free track. The end can either be
   * guaranteed (station, end of line) or temporary (moving train, vertex
   * order/TTD restriction)
   */
  enum class FreeTrackStopType { GuaranteedEnd, TemporaryEnd };

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
   * @param train_positions Train positions
   * @param trains_in_network Trains in network, only these will be considered
   * for collisions on edges
   * @param trains_on_edges
   * @param current_train ID of current train to find the free track end for
   */
  [[nodiscard]] std::pair<std::vector<EdgeSegment>, FreeTrackStopType>
  find_free_track_ahead(
      const std::vector<TrainPosition>&              train_positions,
      const std::unordered_set<size_t>&              trains_in_network,
      const std::vector<std::unordered_set<size_t>>& trains_on_edges,
      size_t                                         current_train) const;

  /**
   *
   * @param free_track Vector of EdgeSegments describing the free track ahead
   * @param initial_speed initial speed at which train is traveling
   * @param train The train to calculate movements for
   * @return Series of train movements to traverse the free track as fast as
   * possible
   */
  [[nodiscard]] static std::vector<TrainMovement>
  calculate_optimal_train_movements(const std::vector<EdgeSegment>& free_track,
                                    double initial_speed, const Train& train);

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

  /**
   *
   * @param movements Train movements to sum time for
   * @return The total time needed to complete all the train movements
   */
  static double
  get_movement_total_time(const std::span<const TrainMovement>& movements);

  static void move_train(TrainPosition& train_position, double& train_velocity,
                         std::vector<TrainMovement>& train_movements, double t);

  /**
   * Calculates the distance traveled in a certain time and changes the
   * train_movement to represent the remaining movement after the time only.
   * @param movement Train movement to partially apply, will be updated
   * in-place to only hold the 'non-consumed' part.
   * @param t The time to apply to the movement
   * @return Distance traveled in time t.
   */
  static double distance_traveled_consume_movement(TrainMovement& movement,
                                                   double         t);

  /**
   * Checks whether a train's exit vertex (from the schedule) is on the train's
   * current path.
   * @param train_id Train ID
   * @return true if the exit vertex appears at the end of the train's path.
   */
  bool path_includes_exit_vertex(const size_t train_id) const;

  static double squared(double x);
};
} // namespace cda_rail::simulator
