#include "simulator/EventSimulator.hpp"

#include "CustomExceptions.hpp"
#include "Definitions.hpp"
#include "EOMHelper.hpp"
#include "GeneralHelper.hpp"
#include "plog/Log.h"
#include "probleminstances/GeneralPerformanceOptimizationInstance.hpp"
#include "simulator/GeneralSimulator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <future>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <span>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// ----------------
// CONSTRUCTOR

cda_rail::simulator::EventSimulator::EventSimulator(
    cda_rail::instances::GeneralPerformanceOptimizationInstance const& instance,
    std::vector<cda_rail::index_set> ttd_sections)
    : GeneralSimulator(
          std::make_shared<
              const instances::GeneralPerformanceOptimizationInstance>(
              instance),
          std::move(ttd_sections)) {}

cda_rail::simulator::EventSimulator::EventSimulator(
    cda_rail::instances::GeneralPerformanceOptimizationInstance const& instance,
    std::vector<cda_rail::index_set>    ttd_sections,
    std::vector<cda_rail::index_vector> train_edges,
    std::vector<cda_rail::index_vector> ttd_orders,
    std::vector<cda_rail::index_vector> vertex_orders,
    std::vector<std::vector<double>>    stop_positions)
    : GeneralSimulator(
          std::make_shared<
              const instances::GeneralPerformanceOptimizationInstance>(
              instance),
          std::move(ttd_sections), std::move(train_edges),
          std::move(ttd_orders), std::move(vertex_orders),
          std::move(stop_positions)) {}

// ---------------
// SIMULATE

cda_rail::simulator::SimulatorResults
cda_rail::simulator::EventSimulator::simulate(
    bool late_entry_possible, bool limit_speed_by_leaving_edges,
    bool save_trajectories, bool disappear_at_partial_route_end) const {
  cda_rail::initialize_plog(true);

  // Initialize return values
  auto const&         train_list       = get_instance()->get_const_train_list();
  auto const          number_of_trains = train_list.size();
  std::vector<double> exit_times(
      number_of_trains,
      0); // Initially each train is simulated until
          // time t=0, for t>0 the heuristic is needed.
  std::vector<std::vector<double>> stop_times(number_of_trains);
  std::vector<std::map<double, PosVel>>
      train_trajectories; // time -> {pos, vel}

  // Initialize variables to keep track of positions and velocities
  if (save_trajectories) {
    train_trajectories.clear();
    train_trajectories.resize(number_of_trains);
  }

  std::vector<TrainPosition> train_positions(number_of_trains,
                                             {.rear = -1.0, .front = -1.0});
  std::vector<double> train_velocities(number_of_trains, -1.0); // velocities
  std::vector<std::vector<TrainMovement>> train_movements(
      number_of_trains, std::vector<TrainMovement>());
  std::unordered_set<size_t> trains_in_network;
  std::unordered_set<size_t> trains_left;
  std::unordered_set<size_t> trains_finished_simulating;
  std::vector<double>        vertex_headways(
      get_instance()->get_const_network().number_of_vertices(),
      0.0); // headways for vertices
  train_timestep_queue next_train_timestep;

  // Get first timestep for each train and detect trains that are not scheduled
  //  to enter the network
  for (size_t tr = 0; tr < number_of_trains; ++tr) {
    if (get_train_edges_of_tr(tr).empty()) {
      trains_finished_simulating.insert(tr);
      continue;
    }

    auto const entry_time =
        get_instance()->get_const_schedule(tr).get_entry_time();
    next_train_timestep.push({tr, entry_time});
  }

  // TODO: Vertex headways in network

  auto build_results = [&exit_times, &stop_times, &vertex_headways,
                        &train_trajectories](bool success) {
    return SimulatorResults{.success         = success,
                            .exit_times      = std::move(exit_times),
                            .stop_times      = std::move(stop_times),
                            .vertex_headways = vertex_headways,
                            .train_trajectories =
                                std::move(train_trajectories)};
  };

  const auto trains_on_edges = tr_on_edges();

  double t =
      next_train_timestep.empty() ? 0 : next_train_timestep.top().timestep;

  PLOGV << "Starting simulation...";

  // Simulate until no more trains are left to simulate
  while (!next_train_timestep.empty()) {
    auto [current_train_id, current_timestep] = next_train_timestep.top();
    const auto& current_train = train_list.get_train(current_train_id);
    next_train_timestep.pop();

    PLOGV << "Simulating at timestep " << current_timestep;

    // Advance all trains based on how much time has passed since last timestep
    for (const auto& tr : trains_in_network) {
      move_train(train_positions.at(tr), train_velocities.at(tr),
                 train_movements.at(tr), current_timestep - t);
    }
    t = current_timestep;

    // If train is yet to enter the network enter it now
    // TODO: Vertex headway constraints?
    // TODO: Check vertex order before entering?
    if (train_positions.at(current_train_id).front < 0 &&
        !get_train_edges_of_tr(current_train_id).empty()) {
      train_positions.at(current_train_id) = {
          -train_list.get_train(current_train_id).get_length(), 0};
      train_velocities.at(current_train_id) =
          get_instance()
              ->get_const_schedule(current_train_id)
              .get_initial_velocity();
      trains_in_network.insert(current_train_id);
    }

    // Append train trajectories
    for (const auto train_id : trains_in_network) {
      if (train_trajectories.at(train_id).contains(t))
        continue;
      PosVel pv = {train_positions.at(train_id).front,
                   train_velocities.at(train_id)};
      train_trajectories.at(train_id).emplace(t, pv);
    }

    // Check if train has reached the end of the route
    if (train_positions.at(current_train_id).front -
            train_edge_length(current_train_id) >
        -EPS) {
      trains_in_network.erase(current_train_id);
      trains_finished_simulating.emplace(current_train_id);
      continue;
    }

    // Find how much track ahead of the train is free, i.e. has the current
    //  train scheduled to be the next to use it
    auto [free_track, stop_type] = find_free_track_ahead(
        train_positions, trains_in_network, trains_on_edges, current_train_id);

    auto path_includes_exit = path_includes_exit_vertex(current_train_id);
    auto exit_velocity      = get_instance()
                             ->get_const_schedule(current_train_id)
                             .get_exit_velocity();
    if (path_includes_exit) {
      free_track.emplace_back(0, exit_velocity);
      // Reached end of line, add 'virtual' segment to allow train to leave
      //  network TODO: Is this needed? If yes also change 'reached end of line'
      //  logic
      // free_track.emplace_back(current_train.get_length(), exit_velocity);
    }

    auto next_movements = calculate_optimal_train_movements(
        free_track, train_velocities.at(current_train_id),
        train_list.get_train(current_train_id));

    // Check that train reaches exit target velocity
    if (path_includes_exit &&
        next_movements[next_movements.size() - 1].v != exit_velocity) {
      // Infeasible (I think) because cannot reach exit velocity.
      return build_results(false);
    }

    // TODO: train_movements.empty()?
    auto total_time = get_movement_total_time(next_movements);

    // If the stop at the end of the free track is not guaranteed (i.e. may
    //  disappear by the time needs to brake for it) then recalculate free
    //  track before the train starts braking for it.
    auto next_timestep =
        stop_type != FreeTrackStopType::TemporaryEnd
            ? t + total_time
            : t + total_time - next_movements.at(next_movements.size() - 1).t;

    next_train_timestep.push({current_train_id, next_timestep});
    train_movements.at(current_train_id) = next_movements;
  }

  PLOGV << "Simulation completed successfully.";

  assert(trains_finished_simulating.size() == number_of_trains);
  // Any failed states should have been detected already at this point
  return build_results(true);
}

// ----------------------------
// PRIVATE HELPER FUNCTIONS
// ----------------------------

std::pair<std::vector<cda_rail::simulator::EventSimulator::EdgeSegment>,
          cda_rail::simulator::EventSimulator::FreeTrackStopType>
cda_rail::simulator::EventSimulator::find_free_track_ahead(
    const std::vector<TrainPosition>&              train_positions,
    const std::unordered_set<size_t>&              trains_in_network,
    const std::vector<std::unordered_set<size_t>>& trains_on_edges,
    const size_t                                   current_train) const {
  const auto& current_pos = train_positions.at(current_train);
  const auto  current_edge =
      get_edge_at_position(current_train, current_pos.front);
  const auto& edges = get_train_edges_of_tr(current_train);

  const auto current_edge_index =
      std::ranges::find(edges, current_edge) - edges.begin();

  // For all edges ahead check if the edge is occupied or the following
  //  vertex/ttd must be traversed by a different train first - if either of
  //  these occur the point until which the track is free has been found
  std::vector<EdgeSegment> free_track{};

  for (size_t ei = current_edge_index; ei < edges.size(); ++ei) {
    const auto  edge_id = edges.at(ei);
    const auto& edge    = get_instance()->get_const_network().get_edge(edge_id);
    // TODO: many redundant calls to get_edge_position, calculate milestones
    //  once and reuse instead?
    for (size_t other_tr : trains_in_network) {
      if (other_tr == current_train ||
          !trains_on_edges.at(edge_id).contains(other_tr)) {
        continue;
      }
      // Train may be on edge
      const auto other_tr_edge_end   = get_edge_position(other_tr, edge_id);
      const auto other_tr_edge_start = other_tr_edge_end - edge.length;

      if (train_positions.at(other_tr).rear > other_tr_edge_start &&
          train_positions.at(other_tr).rear < other_tr_edge_end) {
        // Train rear is on the edge
        auto distance_into_edge =
            train_positions.at(other_tr).rear - other_tr_edge_start;
        free_track.push_back({distance_into_edge, edge.max_speed});
        // end of free track
        return {free_track, FreeTrackStopType::TemporaryEnd};
      }
      if (train_positions.at(other_tr).front > other_tr_edge_start &&
          train_positions.at(other_tr).front < other_tr_edge_end) {
        // Train front but not train rear is on the edge TODO

        // TODO: Handle case where other train is coming towards current train
        //  - is this automatically a deadlock since otherwise the previous
        //  vertex shouldn't allow this order?
      }
    }

    // Check whether the following vertex order dictates that a different
    //  train must pass through it first
    auto vertex       = edge.target;
    auto vertex_order = get_vertex_orders_of_vertex(vertex);
    for (const auto other_tr : vertex_order) {
      if (other_tr == current_train) {
        break;
      }
      // Check if train which must pass vertex first has already passed
      // If any train which must pass first has not yet passed, this is the
      //  end of the free track ahead
      auto vertex_position = get_edge_position(other_tr, edge_id);
      if (train_positions.at(other_tr).rear < vertex_position) {
        // end of free track TODO: How close to vertex should train stop
        free_track.push_back({edge.length * 0.9, edge.max_speed});
        return {free_track, FreeTrackStopType::TemporaryEnd};
      }
    }

    // This edge had no train on it - add to vector of free track ahead
    free_track.push_back({edge.length, edge.max_speed});

    // TODO: Check TTDs
    // TODO: also check reversed edge
  }
  return {free_track, FreeTrackStopType::GuaranteedEnd};
}

std::vector<cda_rail::simulator::EventSimulator::TrainMovement>
cda_rail::simulator::EventSimulator::calculate_optimal_train_movements(
    const std::vector<EdgeSegment>& free_track, const double initial_speed,
    const Train& train) {
  // TODO: Does not account for limit_speed_by_leaving_edges = true
  auto current_position = 0.0;
  auto current_speed    = initial_speed;
  auto movements        = std::vector<TrainMovement>{};

  for (int segment_index = 0; segment_index < free_track.size();
       ++segment_index) {
    // Get each segment ahead where the max speed is reduced relative to the
    //  previous segment, this is used to calculate the earliest braking point
    std::vector<MaxSpeedChangePoint> max_speed_reductions{};

    auto previous_max_speed  = -1.0;
    auto cumulative_position = 0.0;
    for (int ahead_segment_index = segment_index;
         ahead_segment_index < free_track.size(); ++ahead_segment_index) {
      const auto& segment = free_track.at(ahead_segment_index);
      if (segment.max_speed < previous_max_speed) {
        max_speed_reductions.push_back(
            {cumulative_position, segment.max_speed});
      }
      cumulative_position += segment.length;
      previous_max_speed = segment.max_speed;
    }

    // Find earliest braking point based on the points where max speed decreases
    // First find the max speed of the earliest braking point
    auto earliest_braking_point_speed_squared =
        std::numeric_limits<double>::max();
    auto earliest_braking_point_position = 0.0;
    auto braking_point_target_speed      = 0.0;

    if (max_speed_reductions.empty()) {
      // Max speed never decreases again, drive till the end as fast as possible
      append_max_permitted_speed_movements(
          movements,
          std::span{free_track}.subspan(segment_index,
                                        free_track.size() - segment_index),
          current_speed, train);
      return movements;
    }

    for (auto& [point_position, new_max_speed] : max_speed_reductions) {
      auto point_braking_point_speed_squared =
          max_speed_squared_two_phase_travel(
              current_speed, new_max_speed, train.get_acceleration(),
              train.get_deceleration(), point_position - current_position);
      if (point_braking_point_speed_squared <
          earliest_braking_point_speed_squared) {
        earliest_braking_point_speed_squared =
            point_braking_point_speed_squared;
        earliest_braking_point_position = point_position;
        braking_point_target_speed      = new_max_speed;
      }
    }

    // Calculate the train movements
    auto braking_point_position = distance_travelled_input_speed_squared(
        squared(current_speed), earliest_braking_point_speed_squared,
        train.get_acceleration());

    for (auto next_segment_index = segment_index;
         next_segment_index < free_track.size(); ++next_segment_index) {
      const auto& next_segment = free_track.at(next_segment_index);
      if (next_segment.length == 0)
        continue;
      const auto max_permitted_speed =
          std::min(next_segment.max_speed, train.get_max_speed());
      auto distance_to_max_permitted_speed = distance_travelled(
          current_speed, max_permitted_speed, train.get_acceleration());

      // Consider different cases based on the travel on this segment
      if (braking_point_position < next_segment.length) {
        if (earliest_braking_point_speed_squared <=
            squared(max_permitted_speed)) {
          // Case 1 - Braking point on edge and v <= permitted_speed
          auto earliest_braking_point_speed =
              std::sqrt(earliest_braking_point_speed_squared);
          auto distance_to_peak =
              distance_travelled(current_speed, earliest_braking_point_speed,
                                 train.get_acceleration());
          movements.emplace_back(current_speed, earliest_braking_point_speed,
                                 train.get_acceleration(), distance_to_peak);
          auto remaining_track      = next_segment.length - distance_to_peak;
          auto future_segment_index = next_segment_index;
          while (free_track.at(future_segment_index).max_speed >
                 braking_point_target_speed) {
            remaining_track += free_track.at(future_segment_index).length;
            ++future_segment_index;
            // These edge segments are being considered here, so they can be
            //  'skipped' for future iterations of the outer loops
            ++segment_index;
            ++next_segment_index;
          }
          movements.emplace_back(earliest_braking_point_speed,
                                 braking_point_target_speed,
                                 -train.get_deceleration(), remaining_track);
          current_position += braking_point_position;
          current_speed = braking_point_target_speed;
          break;
        }
        // Case 2 - Braking point on edge and v > permitted_speed
        movements.emplace_back(current_speed, max_permitted_speed,
                               train.get_acceleration(),
                               distance_to_max_permitted_speed);
        // Calculate new braking point from max permitted speed
        auto distance_to_brake_from_max_permitted =
            distance_travelled(max_permitted_speed, braking_point_target_speed,
                               -train.get_deceleration());
        auto coast_distance = earliest_braking_point_position -
                              distance_to_max_permitted_speed -
                              distance_to_brake_from_max_permitted;
        movements.emplace_back(max_permitted_speed, max_permitted_speed, 0.0,
                               coast_distance);
        movements.emplace_back(max_permitted_speed, braking_point_target_speed,
                               -train.get_deceleration(),
                               distance_to_brake_from_max_permitted);
        auto future_segment_index = next_segment_index + 1;
        while (free_track.at(future_segment_index).max_speed >
               braking_point_target_speed) {
          ++future_segment_index;
          // These edge segments are being considered here, so they can be
          //  'skipped' for future iterations of the outer loops
          ++segment_index;
          ++next_segment_index;
        }
        current_position += braking_point_position;
        current_speed = braking_point_target_speed;
        break;
      }
      if (distance_to_max_permitted_speed < next_segment.length) {
        // Case 3 - Braking point beyond edge and will exceed current speed
        //  limit
        movements.emplace_back(current_speed, max_permitted_speed,
                               train.get_acceleration(),
                               distance_to_max_permitted_speed);
        auto remaining_track =
            next_segment.length - distance_to_max_permitted_speed;
        movements.emplace_back(max_permitted_speed, max_permitted_speed, 0.0,
                               remaining_track);
        current_position += next_segment.length;
        current_speed = max_permitted_speed;
        break;
      }
      // Case 4 - Braking point on edge and won't exceed current speed limit
      auto speed_reached = std::sqrt(final_speed_squared(
          current_speed, train.get_acceleration(), next_segment.length));
      movements.emplace_back(current_speed, speed_reached,
                             train.get_acceleration(), next_segment.length);
      current_position += next_segment.length;
      current_speed = speed_reached;
      // No break - continue with next segment without recalculating v^2
    }
  }
  return movements;
}

void cda_rail::simulator::EventSimulator::append_max_permitted_speed_movements(
    std::vector<TrainMovement>&  movements,
    std::span<const EdgeSegment> free_track, double initial_speed,
    const Train& train) {
  auto current_speed = initial_speed;
  for (const auto& segment : free_track) {
    if (segment.length == 0)
      continue;
    auto max_permitted_speed =
        std::min(segment.max_speed, train.get_max_speed());
    auto distance_to_max_permitted_speed = distance_travelled(
        current_speed, max_permitted_speed, train.get_acceleration());
    if (distance_to_max_permitted_speed > segment.length) {
      auto speed_reached = std::sqrt(final_speed_squared(
          current_speed, train.get_acceleration(), segment.length));
      movements.emplace_back(current_speed, speed_reached,
                             train.get_acceleration(), segment.length);
      current_speed = speed_reached;
    } else {
      if (current_speed != max_permitted_speed) {
        movements.emplace_back(current_speed, max_permitted_speed,
                               train.get_acceleration(),
                               distance_to_max_permitted_speed);
      }
      movements.emplace_back(max_permitted_speed, max_permitted_speed, 0.0,
                             segment.length - distance_to_max_permitted_speed);
      current_speed = max_permitted_speed;
    }
  }
}

double cda_rail::simulator::EventSimulator::get_movement_total_time(
    const std::span<const TrainMovement>& movements) {
  auto total_time = 0.0;
  for (const auto& movement : movements) {
    total_time += movement.t;
  }
  return total_time;
}

void cda_rail::simulator::EventSimulator::move_train(
    TrainPosition& train_position, double& train_velocity,
    std::vector<TrainMovement>& train_movements, const double t) {
  auto   movement                = 0.0;
  auto   time_passed             = 0.0;
  double final_speed             = train_velocity;
  auto   fully_applied_movements = 0;
  for (auto& move : train_movements) {
    if (t - (time_passed + move.t) > -EPS) {
      movement += move.s;
      time_passed += move.t;
      fully_applied_movements++;
      final_speed = move.v;
      if (time_passed - t > -EPS)
        break;
    } else {
      // Only partially apply movement
      movement += distance_traveled_consume_movement(move, t - time_passed);
      final_speed = move.u;
      break;
    }
  }
  train_position.front += movement;
  train_position.rear += movement;
  train_velocity = final_speed;
  train_movements.erase(train_movements.begin(),
                        train_movements.begin() + fully_applied_movements);
}

double cda_rail::simulator::EventSimulator::distance_traveled_consume_movement(
    TrainMovement& movement, const double t) {
  assert(t < movement.t);
  const auto speed_reached = final_speed_time(movement.u, movement.a, t);
  const auto distance_traveled =
      distance_travelled_time(movement.u, speed_reached, t);
  movement.u = speed_reached;
  movement.s -= distance_traveled;
  movement.t -= t;
  return distance_traveled;
}

bool cda_rail::simulator::EventSimulator::path_includes_exit_vertex(
    const size_t train_id) const {
  const auto& train =
      get_instance()->get_const_train_list().get_train(train_id);
  const auto& path = get_train_edges_of_tr(train_id);
  if (path.empty()) {
    return false;
  }
  const auto  last_edge_id = path[path.size() - 1];
  const auto& last_edge =
      get_instance()->get_const_network().get_edge(last_edge_id);
  const auto exit_vertex_id =
      get_instance()->get_const_schedule(train_id).get_exit_vertex();
  return last_edge.target == exit_vertex_id;
}

double cda_rail::simulator::EventSimulator::squared(double x) { return x * x; }
