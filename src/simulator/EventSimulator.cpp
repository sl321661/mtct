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
  std::vector<int>    train_next_stop_indices(number_of_trains, 0);
  std::vector<std::vector<TrainMovement>> train_movements(
      number_of_trains, std::vector<TrainMovement>());
  std::unordered_set<size_t> trains_in_network;
  std::unordered_set<size_t> trains_finished_simulating;
  std::vector<double>        vertex_headways(
      get_instance()->get_const_network().number_of_vertices(),
      0.0); // headways for vertices
  train_timestep_queue upcoming_train_timesteps;
  train_dependencies   train_dependencies(number_of_trains,
                                          std::vector<TrainDependence>{});

  // Get first timestep for each train, check for schedules with no stops, and
  //  detect trains that are not scheduled to enter the network
  for (size_t tr = 0; tr < number_of_trains; ++tr) {
    if (get_stop_positions_of_tr(tr).empty()) {
      train_next_stop_indices.at(tr) = -1;
    }

    if (get_train_edges_of_tr(tr).empty()) {
      trains_finished_simulating.insert(tr);
      continue;
    }

    auto const entry_time =
        get_instance()->get_const_schedule(tr).get_entry_time();
    upcoming_train_timesteps.push({tr, entry_time});
  }

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

  double t = upcoming_train_timesteps.empty()
                 ? 0
                 : upcoming_train_timesteps.top().timestep;

  PLOGV << "Starting simulation...";
  // TODO: Vertex headways not considered
  // TODO: Stations not considered

  // Simulate until no more trains are left to simulate
  while (!upcoming_train_timesteps.empty()) {
    auto [current_train_id, current_timestep] = upcoming_train_timesteps.top();
    const auto& current_train = train_list.get_train(current_train_id);
    upcoming_train_timesteps.pop();

    PLOGV << "Simulating at timestep " << current_timestep;

    // Advance all trains based on how much time has passed since last timestep
    for (const auto& tr : trains_in_network) {
      move_train(train_positions.at(tr), train_velocities.at(tr),
                 train_movements.at(tr), current_timestep - t);
    }
    t = current_timestep;

    // If train is yet to enter the network enter it now
    if (!trains_in_network.contains(current_train_id)) {
      // Train still needs to be entered
      // Check if train can be entered now, if not add a dependence to the train
      //  which is in the way
      const auto dependency =
          can_enter_train_or_get_dependency(train_positions, current_train_id);
      if (!dependency.has_value()) {
        enter_train(trains_in_network, train_positions, train_velocities,
                    current_train_id);
      } else {
        // Can't enter train due to vertex or TTD orders
        if (!late_entry_possible) {
          return build_results(false);
        }
        // Check if dependency will be cleared soon under already calculated
        //  movements
        const auto cleared = get_time_when_dependency_cleared(
            dependency.value(), train_movements.at(dependency->leading_tr),
            train_positions.at(dependency->leading_tr));
        if (cleared.has_value()) {
          // Time when dependency is cleared is already known
          // Add time as current train's next timestep
          upcoming_train_timesteps.push({current_train_id, cleared.value()});
        } else {
          train_dependencies.at(dependency->leading_tr)
              .push_back(dependency.value());
        }
        continue;
      }
    }

    // TODO: Need to be saved more often
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
      // TODO: disappear_at_partial_route_end not considered
      trains_in_network.erase(current_train_id);
      trains_finished_simulating.emplace(current_train_id);
      assert(train_dependencies.at(current_train_id).empty());
      continue;
    }

    // Find how much track ahead of the train is free, i.e. has the current
    //  train scheduled to be the next to use it
    auto [free_track, stop_type, track_dependency] = find_free_track_ahead(
        train_positions, trains_in_network, train_next_stop_indices,
        trains_on_edges, current_train_id);

    if (stop_type == FreeTrackStopType::StationStop ||
        stop_type == FreeTrackStopType::TemporaryEnd) {
      free_track.push_back({0, 0});
    }

    auto free_path_includes_exit = stop_type == FreeTrackStopType::LineEnd &&
                                   path_includes_exit_vertex(current_train_id);
    auto exit_velocity = get_instance()
                             ->get_const_schedule(current_train_id)
                             .get_exit_velocity();
    if (free_path_includes_exit) {
      free_track.emplace_back(0, exit_velocity);
    }

    // TODO: free_track total length == 0?

    auto next_movements = calculate_optimal_train_movements(
        free_track, train_velocities.at(current_train_id),
        train_list.get_train(current_train_id));

    // Check that train reaches exit target velocity
    if (free_path_includes_exit &&
        !approx_equal(next_movements[next_movements.size() - 1].v,
                      exit_velocity)) {
      // Infeasible (I think) because cannot reach exit velocity.
      return build_results(false);
    }

    // Handle service at station
    if (stop_type == FreeTrackStopType::StationStop) {
      const auto  si   = train_next_stop_indices.at(current_train_id);
      const auto& stop = get_instance()
                             ->get_const_schedule(current_train_id)
                             .get_stops()
                             .at(si);
      const auto service_start   = t + get_movement_total_time(next_movements);
      const auto time_to_service = stop.get_service_time() - service_start;
      const auto time_at_station =
          time_to_service > 0 ? time_to_service + stop.get_service_duration()
                              : stop.get_service_duration();
      next_movements.emplace_back(0.0, 0.0, 0.0, 0.0, time_at_station);

      if (get_stop_positions_of_tr(current_train_id).size() ==
          train_next_stop_indices.at(current_train_id) + 1) {
        // No more stops
        train_next_stop_indices.at(current_train_id) = -1;
      } else {
        train_next_stop_indices.at(current_train_id)++;
      }
    }

    auto total_time = get_movement_total_time(next_movements);

    // If the stop at the end of the free track is not guaranteed (i.e. may
    //  disappear by the time needs to brake for it) then recalculate free
    //  track before the train starts braking for it.
    auto next_timestep =
        stop_type != FreeTrackStopType::TemporaryEnd
            ? t + total_time
            : t + total_time - next_movements.at(next_movements.size() - 1).t;

    upcoming_train_timesteps.push({current_train_id, next_timestep});
    train_movements.at(current_train_id) = next_movements;

    // Check whether new movements clear another train's dependency
    for (auto d_i = 0; d_i < train_dependencies.at(current_train_id).size();
         d_i++) {
      const auto& dependency = train_dependencies.at(current_train_id).at(d_i);
      const auto  cleared_at = get_time_when_dependency_cleared(
          dependency, next_movements, train_positions.at(current_train_id));
      if (!cleared_at.has_value()) {
        continue;
      }
      upcoming_train_timesteps.push(
          {dependency.dependent_tr, t + cleared_at.value()});
      train_dependencies.at(current_train_id)
          .erase(train_dependencies.at(current_train_id).begin() + d_i);
    }
    // TODO: check for deadlock with dependencies
  }

  PLOGV << "Simulation completed successfully.";

  // TODO: Since trains may not enter the network if blocked and in this case
  //  never finish simulating but also never have a timestep anymore  this is no
  //  longer given and should be checked instead of asserted
  assert(trains_finished_simulating.size() == number_of_trains);
  // Any failed states should have been detected already at this point
  return build_results(true);
}

// ----------------------------
// PRIVATE HELPER FUNCTIONS
// ----------------------------

std::optional<cda_rail::simulator::EventSimulator::TrainDependence>
cda_rail::simulator::EventSimulator::can_enter_train_or_get_dependency(
    const std::vector<TrainPosition>& train_positions,
    const size_t                      train_id) const {
  // TODO: Check braking path of train with initial velocity - must also be free
  // Check entry vertex order
  const auto entry_vertex_id =
      get_instance()->get_const_schedule(train_id).get_entry_vertex();
  const auto& entry_vertex =
      get_instance()->get_const_network().get_vertex(entry_vertex_id);
  const auto& vertex_order = get_vertex_orders_of_vertex(entry_vertex);
  for (const auto other_tr : vertex_order) {
    if (train_id == other_tr) {
      break;
    }
    const auto train_vertex_position =
        get_vertex_position(other_tr, entry_vertex_id);
    if (train_positions.at(other_tr).rear >= train_vertex_position) {
      // Train has already passed this vertex
      continue;
    }
    return {{other_tr, train_id, train_vertex_position}};
  }
  // Check first edge TTD
  const auto& edges      = get_train_edges_of_tr(train_id);
  const auto  first_edge = edges.at(0);
  const auto  ttd        = get_ttd(first_edge);
  if (!ttd.has_value()) {
    return {};
  }
  const auto& ttd_order = get_ttd_orders_of_ttd(ttd.value());
  for (const auto other_tr : ttd_order) {
    if (train_id == other_tr) {
      break;
    }
    const auto train_ttd_position = train_ttd_end_position(
        train_positions.at(other_tr), other_tr, ttd.value());
    if (train_positions.at(other_tr).rear >= train_ttd_position) {
      // Train has already passed this TTD
      continue;
    }
    return {{other_tr, train_id, train_ttd_position}};
  }
  return {};
}
std::optional<double>
cda_rail::simulator::EventSimulator::get_time_when_dependency_cleared(
    const TrainDependence&            dependency,
    const std::vector<TrainMovement>& leading_movements,
    const TrainPosition&              leading_position) {
  double pos  = leading_position.rear; // Use rear so location is fully cleared
  double time = 0.0;
  for (const auto& movement : leading_movements) {
    if (pos + movement.s >= dependency.required_location) {
      return {time + time_taken(movement.u, movement.v,
                                dependency.required_location - pos)};
    }
    pos += movement.s;
    time += movement.t;
  }
  // Doesn't clear the dependency with movements known at this time
  return {};
}

void cda_rail::simulator::EventSimulator::enter_train(
    std::unordered_set<size_t>& trains_in_network,
    std::vector<TrainPosition>& train_positions,
    std::vector<double>& train_velocities, size_t train_id) const {
  train_positions.at(train_id) = {
      -get_instance()->get_const_train_list().get_train(train_id).get_length(),
      0};
  train_velocities.at(train_id) =
      get_instance()->get_const_schedule(train_id).get_initial_velocity();
  trains_in_network.insert(train_id);
}

std::tuple<std::vector<cda_rail::simulator::EventSimulator::EdgeSegment>,
           cda_rail::simulator::EventSimulator::FreeTrackStopType,
           std::optional<std::pair<
               cda_rail::simulator::EventSimulator::TrainDependence,
               cda_rail::simulator::EventSimulator::FreeTrackDependencyType>>>
cda_rail::simulator::EventSimulator::find_free_track_ahead(
    const std::vector<TrainPosition>&              train_positions,
    const std::unordered_set<size_t>&              trains_in_network,
    const std::vector<int>&                        next_stop_indices,
    const std::vector<std::unordered_set<size_t>>& trains_on_edges,
    const size_t                                   current_train) const {
  auto       front_pos       = train_positions.at(current_train).front;
  const auto current_edge_id = get_edge_at_position(current_train, front_pos);
  // Distance into current edge is needed so only the remaining part of the edge
  //  is considered for the free track ahead
  const auto current_edge_length =
      get_instance()->get_const_network().get_edge(current_edge_id).length;
  auto current_edge_start_position =
      get_edge_position(current_train, current_edge_id) - current_edge_length;
  auto distance_into_current_edge = front_pos - current_edge_start_position;
  const auto& edges               = get_train_edges_of_tr(current_train);

  auto current_edge_index =
      std::ranges::find(edges, current_edge_id) - edges.begin();

  // For all edges ahead check if the edge is occupied or the following
  //  vertex/ttd must be traversed by a different train first - if either of
  //  these occur the point until which the track is free has been found
  std::vector<EdgeSegment> free_track{};

  // TODO: many redundant calls to get_edge_position, calculate milestones
  //  once and reuse instead?

  for (size_t ei = current_edge_index; ei < edges.size(); ++ei) {
    // TODO: Also check reversed edge
    const auto  edge_id = edges.at(ei);
    const auto& edge    = get_instance()->get_const_network().get_edge(edge_id);
    // Check edge TTD
    const auto ttd = get_ttd(edge_id);
    if (!ttd.has_value()) {
      const auto& ttd_order = get_ttd_orders_of_ttd(ttd.value());
      for (const auto other_tr : ttd_order) {
        if (current_train == other_tr) {
          break;
        }
        const auto train_ttd_position = train_ttd_end_position(
            train_positions.at(other_tr), other_tr, ttd.value());
        if (train_positions.at(other_tr).rear >= train_ttd_position) {
          // Train has already passed this TTD
          continue;
        }
        // TTD is blocked, train shouldn't already be on edge
        assert(distance_into_current_edge == 0.0);
        return {free_track,
                FreeTrackStopType::TemporaryEnd,
                {{{other_tr, current_train, train_ttd_position},
                  FreeTrackDependencyType::Pass}}};
      }
    }

    // Check if there is a train on the edge
    for (const auto other_tr_id : trains_in_network) {
      auto other_tr =
          get_instance()->get_const_train_list().get_train(other_tr_id);
      if (other_tr_id == current_train ||
          !trains_on_edges.at(edge_id).contains(other_tr_id)) {
        continue;
      }
      // Train may be on edge
      const auto other_tr_edge_end   = get_edge_position(other_tr_id, edge_id);
      const auto other_tr_edge_start = other_tr_edge_end - edge.length;

      if (other_tr_edge_start >= train_positions.at(other_tr_id).rear &&
          other_tr_edge_start <= train_positions.at(other_tr_id).front) {
        // Some part of train is on the start of the edge
        return {free_track,
                FreeTrackStopType::TemporaryEnd,
                {{{other_tr_id, current_train,
                   other_tr_edge_start + other_tr.get_length()},
                  FreeTrackDependencyType::Pass}}};
      }

      if (train_positions.at(other_tr_id).rear >= other_tr_edge_start &&
          train_positions.at(other_tr_id).rear <= other_tr_edge_end) {
        // Train rear is on the edge
        auto other_tr_distance_into_edge =
            train_positions.at(other_tr_id).rear - other_tr_edge_start;
        free_track.push_back(
            {other_tr_distance_into_edge - distance_into_current_edge,
             edge.max_speed});
        // end of free track
        return {free_track,
                FreeTrackStopType::TemporaryEnd,
                {{{other_tr_id, current_train,
                   other_tr_edge_start + other_tr.get_length()},
                  FreeTrackDependencyType::Follow}}};
      }
    }

    // No other train on edge - check if end of edge is a scheduled stop for
    //  this train
    if (auto si = next_stop_indices.at(current_train); si != -1) {
      auto next_stop_pos = get_stop_positions_of_tr(current_train).at(si);
      if (front_pos + edge.length - next_stop_pos > -EPS) {
        // TODO: stop distance from edge end?
        free_track.emplace_back(edge.length - distance_into_current_edge,
                                edge.max_speed);
        return {free_track, FreeTrackStopType::StationStop, {}};
      }
    }

    // Check whether the following vertex order dictates that a different
    //  train must pass through it first
    auto vertex_id    = edge.target;
    auto vertex_order = get_vertex_orders_of_vertex(vertex_id);
    for (const auto other_tr : vertex_order) {
      if (other_tr == current_train) {
        break;
      }
      // Check if train which must pass vertex first has already passed
      // If any train which must pass first has not yet passed, this is the
      //  end of the free track ahead
      auto vertex_position = get_vertex_position(other_tr, vertex_id);
      if (train_positions.at(other_tr).rear < vertex_position) {
        free_track.push_back(
            {edge.length - distance_into_current_edge, edge.max_speed});
        // end of free track
        return {free_track,
                FreeTrackStopType::TemporaryEnd,
                {{{other_tr, current_train, vertex_position},
                  FreeTrackDependencyType::Pass}}};
      }
    }

    // This edge had no train on it - add to vector of free track ahead
    // If this is the first edge being evaluated, only add the amount the train
    //  has not already traversed
    free_track.push_back(
        {edge.length - distance_into_current_edge, edge.max_speed});

    front_pos += edge.length;
    // For all edges after the first one the current train cannot have partially
    //  traversed it already
    distance_into_current_edge = 0.0;

    // TODO: also check reversed edge
  }
  return {free_track, FreeTrackStopType::LineEnd, {}};
}

std::vector<cda_rail::simulator::EventSimulator::TrainMovement>
cda_rail::simulator::EventSimulator::calculate_optimal_train_movements(
    std::vector<EdgeSegment>& free_track, const double initial_speed,
    const Train& train) {
  // TODO: Does not account for limit_speed_by_leaving_edges = true
  limit_segment_speeds(free_track, train.get_max_speed());

  auto current_position = 0.0;
  auto current_speed    = initial_speed;
  auto movements        = std::vector<TrainMovement>{};

  for (int segment_index = 0; segment_index < free_track.size();
       ++segment_index) {
    auto next_speed_restriction = find_first_speed_restriction(
        free_track, segment_index, train, current_speed);

    if (!next_speed_restriction.has_value()) {
      // Max speed never decreases again, drive till the end as fast as possible
      append_max_permitted_speed_movements(
          movements,
          std::span{free_track}.subspan(segment_index,
                                        free_track.size() - segment_index),
          current_speed, train);
      return movements;
    }

    auto [speed_restriction, restriction_peak_speed_squared] =
        next_speed_restriction.value();

    // Calculate the train movements

    auto braking_point_position = distance_travelled_input_speed_squared(
        squared(current_speed), restriction_peak_speed_squared,
        train.get_acceleration());

    for (auto next_segment_index = segment_index;
         next_segment_index < free_track.size(); ++next_segment_index) {
      const auto& next_segment = free_track.at(next_segment_index);
      if (next_segment.length == 0)
        continue;
      const auto max_permitted_speed =
          std::min(next_segment.max_speed, train.get_max_speed());
      auto distance_to_max_speed = distance_travelled(
          current_speed, next_segment.max_speed, train.get_acceleration());

      // Consider different cases based on the travel on this segment
      if (braking_point_position < next_segment.length) {
        if (restriction_peak_speed_squared <= squared(next_segment.max_speed)) {
          // Case 1 - Braking point on edge and v <= permitted_speed
          auto earliest_braking_point_speed =
              std::sqrt(restriction_peak_speed_squared);
          auto distance_to_peak =
              distance_travelled(current_speed, earliest_braking_point_speed,
                                 train.get_acceleration());
          movements.emplace_back(current_speed, earliest_braking_point_speed,
                                 train.get_acceleration(), distance_to_peak);
          auto remaining_track       = next_segment.length - distance_to_peak;
          auto skipped_segment_index = next_segment_index;
          while (free_track.at(skipped_segment_index).max_speed >
                 speed_restriction.new_max_speed) {
            remaining_track += free_track.at(skipped_segment_index).length;
            ++skipped_segment_index;
            // These edge segments are being considered here, so they can be
            //  'skipped' for future iterations of the outer loops
            ++segment_index;
            ++next_segment_index;
          }
          movements.emplace_back(earliest_braking_point_speed,
                                 speed_restriction.new_max_speed,
                                 -train.get_deceleration(), remaining_track);
          current_position += distance_to_peak + remaining_track;
          current_speed = speed_restriction.new_max_speed;
          break;
        }
        // Case 2 - Braking point on edge and v > permitted_speed
        if (!approx_equal(current_speed, next_segment.max_speed))
          movements.emplace_back(current_speed, next_segment.max_speed,
                                 train.get_acceleration(),
                                 distance_to_max_speed);
        // Calculate new braking point from max permitted speed
        auto distance_to_brake_from_max_permitted = distance_travelled(
            next_segment.max_speed, speed_restriction.new_max_speed,
            -train.get_deceleration());
        auto coast_distance = speed_restriction.position -
                              distance_to_max_speed -
                              distance_to_brake_from_max_permitted;
        movements.emplace_back(next_segment.max_speed, next_segment.max_speed,
                               0.0, coast_distance);
        movements.emplace_back(
            next_segment.max_speed, speed_restriction.new_max_speed,
            -train.get_deceleration(), distance_to_brake_from_max_permitted);
        auto skipped_segment_index = next_segment_index + 1;
        auto skipped_distance      = 0.0;
        while (free_track.at(skipped_segment_index).max_speed >
               speed_restriction.new_max_speed) {
          skipped_distance += free_track.at(skipped_segment_index).length;
          ++skipped_segment_index;
          // These edge segments are being considered here, so they can be
          //  'skipped' for future iterations of the outer loops
          ++segment_index;
          ++next_segment_index;
        }
        current_position += next_segment.length + skipped_distance;
        current_speed = speed_restriction.new_max_speed;
        break;
      }
      if (distance_to_max_speed < next_segment.length) {
        // Case 3 - Braking point beyond edge and will exceed current speed
        //  limit
        movements.emplace_back(current_speed, next_segment.max_speed,
                               train.get_acceleration(), distance_to_max_speed);
        auto remaining_track = next_segment.length - distance_to_max_speed;
        movements.emplace_back(next_segment.max_speed, next_segment.max_speed,
                               0.0, remaining_track);
        current_position += next_segment.length;
        current_speed = next_segment.max_speed;
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

std::optional<
    std::pair<cda_rail::simulator::EventSimulator::MaxSpeedChangePoint, double>>
cda_rail::simulator::EventSimulator::find_first_speed_restriction(
    const std::vector<EdgeSegment>& free_track, int init_segment_index,
    const Train& train, double current_speed) {
  // Find all points where the speed limit is reduced
  std::vector<MaxSpeedChangePoint> max_speed_reductions{};
  auto                             previous_max_speed  = -1.0;
  auto                             cumulative_position = 0.0;
  for (int segment_index = init_segment_index;
       segment_index < free_track.size(); ++segment_index) {
    const auto& segment = free_track.at(segment_index);
    if (segment.max_speed < previous_max_speed) {
      max_speed_reductions.emplace_back(cumulative_position, segment.max_speed);
    }
    cumulative_position += segment.length;
    previous_max_speed = segment.max_speed;
  }

  if (max_speed_reductions.empty())
    return {};

  // Of the points where the speed limit is reduced, find the first one which
  //   needs to be considered by this train
  auto peak_speed_squared      = std::numeric_limits<double>::max();
  auto restriction_position    = 0.0;
  auto restriction_speed_limit = 0.0;
  for (auto& [point_position, new_max_speed] : max_speed_reductions) {
    auto point_peak_speed_squared = max_speed_squared_two_phase_travel(
        current_speed, new_max_speed, train.get_acceleration(),
        train.get_deceleration(), point_position);
    if (point_peak_speed_squared < peak_speed_squared) {
      peak_speed_squared      = point_peak_speed_squared;
      restriction_position    = point_position;
      restriction_speed_limit = new_max_speed;
    }
  }

  return {
      {{restriction_position, restriction_speed_limit}, peak_speed_squared}};
}

void cda_rail::simulator::EventSimulator::limit_segment_speeds(
    std::vector<EdgeSegment>& segments, double speed_limit) {
  for (auto& s : segments) {
    if (s.max_speed > speed_limit) {
      s.max_speed = speed_limit;
    }
  }
}

double cda_rail::simulator::EventSimulator::get_shared_track_ahead_distance(
    TrainPosition& tr1_pos, size_t tr1, size_t tr2) const {
  const auto& tr1_edges = get_train_edges_of_tr(tr1);
  const auto& tr2_edges = get_train_edges_of_tr(tr2);

  // Find start edge
  size_t tr1_edge_index     = 0;
  size_t tr1_edge_id        = 0;
  auto   distance_into_edge = 0.0;
  auto   current_pos        = 0.0;
  for (auto i = 0; i < tr1_edges.size(); ++i) {
    const auto  edge_id = tr1_edges.at(i);
    const auto& edge    = get_instance()->get_const_network().get_edge(edge_id);
    current_pos += edge.length;
    if (tr1_pos.front < current_pos) {
      tr1_edge_index     = i;
      tr1_edge_id        = edge_id;
      distance_into_edge = tr1_pos.front - current_pos + edge.length;
      break;
    }
  }

  double shared_distance = 0.0;
  auto   it              = std::ranges::find(tr2_edges, tr1_edge_id);
  if (it == tr2_edges.end()) {
    // No (more) track is shared between trains
    return 0.0;
  }
  auto tr2_edge_index = std::distance(tr2_edges.begin(), it);
  for (auto i = tr2_edge_index; i < tr2_edges.size(); ++i) {
    if (tr1_edges.at(tr1_edge_index) != tr2_edges.at(i)) {
      return shared_distance;
    }
    tr1_edge_index++;

    if (tr1_edge_id == tr2_edges.at(i)) {
      // First shared edge, only add distance that has not been traversed by tr1
      shared_distance += distance_into_edge;
    } else {
      const auto& edge =
          get_instance()->get_const_network().get_edge(tr2_edges.at(i));
      shared_distance += edge.length;
    }
  }
  return shared_distance;
}

void cda_rail::simulator::EventSimulator::append_max_permitted_speed_movements(
    std::vector<TrainMovement>&  movements,
    std::span<const EdgeSegment> free_track, double initial_speed,
    const Train& train) {
  auto current_speed = initial_speed;
  for (const auto& segment : free_track) {
    if (segment.length == 0)
      continue;
    auto distance_to_max_permitted_speed = distance_travelled(
        current_speed, segment.max_speed, train.get_acceleration());
    if (distance_to_max_permitted_speed > segment.length) {
      auto speed_reached = std::sqrt(final_speed_squared(
          current_speed, train.get_acceleration(), segment.length));
      movements.emplace_back(current_speed, speed_reached,
                             train.get_acceleration(), segment.length);
      current_speed = speed_reached;
    } else {
      if (!approx_equal(current_speed, segment.max_speed)) {
        movements.emplace_back(current_speed, segment.max_speed,
                               train.get_acceleration(),
                               distance_to_max_permitted_speed);
      }
      movements.emplace_back(segment.max_speed, segment.max_speed, 0.0,
                             segment.length - distance_to_max_permitted_speed);
      current_speed = segment.max_speed;
    }
  }
}

void cda_rail::simulator::EventSimulator::append_movement_if_not_zero(
    std::vector<TrainMovement>& movements, const TrainMovement& m) {
  if (m.s == 0)
    return;
  movements.push_back(m);
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

double cda_rail::simulator::EventSimulator::get_vertex_position(
    const size_t train_id, const size_t vertex_id) const {
  const auto& edges    = get_train_edges_of_tr(train_id);
  auto        position = 0.0;
  if (!edges.empty() &&
      get_instance()->get_const_network().get_edge(edges[0]).source ==
          vertex_id) {
    return position;
  }
  for (const auto edge_id : edges) {
    const auto& edge = get_instance()->get_const_network().get_edge(edge_id);
    position += edge.length;
    if (edge.target == vertex_id) {
      return position;
    }
  }

  const auto& vertex =
      get_instance()->get_const_network().get_vertex(vertex_id);
  throw exceptions::ConsistencyException(
      "Vertex " + vertex.name + " not found in train " +
      std::to_string(train_id) + "'s route.");
}

bool cda_rail::simulator::EventSimulator::path_includes_exit_vertex(
    const size_t train_id) const {
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

cda_rail::index_set
cda_rail::simulator::EventSimulator::get_passed_edges_of_train(
    const TrainPosition& train_position, const size_t train_id) const {
  // TODO: Ensure trains which have exited the network have a position beyond
  //  the last edge so all edges are returned in this function
  const auto& edges = get_train_edges_of_tr(train_id);
  index_set   passed_edges;
  double      distance = 0.0;
  for (const auto edge_id : edges) {
    distance += get_instance()->get_const_network().get_edge(edge_id).length;
    if (distance < train_position.rear)
      passed_edges.insert(edge_id);
    else
      break;
  }
  return passed_edges;
}

bool cda_rail::simulator::EventSimulator::train_has_passed_vertex(
    const TrainPosition& train_position, const size_t train_id,
    const size_t vertex_id) const {
  const auto pos = get_vertex_position(train_id, vertex_id);
  return train_position.rear > pos;
}

double cda_rail::simulator::EventSimulator::train_ttd_end_position(
    const TrainPosition& train_position, size_t train_id, size_t ttd_id) const {
  const auto& train_edges      = get_train_edges_of_tr(train_id);
  const auto& ttd_edges        = get_ttd_sections().at(ttd_id);
  auto        last_position    = -1.0;
  auto        current_position = 0.0;
  for (const auto edge_id : train_edges) {
    const auto& edge = get_instance()->get_const_network().get_edge(edge_id);
    current_position += edge.length;
    if (std::ranges::find(ttd_edges, edge_id) != ttd_edges.end()) {
      // Edge is part of ttd
      last_position = current_position;
    }
  }
  if (last_position == -1.0) {
    throw exceptions::InvalidInputException(
        "Given train does not cross given ttd.");
  }
  return last_position;
}

bool cda_rail::simulator::EventSimulator::train_has_passed_ttd(
    const TrainPosition& train_position, const size_t train_id,
    const size_t ttd_id) const {
  const auto& ttd_edges   = get_ttd_sections().at(ttd_id);
  const auto& train_edges = get_train_edges_of_tr(train_id);
  const auto& edges_passed =
      get_passed_edges_of_train(train_position, train_id);
  for (const auto edge_id : ttd_edges) {
    if (std::ranges::find(train_edges, edge_id) != train_edges.end() &&
        std::ranges::find(edges_passed, edge_id) == edges_passed.end()) {
      // Edge is traversed by train at some point but has not yet been passed
      return false;
    }
  }
  return true;
}

double cda_rail::simulator::EventSimulator::squared(double x) { return x * x; }
