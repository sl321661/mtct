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
  std::unordered_set<size_t> trains_in_network;
  std::unordered_set<size_t> trains_left;
  std::unordered_set<size_t> trains_finished_simulating;
  std::vector<double>        vertex_headways(
      get_instance()->get_const_network().number_of_vertices(),
      0.0); // headways for vertices
  train_timestep_queue next_train_timestep;

  // Get first timestep for each train and detect trains that are not scheduled
  // to enter the network
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

  PLOGV << "Starting simulation from time " << t;

  // Simulate until no more trains are left to simulate
  while (!next_train_timestep.empty()) {
    auto current_train    = next_train_timestep.top().tr;
    auto current_timestep = next_train_timestep.top().timestep;
    next_train_timestep.pop();
    // TODO: Advance all trains based on t from previous timestep and new t

    // If train is yet to enter the network enter it now
    // TODO: Vertex headway constraints?
    // TODO: Check vertex order here?
    if (train_positions.at(current_train).front < 0 &&
        get_train_edges_of_tr(current_train).size() > 0) {
      train_positions.at(current_train) = {
          -train_list.get_train(current_train).get_length(), 0};
      train_velocities.at(current_train) =
          get_instance()
              ->get_const_schedule(current_train)
              .get_initial_velocity();
      trains_in_network.insert(current_train);
    }

    // Find how much track ahead of the train is free with the current train
    //  scheduled to be the next to use it
    // find_free_track_end(current_train, train_positions, trains_in_network,
    //                     trains_on_edges);
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
      const auto& edge = get_instance()->get_const_network().get_edge(edge_id);
      auto        found_end = false;
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
          found_end = true;
          break;
        }
        if (train_positions.at(other_tr).front > other_tr_edge_start &&
            train_positions.at(other_tr).front < other_tr_edge_end) {
          // Train front but not train rear is on the edge TODO

          // TODO: Handle case where other train is coming towards current train
          //  - is this automatically a deadlock since otherwise the previous
          //  vertex shouldn't allow this order?
        }
      }
      if (found_end) {
        // Found end of free track
        break;
      }
      // This edge had no train on it - add to vector of free track ahead

      // Check whether the following vertex order dictates that a different
      //  train must pass through it first
      auto vertex_occupied = false;
      auto vertex          = edge.target;
      auto vertex_order    = get_vertex_orders_of_vertex(vertex);
      for (const auto other_tr : vertex_order) {
        if (other_tr == current_train) {
          break;
        }
        // Check if train which must pass vertex first has already passed
        // If any train which must pass first has not yet passed, this is the
        //  end of the free track ahead
        auto vertex_position = get_edge_position(other_tr, edge_id);
        if (train_positions.at(other_tr).rear < vertex_position) {
          vertex_occupied = true;
          break;
        }
      }

      if (vertex_occupied) {
        // TODO: Where is the eps added again later on?
        free_track.push_back({edge.length - EPS, edge.max_speed});
        break;
      }
      free_track.push_back({edge.length, edge.max_speed});

      // TODO: Check TTDs
      // TODO: also check reversed edge
    }
  }

  assert(trains_finished_simulating.size() == number_of_trains);
  // Any failed states should have been detected already at this point
  return build_results(true);

  // int cycles_without_movement = 0;
  // while (cycles_without_movement < CycleLimit) {
  // PLOGV << "----------------------------";
  // PLOGV << "Current time: " << t;
  //
  // bool movement_detected = false;
  //
  // // Blocked Vertices, i.e., headway > t
  // cda_rail::index_set blocked_vertices;
  // for (size_t vertex = 0; vertex < vertex_headways.size(); ++vertex) {
  //   if (vertex_headways.at(vertex) > t) {
  //     blocked_vertices.insert(vertex);
  //   }
  // }
  //
  // // Move existing trains
  // for (const auto& tr : trains_in_network) {
  //   const auto& train_object = train_list.get_train(tr);
  //
  //   if (trains_finished_simulating.contains(tr) ||
  //       t < tr_stop_until.at(tr) + dt) {
  //     PLOGV << train_object.get_name() << " skipped.";
  //     continue;
  //   }
  //
  //   // Calculate MA
  //
  //   const auto tr_ma_data = get_ma_and_maxv(
  //       tr, train_velocities, tr_next_stop_id.at(tr), t, dt,
  //       blocked_vertices, train_positions, trains_in_network, trains_left,
  //       trains_on_edges, limit_speed_by_leaving_edges);
  //   PLOGV << train_object.get_name() << " positioned at "
  //         << train_positions.at(tr).front
  //         << " has MA: " << train_positions.at(tr).front + tr_ma_data.ma
  //         << " and max velocity: " << tr_ma_data.max_v;
  //
  //   auto tr_new_speed =
  //       std::min(tr_ma_data.max_v,
  //                get_v1_from_ma(train_velocities.at(tr), tr_ma_data.ma,
  //                               train_object.get_deceleration(), dt));
  //   if (tr_new_speed < V_MIN) {
  //     tr_new_speed = 0.0; // Train is stopped
  //   }
  //
  //   PLOGV << "tr_new_speed = " << tr_new_speed;
  //
  //   // Move trains
  //   if (move_train(tr, train_velocities.at(tr), tr_new_speed, tr_ma_data.ma,
  //                  dt, train_positions)) {
  //     movement_detected = true;
  //
  //     auto const tr_route_len = train_edge_length(tr);
  //     auto const last_edge_leaves_network =
  //         get_instance()
  //             ->get_const_network()
  //             .get_edge(get_train_edges_of_tr(tr).back())
  //             .target ==
  //         get_instance()->get_const_schedule(tr).get_exit_vertex();
  //     if (!last_edge_leaves_network &&
  //         tr_route_len < train_positions.at(tr).front + GRB_EPS) {
  //       PLOGV << "Train " << train_object.get_name()
  //             << " overshot the end of its route at "
  //             << train_positions.at(tr).front << " and is now stopped.";
  //       train_positions.at(tr).front = tr_route_len;
  //       tr_new_speed                 = 0.0;
  //     }
  //   }
  //   train_velocities.at(tr) = tr_new_speed;
  //   PLOGV << "At time " << t << ", " << train_object.get_name()
  //         << " moved to " << train_positions.at(tr).front << " with speed "
  //         << tr_new_speed << " and MA "
  //         << train_positions.at(tr).front +
  //                cda_rail::braking_distance(tr_new_speed,
  //                                           train_object.get_deceleration());
  //
  //   if (save_trajectories) {
  //     train_trajectories.at(tr)[t] = {.pos = train_positions.at(tr).front,
  //                                     .vel = train_velocities.at(tr)};
  //   }
  // }
  //
  // // Update rear positions of trains
  // update_rear_positions(train_positions);
  //
  // cda_rail::index_vector trains_to_remove;
  // for (const auto& tr : trains_in_network) {
  //   if (trains_finished_simulating.contains(tr)) {
  //     continue;
  //   }
  //
  //   // Remove trains that have left the network
  //   const auto tr_status = tr_reached_end(tr, train_positions);
  //   if (tr_status == DestinationType::Network) {
  //     trains_to_remove.emplace_back(tr);
  //     trains_left.insert(tr);
  //     trains_finished_simulating.insert(tr);
  //     exit_times.at(tr) = t;
  //     const auto& exit_vertex_idx =
  //         get_instance()->get_const_schedule(tr).get_exit_vertex();
  //     const auto& exit_vertex =
  //         get_instance()->get_const_network().get_vertex(exit_vertex_idx);
  //     vertex_headways.at(exit_vertex_idx) = t + exit_vertex.headway;
  //     PLOGV << "At time " << t << ", " << train_list.get_train(tr).get_name()
  //           << " left the network.";
  //   } else if (tr_status == DestinationType::Edge) {
  //     trains_finished_simulating.insert(tr);
  //     exit_times.at(tr) = t;
  //     PLOGV << "At time " << t << ", " << train_list.get_train(tr).get_name()
  //           << " reached the end of its route on an edge within the
  //           network.";
  //     if (disappear_at_partial_route_end) {
  //       trains_to_remove.emplace_back(tr);
  //       trains_left.insert(tr);
  //       PLOGV << "At time " << t << ", "
  //             << train_list.get_train(tr).get_name()
  //             << " disappeared at the end of its route.";
  //     }
  //   } else if (tr_status == DestinationType::Station) {
  //     assert(tr_next_stop_id.at(tr).has_value());
  //     const auto& last_stop =
  //         get_instance()->get_const_schedule(tr).get_stops().at(
  //             tr_next_stop_id.at(tr).value());
  //     exit_times.at(tr) = std::max(t + last_stop.get_service_duration(),
  //                                  last_stop.get_earliest_departure());
  //     stop_times.at(tr).emplace_back(t);
  //     tr_next_stop_id.at(tr) = {};
  //     trains_finished_simulating.insert(tr);
  //     PLOGV << "At time " << t << ", " << train_list.get_train(tr).get_name()
  //           << " reached the end of its route at station "
  //           << last_stop.get_station().name << ", stopping until "
  //           << exit_times.at(tr);
  //   } else {
  //     // Train is still in the network
  //     // Update stop information if a train has reached its next stop
  //     if (tr_next_stop_id.at(tr).has_value()) {
  //       if (train_velocities.at(tr) < V_MIN &&
  //           (train_positions.at(tr).front >=
  //            get_stop_positions_of_tr(tr).at(tr_next_stop_id.at(tr).value())
  //            -
  //                STOP_TOLERANCE)) {
  //         // Train has reached its next stop
  //         stop_times.at(tr).push_back(t);
  //         const auto& tr_stops =
  //             get_instance()->get_const_schedule(tr).get_stops();
  //         const auto& stop_info =
  //         tr_stops.at(tr_next_stop_id.at(tr).value()); tr_stop_until.at(tr) =
  //             std::max(t + stop_info.get_service_duration(),
  //                      stop_info.get_earliest_departure());
  //         PLOGV << "At time " << t << ", "
  //               << train_list.get_train(tr).get_name()
  //               << " reached its next stop at "
  //               << stop_info.get_station().name << ", stopping until "
  //               << tr_stop_until.at(tr);
  //
  //         if (tr_next_stop_id.at(tr).value() + 1 <
  //             get_stop_positions_of_tr(tr).size()) {
  //           // Update next stop ID
  //           tr_next_stop_id.at(tr) = tr_next_stop_id.at(tr).value() + 1;
  //           PLOGV << "Next stop: "
  //                 << tr_stops.at(tr_next_stop_id.at(tr).value())
  //                        .get_station()
  //                        .name;
  //         } else {
  //           // No more stops scheduled
  //           tr_next_stop_id.at(tr) = std::nullopt;
  //           PLOGV << "No more stops scheduled";
  //         }
  //       }
  //     }
  //   }
  // }
  //
  // // Remove trains that have left the network
  // for (const auto& tr : trains_to_remove) {
  //   trains_in_network.erase(tr);
  // }
  //
  // // Check for new trains entering the network
  // const auto [tr_to_enter_success, tr_to_enter] = get_entering_trains(
  //     t, trains_in_network, trains_left, trains_finished_simulating,
  //     late_entry_possible, dt);
  // if (!tr_to_enter_success) {
  //   PLOGV
  //       << "Simulation failed: Not all trains can enter the network at time "
  //       << t;
  //   return build_results(false);
  // }
  // for (const auto& tr : tr_to_enter) {
  //   const auto& train_schedule = get_instance()->get_const_schedule(tr);
  //   const auto& entry_vertex =
  //   get_instance()->get_const_network().get_vertex(
  //       train_schedule.get_entry_vertex());
  //   if (vertex_headways.at(train_schedule.get_entry_vertex()) > t) {
  //     PLOGV << "At time " << t << ", " << train_list.get_train(tr).get_name()
  //           << " cannot enter the network at " << entry_vertex.name
  //           << " due to vertex headway constraints until time "
  //           << vertex_headways.at(train_schedule.get_entry_vertex());
  //   } else if (!is_ok_to_enter(tr, train_positions, train_velocities,
  //                              trains_in_network, trains_on_edges)) {
  //     PLOGV << "At time " << t << ", " << train_list.get_train(tr).get_name()
  //           << " cannot enter the network at " << entry_vertex.name
  //           << " due to moving authority constraints constraints.";
  //   } else {
  //     trains_in_network.insert(tr);
  //     vertex_headways.at(train_schedule.get_entry_vertex()) =
  //         t + entry_vertex.headway;
  //     train_positions.at(tr)  = {.rear =
  //                                    -train_list.get_train(tr).get_length(),
  //                                .front = 0.0}; // Initialize positions
  //     train_velocities.at(tr) = train_schedule.get_initial_velocity();
  //     if (!get_stop_positions_of_tr(tr).empty()) {
  //       tr_next_stop_id.at(tr) = 0;
  //     }
  //     movement_detected = true;
  //     PLOGV << "At time " << t << ", " << train_list.get_train(tr).get_name()
  //           << " entered the network at " << entry_vertex.name;
  //     PLOGV << "New entry blocked until time "
  //           << vertex_headways.at(train_schedule.get_entry_vertex());
  //     if (save_trajectories) {
  //       train_trajectories.at(tr)[t] = {.pos = train_positions.at(tr).front,
  //                                       .vel = train_velocities.at(tr)};
  //     }
  //   }
  // }
  //
  // // Check if all trains have reached their destination
  // if (trains_finished_simulating.size() == number_of_trains) {
  //   PLOGV << "All trains have reached their destination at time " << t;
  //   return build_results(true);
  // }
  //
  // cycles_without_movement =
  //     movement_detected ? 0 : cycles_without_movement + 1;
  //
  // // Check if there might be a deadlock situation
  // if (!movement_detected) {
  //   // There might be a deadlock situation if none of the constraints depend
  //   // on time
  //   PLOGV << "No movement detected at time " << t;
  //   bool reason_found = false;
  //   if (std::ranges::any_of(vertex_headways, [t](int vertex_headway) {
  //         return vertex_headway > t;
  //       })) {
  //     PLOGV << "Vertex headway constraint prevents movement.";
  //     reason_found = true;
  //   } else if (std::ranges::any_of(tr_stop_until, [t, dt](int stop_time) {
  //                return stop_time + dt > t;
  //              })) {
  //     PLOGV << "Train stop constraint prevents movement.";
  //     reason_found = true;
  //   } else {
  //     for (size_t tr = 0; tr < train_velocities.size(); ++tr) {
  //       if (!trains_finished_simulating.contains(tr)) {
  //         if (trains_in_network.contains(tr)) {
  //           if ((get_instance()
  //                    ->get_const_network()
  //                    .get_edge(get_train_edges_of_tr(tr).back())
  //                    .target ==
  //                get_instance()->get_const_schedule(tr).get_exit_vertex()) &&
  //               get_instance()->get_const_schedule(tr).get_exit_time() > t) {
  //             PLOGV << train_list.get_train(tr).get_name()
  //                   << " is blocked by earliest exit.";
  //             reason_found = true;
  //             break;
  //           }
  //         } else {
  //           if (get_instance()->get_const_schedule(tr).get_entry_time() > t)
  //           {
  //             PLOGV << train_list.get_train(tr).get_name()
  //                   << " is blocked by earliest entry.";
  //             reason_found = true;
  //             break;
  //           }
  //         }
  //       }
  //     }
  //   }
  //   if (!reason_found) {
  //     PLOGV << "Trains are in a deadlock situation.";
  //     return build_results(false);
  //   }
  // }
  //
  // // Update time
  // t += dt;
  // }

  // throw std::runtime_error("Simulation failed: Cycle limit reached.");
}

// ----------------------------
// PRIVATE HELPER FUNCTIONS
// ----------------------------

void cda_rail::simulator::EventSimulator::find_free_track_end(
    const size_t current_train,
    const std::vector<cda_rail::simulator::EventSimulator::TrainPosition>&
                                                   train_positions,
    const std::unordered_set<size_t>&              trains_in_network,
    const std::vector<std::unordered_set<size_t>>& trains_on_edges) const {
  // auto current_pos  = train_positions.at(current_train);
  // auto current_edge = get_edge_at_position(current_train, current_pos.front);
  // auto edges        = get_train_edges_of_tr(current_train);
  //
  // auto current_edge_index =
  //     std::ranges::find(edges, current_edge) - edges.begin();
  //
  // // For all edges ahead check if the edge is occupied or the following
  // vertex
  // //  must be traversed by a different train first - if either of these occur
  // //  the point until which the track is free has been found
  // for (size_t ei = current_edge_index; ei < edges.size(); ++ei) {
  //   auto edge = get_instance()->get_const_network().get_edge(ei);
  //   for (size_t other_tr : trains_in_network) {
  //     if (other_tr == current_train ||
  //         !trains_on_edges.at(ei).contains(other_tr)) {
  //       continue;
  //     }
  //     // Train may be on edge
  //     auto milestones          = edge_milestones(other_tr);
  //     auto other_tr_edge_end   = get_edge_position(other_tr, ei);
  //     auto other_tr_edge_start = other_tr_edge_end - edge.length;
  //
  //     if (train_positions.at(other_tr).rear > other_tr_edge_start &&
  //         train_positions.at(other_tr).rear < other_tr_edge_end) {
  //       // Train rear is on the edge
  //
  //     } else if (train_positions.at(other_tr).front > other_tr_edge_start &&
  //                train_positions.at(other_tr).front < other_tr_edge_end) {
  //       // Train front is on the edge
  //
  //       // TODO: Handle case where other train is coming towards current
  //       train
  //       //  - is this automatically a deadlock since otherwise the previous
  //       //  vertex shouldn't allow this order?
  //     }
  //   }
  //
  //   auto vertex       = edge.target;
  //   auto vertex_order = get_vertex_orders_of_vertex(vertex);
  //   // TODO: Check TTDs
  //   // TODO: also check reversed edge
}
std::vector<cda_rail::simulator::EventSimulator::TrainMovement>
cda_rail::simulator::EventSimulator::calculate_optimal_train_movements(
    const std::vector<EdgeSegment>& free_track, const double initial_speed,
    const Train& train) const {
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

    if (max_speed_reductions.size() == 0) {
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
      const auto  max_permitted_speed =
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
          movements.push_back({current_speed, earliest_braking_point_speed,
                               train.get_acceleration(), distance_to_peak});
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
          movements.push_back({earliest_braking_point_speed,
                               braking_point_target_speed,
                               train.get_deceleration(), remaining_track});
          current_position += braking_point_position;
          current_speed = braking_point_target_speed;
          break;
        }
        // Case 2 - Braking point on edge and v > permitted_speed
        movements.push_back({current_speed, max_permitted_speed,
                             train.get_acceleration(),
                             distance_to_max_permitted_speed});
        // Calculate new braking point from max permitted speed
        auto distance_to_brake_from_max_permitted =
            distance_travelled(max_permitted_speed, braking_point_target_speed,
                               -train.get_deceleration());
        auto coast_distance = earliest_braking_point_position -
                              distance_to_max_permitted_speed -
                              distance_to_brake_from_max_permitted;
        movements.push_back(
            {max_permitted_speed, max_permitted_speed, 0.0, coast_distance});
        movements.push_back({max_permitted_speed, braking_point_target_speed,
                             train.get_deceleration(),
                             distance_to_brake_from_max_permitted});
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
        movements.push_back({current_speed, max_permitted_speed,
                             train.get_acceleration(),
                             distance_to_max_permitted_speed});
        auto remaining_track =
            next_segment.length - distance_to_max_permitted_speed;
        movements.push_back(
            {max_permitted_speed, max_permitted_speed, 0.0, remaining_track});
        current_position += next_segment.length;
        current_speed = max_permitted_speed;
        break;
      }
      // Case 4 - Braking point on edge and won't exceed current speed limit
      auto speed_reached = std::sqrt(final_speed_squared(
          current_speed, train.get_acceleration(), next_segment.length));
      movements.push_back({current_speed, speed_reached,
                           train.get_acceleration(), next_segment.length});
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
    auto max_permitted_speed =
        std::min(segment.max_speed, train.get_max_speed());
    auto distance_to_max_permitted_speed = distance_travelled(
        current_speed, max_permitted_speed, train.get_acceleration());
    if (distance_to_max_permitted_speed > segment.length) {
      auto speed_reached = std::sqrt(final_speed_squared(
          current_speed, train.get_acceleration(), segment.length));
      movements.push_back({current_speed, speed_reached,
                           train.get_acceleration(), segment.length});
      current_speed = speed_reached;
    } else {
      if (current_speed != max_permitted_speed) {
        movements.push_back({current_speed, max_permitted_speed,
                             train.get_acceleration(),
                             distance_to_max_permitted_speed});
      }
      movements.push_back({max_permitted_speed, max_permitted_speed, 0.0,
                           segment.length - distance_to_max_permitted_speed});
      current_speed = max_permitted_speed;
    }
  }
}

double cda_rail::simulator::EventSimulator::squared(double x) { return x * x; }
