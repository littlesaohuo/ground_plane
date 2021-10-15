#include "linefit.h"

Bin::Bin() : has_point_(false), min_z_(numeric_limits<double>::max()) {}

void Bin::AddPoint(const double& d, const double& z,
                   const pair<int, int>& coordinate) {
  has_point_ = true;
  if (z < min_z_) {
    min_z_range_ = d;
    min_z_ = z;
    coordinate_ = coordinate;
  }
}

Bin::MinZPoint Bin::GetMinZPoint() {
  MinZPoint point;
  if (has_point_) {
    point.d = min_z_range_;
    point.z = min_z_;
    point.coordinate = coordinate_;
  }
  return point;
}

Segment::Segment(const GroundSegmentationParams& params)
    : params_(params),
      bins_(params.n_bins),
      start_bin_index_(numeric_limits<int>::max()),
      end_bin_index_(numeric_limits<int>::min()) {}

void Segment::FitSegmentLines() {
  if (start_bin_index_ > end_bin_index_) return;
  lines_.clear();
  auto line_start = bins_.begin() + start_bin_index_;
  auto line_end = bins_.begin() + end_bin_index_ + 1;
  bool is_long_line = false;
  double cur_ground_height = -params_.sensor_height;
  list<Bin::MinZPoint> current_line_points(1, line_start->GetMinZPoint());
  LocalLine cur_line = make_pair(0, 0);
  for (auto line_iter = line_start + 1; line_iter != line_end; ++line_iter) {
    if (line_iter->HasPoint()) {
      Bin::MinZPoint cur_point = line_iter->GetMinZPoint();
      if (cur_point.d - current_line_points.back().d > params_.long_threshold)
        is_long_line = true;
      if (current_line_points.size() >= 2) {
        // Get expected z value to possibly reject far away points.
        double expected_z = numeric_limits<double>::max();
        if (is_long_line && current_line_points.size() > 2) {
          expected_z = cur_line.first * cur_point.d + cur_line.second;
        }
        current_line_points.push_back(cur_point);
        cur_line = FitLocalLine(current_line_points);
        const double error = GetMeanError(current_line_points, cur_line);
        // Check if not a good line.
        if (error > params_.max_error_square ||
            fabs(cur_line.first) > params_.max_slope ||
            is_long_line &&
                fabs(expected_z - cur_point.z) > params_.max_long_height) {
          // Add line until previous point as ground.
          current_line_points.pop_back();
          // Don't let lines with 2 base points through.
          if (current_line_points.size() >= 3) {
            const LocalLine new_line = FitLocalLine(current_line_points);
            lines_.push_back(LocalLineToLine(new_line, current_line_points));
            cur_ground_height =
                new_line.first * current_line_points.back().d + new_line.second;
          }
          // Start new line.
          is_long_line = false;
          current_line_points.erase(current_line_points.begin(),
                                    --current_line_points.end());
          --line_iter;
        }
        // Good line, continue.
        else {
        }
      } else {
        // Not enough points.
        if (cur_point.d - current_line_points.back().d <
                params_.long_threshold &&
            fabs(current_line_points.back().z - cur_ground_height) <
                params_.max_start_height) {
          // Add point if valid.
          current_line_points.push_back(cur_point);
        } else {
          // Start new line.
          current_line_points.clear();
          current_line_points.push_back(cur_point);
        }
      }
    }
  }
  // Add last line.
  if (current_line_points.size() > 2) {
    const LocalLine new_line = FitLocalLine(current_line_points);
    lines_.push_back(LocalLineToLine(new_line, current_line_points));
  }
}

double Segment::VerticalDistanceToLine(const double& d, const double& z) {
  const double kMargin = params_.line_expand_range;
  double distance = -1;
  for (auto it = lines_.begin(); it != lines_.end(); ++it) {
    if (it->first.d - kMargin < d && it->second.d + kMargin > d) {
      const double delta_z = it->second.z - it->first.z;
      const double delta_d = it->second.d - it->first.d;
      const double expected_z =
          (d - it->first.d) / delta_d * delta_z + it->first.z;
      distance = fabs(z - expected_z);
    }
  }
  return distance;
}

double Segment::GetMaxError(const list<Bin::MinZPoint>& points,
                            const LocalLine& line) {
  double max_error = 0;
  for (auto it = points.begin(); it != points.end(); ++it) {
    const double residual = (line.first * it->d + line.second) - it->z;
    const double error = residual * residual;
    if (error > max_error) max_error = error;
  }
  return max_error;
}

double Segment::GetMeanError(const list<Bin::MinZPoint>& points,
                             const LocalLine& line) {
  double error_sum = 0;
  for (auto it = points.begin(); it != points.end(); ++it) {
    const double residual = (line.first * it->d + line.second) - it->z;
    error_sum += residual * residual;
  }
  return error_sum / points.size();
}

Segment::LocalLine Segment::FitLocalLine(const list<Bin::MinZPoint>& points) {
  const unsigned int n_points = points.size();
  Eigen::MatrixXd X(n_points, 2);
  Eigen::VectorXd Y(n_points);
  unsigned int counter = 0;
  for (auto iter = points.begin(); iter != points.end(); ++iter) {
    X(counter, 0) = iter->d;
    X(counter, 1) = 1;
    Y(counter) = iter->z;
    ++counter;
  }
  Eigen::VectorXd result = X.colPivHouseholderQr().solve(Y);
  LocalLine line_result;
  line_result.first = result(0);
  line_result.second = result(1);
  return line_result;
}

Segment::Line Segment::LocalLineToLine(
    const LocalLine& local_line, const list<Bin::MinZPoint>& line_points) {
  Line line;
  const double first_d = line_points.front().d;
  const double second_d = line_points.back().d;
  const double first_z = local_line.first * first_d + local_line.second;
  const double second_z = local_line.first * second_d + local_line.second;
  line.first.z = first_z;
  line.first.d = first_d;
  line.second.z = second_z;
  line.second.d = second_d;
  return line;
}

GroundSegmenter::GroundSegmenter(const GroundSegmentationParams& params)
    : params_(params) {
  segment_step_ = 2 * M_PI / params.n_segments;
  bin_step_ = (params.r_max - params.r_min) / params.n_bins;
}

void GroundSegmenter::Excute(const pcl::PointCloud<pcl::PointXYZ>& cloud,
                             vector<uint8_t>& labels) {
  labels.clear();
  labels.resize(cloud.size(), 0);
  segments_.clear();
  segments_.resize(params_.n_segments, Segment(params_));
  segment_coordinates_.clear();
  segment_coordinates_.resize(cloud.size());

  InsertPoints(cloud, labels);
  GetLines();
  AssignCluster(labels);
}

void GroundSegmenter::InsertPoints(const pcl::PointCloud<pcl::PointXYZ>& cloud,
                                   vector<uint8_t>& labels) {
  for (size_t i = 0; i < cloud.size(); i++) {
    const auto& point = cloud[i];
    const double range = sqrt(point.x * point.x + point.y * point.y);
    if (range <= params_.r_min || range >= params_.r_max) {
      labels.at(i) = 3;
      continue;
    }
    const double angle = atan2(point.y, point.x);
    int segment_index = floor((angle + M_PI) / segment_step_);
    segment_index =
        ((segment_index >= params_.n_segments) ? params_.n_segments - 1
                                               : segment_index);
    segment_index = ((segment_index < 0) ? 0 : segment_index);
    int bin_index = floor((range - params_.r_min) / bin_step_);
    bin_index = (bin_index >= params_.n_bins) ? params_.n_bins - 1 : bin_index;
    bin_index = ((bin_index < 0) ? 0 : bin_index);
    pair<int, int> coordinate(segment_index, bin_index);
    segments_[segment_index][bin_index].AddPoint(range, point.z, coordinate);
    segment_coordinates_[i] = Bin::MinZPoint(range, point.z, coordinate);
    if (bin_index < segments_[segment_index].start_bin_index_) {
      segments_[segment_index].start_bin_index_ = bin_index;
    }
    if (bin_index > segments_[segment_index].end_bin_index_) {
      segments_[segment_index].end_bin_index_ = bin_index;
    }
  }
}

void GroundSegmenter::GetLines() {
  for (size_t i = 0; i < params_.n_segments; i++) {
    segments_[i].FitSegmentLines();
  }
}

void GroundSegmenter::AssignCluster(vector<uint8_t>& labels) {
  for (size_t i = 0; i < labels.size(); i++) {
    const Bin::MinZPoint& point_2d = segment_coordinates_[i];
    const int segment_index = point_2d.coordinate.first;
    if (segment_index >= 0) {
      double dist = segments_[segment_index].VerticalDistanceToLine(point_2d.d,
                                                                    point_2d.z);
      // Search neighboring segments.
      int steps = 1;
      while (dist == -1 && steps * segment_step_ < params_.line_search_angle) {
        // Fix indices that are out of bounds.
        int index_1 = segment_index + steps;
        while (index_1 >= params_.n_segments) index_1 -= params_.n_segments;
        int index_2 = segment_index - steps;
        while (index_2 < 0) index_2 += params_.n_segments;
        // Get distance to neighboring lines.
        const double dist_1 =
            segments_[index_1].VerticalDistanceToLine(point_2d.d, point_2d.z);
        const double dist_2 =
            segments_[index_2].VerticalDistanceToLine(point_2d.d, point_2d.z);
        // Select larger distance if both segments return a valid distance.
        if (dist_1 > dist) {
          dist = dist_1;
        }
        if (dist_2 > dist) {
          dist = dist_2;
        }
        ++steps;
      }
      if (dist < params_.max_dist_to_line && dist != -1) {
        labels.at(i) = 1;
      } else {
        labels.at(i) = 3;
      }
    }
  }
}
