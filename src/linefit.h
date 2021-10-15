#ifndef LINEFIT_H
#define LINEFIT_H

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <atomic>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <vector>

using namespace std;

struct GroundSegmentationParams {
  // Minimum range of labels.
  double r_min = 0.0;
  // Maximum range of labels.
  double r_max = 71.0;
  // Number of radial bins.
  int n_bins = 85;
  // Number of angular segments.
  int n_segments = 85;
  // Maximum distance to a ground line to be classified as ground.
  double max_dist_to_line = 0.15;
  // Max slope to be considered ground line.
  double max_slope = 0.15;
  // Max error for line fit.
  double max_error_square = 0.0015;
  // Distance at which points are considered far from each other.
  double long_threshold = 5.0;
  // Maximum slope for
  double max_long_height = 0.4;
  // Maximum heigh of starting line to be labelled ground.
  double max_start_height = 0.4;
  // The distance to extend outward when judging whether the point is online.
  double line_expand_range = 0.15;
  // How far to search for a line in angular direction [rad].
  double line_search_angle = 0.5;
  // Height of sensor above ground.
  double sensor_height = 0.43;
};

class Bin {
 public:
  struct MinZPoint {
    MinZPoint() : d(0.0), z(0.0), coordinate(-1, -1) {}
    MinZPoint(const double& d, const double& z,
              const pair<int, int>& coordinate)
        : d(d), z(z), coordinate(coordinate) {}

    double d;
    double z;
    pair<int, int> coordinate;
  };

 public:
  Bin();
  void AddPoint(const double& d, const double& z,
                const pair<int, int>& coordinate);
  MinZPoint GetMinZPoint();
  inline bool HasPoint() { return has_point_; }

 private:
  bool has_point_;
  double min_z_;
  double min_z_range_;
  pair<int, int> coordinate_;
};

class Segment {
 public:
  typedef pair<Bin::MinZPoint, Bin::MinZPoint> Line;
  typedef pair<double, double> LocalLine;

  int start_bin_index_;
  int end_bin_index_;

 public:
  Segment(const GroundSegmentationParams& params);

  void FitSegmentLines();

  double VerticalDistanceToLine(const double& d, const double& z);

  inline Bin& operator[](const size_t& index) { return bins_[index]; }

 private:
  GroundSegmentationParams params_;

  vector<Bin> bins_;

  list<Line> lines_;

  double GetMaxError(const list<Bin::MinZPoint>& points, const LocalLine& line);

  double GetMeanError(const list<Bin::MinZPoint>& points,
                      const LocalLine& line);

  LocalLine FitLocalLine(const list<Bin::MinZPoint>& points);

  Line LocalLineToLine(const LocalLine& local_line,
                       const list<Bin::MinZPoint>& line_points);
};

class GroundSegmenter {
 public:
  GroundSegmenter(
      const GroundSegmentationParams& params = GroundSegmentationParams());

  void Excute(const pcl::PointCloud<pcl::PointXYZ>& cloud,
              vector<uint8_t>& labels);

 private:
  GroundSegmentationParams params_;

  double segment_step_;
  double bin_step_;

  // Access with segments_[segment][bin].
  vector<Segment> segments_;

  // 2D coordinates (d, z) of every point in its respective segment.
  vector<Bin::MinZPoint> segment_coordinates_;

  void InsertPoints(const pcl::PointCloud<pcl::PointXYZ>& cloud,
                    vector<uint8_t>& labels);

  void GetLines();

  void AssignCluster(vector<uint8_t>& labels);
};

#endif  // LINEFIT_H
