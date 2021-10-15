#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <fstream>
#include <iostream>
#include <memory>

#include "linefit.h"

using namespace std;

void ReadPoints(const string&, pcl::PointCloud<pcl::PointXYZ>&, const int& = 6);
void AreaPickingEventOccurred(const pcl::visualization::AreaPickingEvent&,
                              void*);
void save_bin(const string& path, const string& out_path, const vector<int> label);
pcl::PointCloud<pcl::PointXYZ>::Ptr points_ptr(
    new pcl::PointCloud<pcl::PointXYZ>);

int main() {
  string point_path = "/home/fzx/code/center-point/data/kitti/training/velodyne/";
  string outPpath = "/home/fzx/code/center-point/data/kitti/training/velodyne_no_gp/";
  int total_bin_num = 13976;
  for (int num = 0; num < total_bin_num; num++)
    { 
      char fl[256];
      sprintf(fl,"%06d.bin", num);
      string cur_point_path = point_path+string(fl);
      ReadPoints(cur_point_path, *points_ptr);
      cout << "Point size: " << points_ptr->size() << endl;

      shared_ptr<GroundSegmenter> gs_ptr(new GroundSegmenter());
      vector<uint8_t> labels;
      gs_ptr->Excute(*points_ptr, labels);
      pcl::PointCloud<pcl::PointXYZ>::Ptr gnd_points_ptr(
          new pcl::PointCloud<pcl::PointXYZ>);
      gnd_points_ptr->reserve(labels.size());
      pcl::PointCloud<pcl::PointXYZ>::Ptr obs_points_ptr(
          new pcl::PointCloud<pcl::PointXYZ>);
      obs_points_ptr->reserve(labels.size());
      vector<int>objects_label;
      for (int i = 0; i < labels.size(); i++)
      {
        if (labels[i] == 1)
        {
          gnd_points_ptr->points.push_back(points_ptr->at(i));
        }
        else if (labels[i] == 3)
        {
          objects_label.push_back(i);
          obs_points_ptr->points.push_back(points_ptr->at(i));
        }
      }
      string out_file_path = outPpath + string(fl);
      save_bin(cur_point_path, out_file_path, objects_label);
    }

  // cout << "Gnd size: " << gnd_points_ptr->size()
  //      << " Obs size: " << obs_points_ptr->size()
  //      << " Total size: " << gnd_points_ptr->size() + obs_points_ptr->size()
  //      << endl;

  // pcl::visualization::PCLVisualizer::Ptr viewer_ptr(
  //     new pcl::visualization::PCLVisualizer("3D"));
  // viewer_ptr->addCoordinateSystem();
  // viewer_ptr->registerAreaPickingCallback(AreaPickingEventOccurred,
  //                                         (void*)&(points_ptr));

  // pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> gnd_color(
  //     gnd_points_ptr, 0, 255, 0);
  // viewer_ptr->addPointCloud(gnd_points_ptr, gnd_color, "gnd");
  // pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> obs_color(
  //     obs_points_ptr, 255, 0, 0);
  // viewer_ptr->addPointCloud(obs_points_ptr, obs_color, "obs");

  // while (!viewer_ptr->wasStopped()) {
  //   viewer_ptr->spinOnce(100);
  // }

  return 0;
}

void save_bin(const string& in_path,const string& out_path,  const vector<int> label)
{
  ifstream fin(in_path.c_str(), ios::binary);

  ofstream wf(out_path.c_str(), ios::binary);

  int num_feature = 6;
  fin.seekg(0, ios::end);
  size_t num_points = fin.tellg() / (num_feature * sizeof(float));
  fin.seekg(0, ios::beg);
  vector<float> values(num_feature * num_points);
  fin.read((char*)&values[0], num_feature * num_points * sizeof(float));
  fin.close();
  for (size_t ni = 0; ni < label.size(); ni++) {
    size_t index = num_feature * label[ni];
    wf.write((char*)&values[index], num_feature*sizeof(float));
  }
  wf.close();
}

void ReadPoints(const string& path, pcl::PointCloud<pcl::PointXYZ>& clout_out,
                const int& size) {
  clout_out.clear();
  ifstream fin(path.c_str(), ios::binary);
  if (!fin.is_open()) {
    cerr << "[ReadPoints]: can't open path: " << path << endl;
    return;
  }
  fin.seekg(0, ios::end);
  size_t num_points = fin.tellg() / (size * sizeof(float));
  fin.seekg(0, ios::beg);
  vector<float> values(size * num_points);
  fin.read((char*)&values[0], size * num_points * sizeof(float));
  fin.close();
  clout_out.resize(num_points);
  for (size_t ni = 0; ni < num_points; ni++) {
    auto& point = clout_out[ni];
    size_t index = size * ni;
    point.x = values[index];
    point.y = values[index + 1];
    point.z = values[index + 2];
  }
}

void AreaPickingEventOccurred(const pcl::visualization::AreaPickingEvent& event,
                              void* args) {
  vector<int> indices;
  if (event.getPointsIndices(indices) == -1) return;
  for (size_t i = 0; i < indices.size(); i++) {
    const auto& pt = points_ptr->at(indices[i]);
    cout << pt.x << " " << pt.y << " " << pt.z << endl;
  }
}
