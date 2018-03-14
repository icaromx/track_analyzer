#include "TROOT.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TError.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TStyle.h"
#include "TString.h"
#include "TVector3.h"
#include "TCanvas.h"
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include "/grid/fermiapp/products/larsoft/eigen/v3_3_3/include/eigen3/Eigen/Dense"

//#include "/usr/local/Cellar/eigen/3.3.4/include/eigen3/Eigen/Dense" //Needed on MACOS
using namespace std;

struct Point {
  float x;
  float y;
  float z;
  float q;
};
struct PCAResults {
  TVector3 centroid;
  pair<TVector3,TVector3> endPoints;
  float length;
  TVector3 eVals;
  vector<TVector3> eVecs;
};
struct TrkPoint{
    double c;
    double x;
    double y;
    double z;
    double q;
};
struct by_y { 
    bool operator()(TrkPoint const &a, TrkPoint const &b) { 
        if(a.y == b.y) return a.x > b.x;
        else return a.y > b.y;
    }
};
struct reverse_by_y { 
    bool operator()(TrkPoint const &a, TrkPoint const &b) { 
        if(a.y == b.y) return a.x < b.x;
        else return a.y < b.y;
    }
};
typedef vector<TrkPoint> track_def;
typedef vector<Point> PointCloud;
void LoadPointCloud(PointCloud &points, const track_def &ord_trk);
PCAResults DoPCA(const PointCloud &points);


/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////MAIN PROGRAM STARTS////////////////////////////////////////

int main(int argc, char **argv){
  std::string line;
  TFile *f_output;
  //char *inputs = argv[1];
  std::ifstream ifs(argv[1]);
  //std::ifstream ifs("../event_762.txt");
  //std::ifstream ifs("event_5253.txt");
  //std::ifstream ifs("../run_sample_filelist.txt");
  f_output = TFile::Open("../results.root","RECREATE");
  
  TH1F *hist_angles = new TH1F("hist_angles","hist_angles;Cos(#theta);Events",100,-1,1);
  TH1F *hist_low_cos = new TH1F("hist_low_cos","Cos(#theta) < 0.8;Cos(#theta);Events",100,-1,1);
  TNtuple *nt_tracks = new TNtuple("tracks","tracks","run:event:c:x:y:z:q");

  system("mkdir ../csv_files/");
  system("mkdir ../track_plots/");


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////OUTPUT COMMENTS//////////////////////////////////////////////////
  cout << "Using only topological cuts with this selection." << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector< std::vector<int> > kept_tracks;

  while(std::getline(ifs, line)){
    gROOT->Reset();
    gErrorIgnoreLevel = kError;

    TString filename;
    filename.Form("%s",line.c_str());    
    TFile *infile = new TFile(filename);

    ////////////////PARAMETERS//////////////////////
    double pe_threshold = 10;
    double vdrift = 0.1101; //cm/us CONSTANT
    double x_anode = 0.0;
    double x_cathode = 256.4;
    double x_anode_th = 3.0;
    double x_cathode_th = 4.0; 

    int res = 20;
   // int track_num = 14;
    double alpha = 6.;

    //Extract Event Metadata
    TTree *Trun = (TTree*)infile->Get("Trun");
    Int_t run_num;
    Int_t ev_num;
    Trun->SetBranchAddress("runNo",&run_num);
    Trun->SetBranchAddress("eventNo",&ev_num);
    Trun->GetEntry(0);
    //cout << "Looking at run " << run_num << " from event " << ev_num << endl;

    //Extract Coordinate information
    TTree *T_charge_cluster = (TTree*)infile->Get("T_charge_cluster_nfc"); 
    Double_t cluster_id;
    Double_t qx;
    Double_t qy;
    Double_t qz;
    Double_t qc;
    
    T_charge_cluster->SetBranchAddress("qx",&qx);
    T_charge_cluster->SetBranchStatus("qx", kTRUE);
    T_charge_cluster->SetBranchAddress("qy",&qy);
    T_charge_cluster->SetBranchStatus("qy", kTRUE);
    T_charge_cluster->SetBranchAddress("qz",&qz);
    T_charge_cluster->SetBranchStatus("qz", kTRUE);
    T_charge_cluster->SetBranchAddress("qc",&qc);
    T_charge_cluster->SetBranchStatus("qc", kTRUE);
    T_charge_cluster->SetBranchAddress("cluster_id", &cluster_id);
    T_charge_cluster->SetBranchStatus("cluster_id", kTRUE);
    int size = T_charge_cluster->GetEntries();
    /////////////////////////////////////////////////////////////
    //Extract Clusters///////////////////////////////////////////
    std::vector<Int_t> clusters;
    Int_t prev_cval;

    for (int i = 0; i < size; ++i){
      T_charge_cluster -> GetEntry(i);
      if (i == 0){
        clusters.push_back(cluster_id);
        prev_cval = cluster_id;
      }else if(prev_cval != cluster_id){
        prev_cval = cluster_id;
        clusters.push_back(cluster_id);
      }
    }

    //Looking at tracks individually
    int num_clusters = clusters.size();

    //Loop through individual clusters
    std::vector<int> event_kept_trks;
    for (int c = 0; c < num_clusters; ++c){
      	int cluster = clusters[c];

      	//if(cluster != 14) continue;
      	track_def trk;
      	//Load every point (cluster id, x, y, z, charge) of a track into the trk object.
      	for (int i = 0; i < size; ++i){
	        T_charge_cluster -> GetEntry(i);
	        //Will only store information for current cluster
	        if(cluster_id != cluster) continue;
	        //if(cluster != 14) continue;
	        
	        TrkPoint tempPoint;
	        tempPoint.c = cluster_id;
	        tempPoint.x = qx;
	        tempPoint.y = qy;
	        tempPoint.z = qz;
	        tempPoint.q = qc;
	        trk.push_back(tempPoint);
	      }
	      //track size has to be larger than the moving window size
	      if(trk.size() < res) continue;
	      

	      //////////////////////////////////////////////////
	      ////////ORDERING ALGORITHM BEGINS/////////////////
	      //Sort track in descending y value
	      std::sort(trk.begin(), trk.end(), by_y());
	      int trk_size = trk.size();
	      track_def points_left;
	      track_def points_gd;

	      for (int i = 1; i < trk_size; ++i){
	        TrkPoint tempPoint;
	        tempPoint.c = trk[i].c;
	        tempPoint.x = trk[i].x;
	        tempPoint.y = trk[i].y;
	        tempPoint.z = trk[i].z;
	        tempPoint.q = trk[i].q;
	        points_left.push_back(tempPoint);
	      }

	      track_def ord_trk;
	      ord_trk.push_back(trk[0]);
	      int pl_size = points_left.size();
	      double old_dist = 10000000;
	      int low_dist_at = -1;
	      double dist;
	      std::vector<int> used;
	      std::vector<int> recheck;
	      int m = 0;
	      int i = 0;
	      used.clear();
	      recheck.clear();
	      double low_ord_y = 2000.;
	      double low_x = 300.;
	      double high_x = -300.;

	    while(pl_size != 0){
	        for (int j = 0; j < pl_size; ++j){
	          dist = sqrt(pow(ord_trk.back().x - points_left[j].x, 2.0) + pow(ord_trk.back().y - points_left[j].y, 2.0) + pow(ord_trk.back().z - points_left[j].z, 2.0));
	          if (dist < old_dist){
	            old_dist = dist;
	            low_dist_at = j;
	          }
	        }
	        TrkPoint tempPoint;
	        tempPoint.c = points_left[low_dist_at].c;
	        tempPoint.x = points_left[low_dist_at].x;
	        tempPoint.y = points_left[low_dist_at].y;
	        tempPoint.z = points_left[low_dist_at].z;
	        tempPoint.q = points_left[low_dist_at].q;
	        if (old_dist > alpha){
	          points_gd.push_back(tempPoint);
	          old_dist = 10000000;
	          points_left.erase (points_left.begin() + low_dist_at);
	          pl_size = points_left.size();
	          i++;
	        }else{
	          ord_trk.push_back(tempPoint);
	          if (tempPoint.y < low_ord_y){
	          	low_ord_y = tempPoint.y;
	          }
	          if (tempPoint.x < low_x){
	          	low_x = tempPoint.x;
	          }
	          if (tempPoint.x > high_x){
	          	high_x = tempPoint.x;
	          }
	          old_dist = 10000000;
	          points_left.erase (points_left.begin() + low_dist_at);
	          pl_size = points_left.size();
	          i = 0;
	        }
	        if (pl_size == 0) break;
	    }

		///////////////////////////////////////////////////////////////////////
		///////////////Finished Ordering Points////////////////////////////////
		double bottom_dist;

		bottom_dist = abs(trk.back().y - low_ord_y);

		if (bottom_dist > 10) continue;


		////////////////////////////////////////////////////////////////////////
		////////////////////////MAKING CSV FILES////////////////////////////////
	    ofstream outfile;
		std::string cluster_string = std::to_string(cluster);
		std::string unique_string = std::to_string(ev_num);
		//string csv_filename = "track_" + std::to_string(cn) + "_" + cluster_string + ".csv";
		string csv_filename = "../csv_files/track_" + cluster_string + "_" + unique_string + ".csv";
		outfile.open(csv_filename);
		for (int i = 0; i < ord_trk.size(); ++i){
			outfile << ord_trk.at(i).x << ", " << ord_trk.at(i).y << endl;
			nt_tracks -> Fill(run_num, ev_num,ord_trk.at(i).c, ord_trk.at(i).x, ord_trk.at(i).y, ord_trk.at(i).z, ord_trk.at(i).q);
		} 
		outfile.close();
	    ///////////////////////////////////////////////////////////////////////
		///////////////Starting Moving Window//////////////////////////////////
	 	
 		int win_size;
 		//for (int w = 6; w < 50; ++w){ //<window iter>
 		//win_size = w;
	 	track_def prev_chunk;
	 	track_def post_chunk;
	 	PointCloud prev_points;
	 	PointCloud post_points;
	 	TrkPoint VertexPoint;

	 	double dotProd, min_ang = 1.;
	 	double ev_lowest = 100000;
	 	double rr = 0.0;
	 	win_size = 40;
	 	int buffer_size = 0;
	 	int effective_window = win_size + buffer_size;
	 	int vertex;
	    if(ord_trk.size() < win_size*2) continue;
	    for (int i = effective_window; i < ord_trk.size() - effective_window; ++i){
			track_def prev_chunk;
		 	track_def post_chunk;
		 	PointCloud prev_points;
		 	PointCloud post_points;    
		 	PCAResults prev_results;
		 	PCAResults post_results;

	     	for (int j = i - effective_window; j < i - buffer_size; ++j){
	     		TrkPoint prev_tempPoint;
		        prev_tempPoint.c = ord_trk[j].c;
		        prev_tempPoint.x = ord_trk[j].x;
		        prev_tempPoint.y = ord_trk[j].y;
		        prev_tempPoint.z = ord_trk[j].z;
		        prev_tempPoint.q = ord_trk[j].q;
		        prev_chunk.push_back(prev_tempPoint);	
	     	}
	     	for (int j = i + buffer_size; j < i + effective_window; ++j){
	     		TrkPoint post_tempPoint;
		        post_tempPoint.c = ord_trk[j].c;
		        post_tempPoint.x = ord_trk[j].x;
		        post_tempPoint.y = ord_trk[j].y;
		        post_tempPoint.z = ord_trk[j].z;
		        post_tempPoint.q = ord_trk[j].q;
		        post_chunk.push_back(post_tempPoint);	
	     	}
	     	LoadPointCloud(prev_points, prev_chunk);
	     	LoadPointCloud(post_points, post_chunk);
	     	prev_results = DoPCA(prev_points);
	     	post_results = DoPCA(post_points);
	     	//cout << prev_results.length << ", " << post_results.length << endl;
	     	rr += prev_results.length/win_size;
	     	dotProd = prev_results.eVecs[0](0)*post_results.eVecs[0](0) + prev_results.eVecs[0](1)*post_results.eVecs[0](1) + prev_results.eVecs[0](2)*post_results.eVecs[0](2);
	     	//cout << rr << ", " << dotProd << endl;
	     	

	     	if(dotProd < abs(min_ang)){
	     		vertex = i + 0;
	     		min_ang = dotProd;
	     		VertexPoint.c = ord_trk[vertex].c;
	     		VertexPoint.x = ord_trk[vertex].x;
	     		VertexPoint.y = ord_trk[vertex].y;
	     		VertexPoint.z = ord_trk[vertex].z;
	     		VertexPoint.q = ord_trk[vertex].c;
	     	}



	    }
	   	//////////////////////////////////////////////////////////////////////
	    double dist_low_vert_y;
	    dist_low_vert_y = abs(VertexPoint.y - low_ord_y);
	    hist_angles -> Fill(min_ang);
	    if(dist_low_vert_y > 10) continue;

	    //if(min_ang < 0.8) cout << ev_num << ", " << run_num << ", " << cluster << ", " << VertexPoint.x << ", " << VertexPoint.y << ", " << min_ang << ", " << dist_low_vert_y << endl;
	    if(min_ang < 0.8) cout << Form("Run = %d, Event = %d, Track = %d", run_num, ev_num, cluster) << endl;

	    hist_angles -> Fill(min_ang);
	    if(min_ang < 0.8) hist_low_cos->Fill(min_ang);
	    if(min_ang < 0.8){
	    	event_kept_trks.push_back(run_num);
	    	event_kept_trks.push_back(ev_num);
	    	event_kept_trks.push_back(cluster);
	    	kept_tracks.push_back(event_kept_trks);
	    }

	 // } <window iter>
   }
   infile->Close();
  }
  
  for (int i = 0; i < kept_tracks.size(); ++i){
	TCanvas c1;
	nt_tracks->Draw("y:x",Form("run == %d && event == %d && c == %d", kept_tracks[i][0], kept_tracks[i][1], kept_tracks[i][2]), "*");
	c1.SaveAs(Form("../track_plots/Run_%d_Ev_%d_Trk_%d.pdf", kept_tracks[i][0], kept_tracks[i][1], kept_tracks[i][2]));
	//c1.SaveAs(Form("track_plots/Run_%d_Ev_%d_Trk_%d.png", kept_tracks[i][0], kept_tracks[i][1], kept_tracks[i][2]));
  }
  f_output->Write();
  f_output->Close();
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////END OF MAIN PROGRAM/////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////START OF FUNCTIONS//////////////////////////////////////

void LoadPointCloud(PointCloud &points, const track_def &ord_trk) {
  for (int i = 0; i < ord_trk.size(); ++i){
    Point tempPoint;
    tempPoint.x = ord_trk.at(i).x;
    tempPoint.y = ord_trk.at(i).y;
    tempPoint.z = ord_trk.at(i).z;
    tempPoint.q = ord_trk.at(i).q;
    
    points.push_back(tempPoint);

  }
  return;
}

PCAResults DoPCA(const PointCloud &points) {
  TVector3 outputCentroid;
  pair<TVector3,TVector3> outputEndPoints;
  float outputLength;
  TVector3 outputEigenValues;
  vector<TVector3> outputEigenVecs;
  float meanPosition[3] = {0., 0., 0.};
  unsigned int nThreeDHits = 0;
  for (unsigned int i = 0; i < points.size(); i++) {
    meanPosition[0] += points[i].x;
    meanPosition[1] += points[i].y;
    meanPosition[2] += points[i].z;
    ++nThreeDHits;
  }
  if (nThreeDHits == 0) {
    PCAResults results;
    return results; 
  }
  const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
  meanPosition[0] /= nThreeDHitsAsFloat;
  meanPosition[1] /= nThreeDHitsAsFloat;
  meanPosition[2] /= nThreeDHitsAsFloat;
  outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);
  float xi2 = 0.0;
  float xiyi = 0.0;
  float xizi = 0.0;
  float yi2 = 0.0;
  float yizi = 0.0;
  float zi2 = 0.0;
  float weightSum = 0.0;
  for (unsigned int i = 0; i < points.size(); i++) {
      const float weight(1.);
      const float x((points[i].x - meanPosition[0]) * weight);
      const float y((points[i].y - meanPosition[1]) * weight);
      const float z((points[i].z - meanPosition[2]) * weight);
      xi2  += x * x;
      xiyi += x * y;
      xizi += x * z;
      yi2  += y * y;
      yizi += y * z;
      zi2  += z * z;
      weightSum += weight * weight;
  }

  Eigen::Matrix3f sig;

  sig << xi2, xiyi, xizi,
         xiyi, yi2, yizi,
         xizi, yizi, zi2;

  sig *= 1.0 / weightSum;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

  typedef std::pair<float,size_t> EigenValColPair;
  typedef std::vector<EigenValColPair> EigenValColVector;

  EigenValColVector eigenValColVector;
  const auto &resultEigenMat(eigenMat.eigenvalues());
  eigenValColVector.emplace_back(resultEigenMat(0), 0);
  eigenValColVector.emplace_back(resultEigenMat(1), 1);
  eigenValColVector.emplace_back(resultEigenMat(2), 2);

  std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;} );

  outputEigenValues = TVector3(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

  const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

  for (const EigenValColPair &pair : eigenValColVector) {
     outputEigenVecs.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
  }

  PCAResults results;

  Eigen::ParametrizedLine<float,3> priAxis(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)),Eigen::Vector3f(outputEigenVecs[0](0),outputEigenVecs[0](1),outputEigenVecs[0](2)));

  Eigen::Vector3f endPoint1(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));
  Eigen::Vector3f endPoint2(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));

  Eigen::Vector3f testPoint;
  Eigen::Vector3f projTestPoint;
  float maxDist1 = -1.0;
  float maxDist2 = -1.0;
  float dist;
  float dotP;
  for (unsigned int i = 0; i < points.size(); i++) {
    testPoint = Eigen::Vector3f(points[i].x,points[i].y,points[i].z);
    projTestPoint = priAxis.projection(testPoint);
    dist = sqrt(pow(projTestPoint(0)-outputCentroid(0),2.0)+pow(projTestPoint(1)-outputCentroid(1),2.0)+pow(projTestPoint(2)-outputCentroid(2),2.0));
    dotP = (projTestPoint(0)-outputCentroid(0))*outputEigenVecs[0](0) + (projTestPoint(1)-outputCentroid(1))*outputEigenVecs[0](1) + (projTestPoint(2)-outputCentroid(2))*outputEigenVecs[0](2);

    if ((dotP < 0.0) && (dist > maxDist1)) {
      endPoint1 = projTestPoint;
      maxDist1 = dist;
    }
    else if ((dotP > 0.0) && (dist > maxDist2)) {
      endPoint2 = projTestPoint;
      maxDist2 = dist;
    }
  }
  outputEndPoints.first = TVector3(endPoint1(0),endPoint1(1),endPoint1(2));
  outputEndPoints.second = TVector3(endPoint2(0),endPoint2(1),endPoint2(2));
  outputLength = sqrt(pow(endPoint2(0)-endPoint1(0),2.0)+pow(endPoint2(1)-endPoint1(1),2.0)+pow(endPoint2(2)-endPoint1(2),2.0));
  results.centroid = outputCentroid;
  results.endPoints = outputEndPoints;
  results.length = outputLength;
  results.eVals = outputEigenValues;
  results.eVecs = outputEigenVecs;
  return results;
}

