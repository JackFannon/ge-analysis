#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "../all_headers.h"
#include "utilities.C"

const int nhists = 1;

int data_type = 0;

std::string source_name[nhists] = {
  "Co60_1",
};

int source_type[nhists] = {
  1,
};

std::string file_name[nhists] = {
  "../Data/2023/20230414_15MeV_Co60-01-00.csv",
};

int run_number[nhists] = {
    9999,
};

TH1D* hists[nhists];

TH1D* hists_calib[nhists];

void calibrate(){
    // Set gStyle options, see utilities.C for the function
    set_style(132);

    // Set max number of bins
    const int nbins = 4096;

    std::ofstream fit_results;
    fit_results.open("calib_fit_results.txt");
    gStyle->SetOptFit();
    gStyle->GetPadTopMargin();
    bool save_fit = true;
    bool K40norm = false;

    // Define vectors for energy and region of interest
    std::vector<double> true_energy[4];
    std::vector<double> roi_low[4];
    std::vector<double> roi_high[4];

    std::vector<double> true_energy_all;
    std::vector<double> ch_mean_all;
    std::vector<double> ch_error_all;
    std::vector<double> res_all;
    std::vector<double> res_error_all;

    std::string roi_filename;
    roi_filename = "roi.txt";

    std::string buffer;

    std::ifstream roi_data;
    roi_data.open(roi_filename.c_str());

    while (!roi_data.eof()) {
        std::getline(roi_data, buffer);
        std::stringstream roi_line;
        double e_true;
        double ch_low;
        double ch_high;
        int type = -1;
        roi_line >> e_true >> ch_low >> ch_high >> type;

        if(type >= 0){
            true_energy[type].push_back(e_true);
            roi_low[type].push_back(ch_low);
            roi_high[type].push_back(ch_high);
        }
    }
}
