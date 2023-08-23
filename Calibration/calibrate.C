#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMath.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include "utilities.C"
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

//====================================================================================================
//======================================== CONFIG OPTIONS ============================================
//====================================================================================================

const std::string DATA_LIST_DIR = "/Users/jack/Software/GeAnalysis/Configs/";
const std::string DATA_LIST = "newGe_list.txt";

// Set the number of bins for the detector
const int nbins = 4096;

// const std::string DATA_LIST = "old_nif_data.txt";
// const std::string DATA_LIST_DIR = "../Data/";
// const std::string DATA_DIR = DATA_LIST_DIR + "CrossCalibration/";

// const std::string ROI_LIST = "cross_calib_roi_list.txt";
// const std::string ROI_LIST_DIR = "./";

// const std::string OUTPUT_ROOT_FILE = "crosscalibold.root";

// const std::string ROI_FILE = "old_detector_roi.txt";

//====================================================================================================
//====================================================================================================
//====================================================================================================

struct region_of_interest {
    std::vector<double> true_energy;
    std::vector<int> lower_bound;
    std::vector<int> upper_bound;
    std::vector<int> isotope_type;
    std::vector<std::string> isotope_name;
};

region_of_interest read_roi_file(std::string roi_filename) {
    // Initalise a region_of_interest struct to store roi values in
    region_of_interest my_roi;

    // Read roi values from the file
    std::string buffer;
    std::ifstream roi_data;

    roi_data.open(roi_filename.c_str());

    if (!roi_data.is_open()) {
        std::cout << "!!!!! Cannot find " << roi_filename << "!!!!!" << std::endl;
    }

    while (!roi_data.eof()) {
        std::getline(roi_data, buffer);
        std::stringstream roi_line(buffer);
        double true_ene;
        double ch_low;
        double ch_high;
        int type;
        std::string name;
        roi_line >> true_ene >> ch_low >> ch_high >> type >> name;
        if (name.empty()) {
            continue;
        }
        my_roi.true_energy.push_back(true_ene);
        my_roi.upper_bound.push_back(ch_high);
        my_roi.lower_bound.push_back(ch_low);
        my_roi.isotope_type.push_back(type);
        my_roi.isotope_name.push_back(name);
    }
    return my_roi;
}

bool calibrate() {
    // Open output ROOT file

    // Histogram vectors
    //     raw counts vs channel number
    std::vector<TH1F *> hist_vect_raw_counts;

    // Loop over the files
    std::ifstream data_list;
    data_list.open((DATA_LIST_DIR + DATA_LIST).c_str());

    if (!data_list.is_open()) {
        std::cerr << "Could not open " << DATA_LIST_DIR + DATA_LIST << std::endl;
        return false;
    }

    std::string data_list_buffer;
    while (!data_list.eof()) {
        // Read the line with the filename from the file
        std::getline(data_list, data_list_buffer);
        std::stringstream data_list_line(data_list_buffer);
        std::string output_file_name;
        std::string data_file_name;
        std::string roi_file_name;
        data_list_line >> output_file_name >> data_file_name >> roi_file_name;
        TFile *output_file = new TFile((output_file_name + ".root").c_str(), "RECREATE");

        // Check that the line is populated before continuing
        if (data_file_name.empty()) {
            continue;
        }

        // Create and fill the histogram of raw counts vs channel
        TH1F *counts_channel_hist = new TH1F("", "", nbins, 0, nbins);
        read_data_into_hist(data_file_name, counts_channel_hist);

        // Get the contents of the region of interest (roi) file
        region_of_interest roi = read_roi_file(roi_file_name);

        // Create vectors to store the fit, the fit mean and the fit error
        std::vector<TF1 *> peak_vect;
        std::vector<double> peak_mean_vect;
        std::vector<double> peak_error_vect;

        // Fit all of the roi peaks
        for (int isotope = 0; isotope < roi.true_energy.size(); isotope++) {
            double peak_mean;
            double peak_error;
            TF1 *peak_fit =
                fit_peak_ge(counts_channel_hist, roi.lower_bound[isotope], roi.upper_bound[isotope], &peak_mean, &peak_error);
            peak_vect.push_back(peak_fit);
            peak_mean_vect.push_back(peak_mean);
            peak_error_vect.push_back(peak_error);

            // Create a canvas and draw the histogram and fit
            TCanvas *isotope_canvas = new TCanvas(roi.isotope_name[isotope].c_str(), "", 600, 600);

            isotope_canvas->cd();

            counts_channel_hist->SetName(roi.isotope_name[isotope].c_str());
            counts_channel_hist->SetAxisRange(0.9 * roi.lower_bound[isotope], 1.1 * roi.upper_bound[isotope]);
            counts_channel_hist->Draw();
            peak_fit->Draw("same");
            isotope_canvas->Write();
            output_file->Write();
            delete isotope_canvas;
        }

        // Plot the calibrated energy graph
        TCanvas *calibration_canvas = new TCanvas("Calibration", "", 1200, 600);
        calibration_canvas->Divide(1, 2);

        // Setup pad for the linear regression plot
        TVirtualPad *linear_pad = calibration_canvas->cd(1);
        TH1F *linear_hist = linear_pad->DrawFrame(0, 0, 10000, nbins, ";True energy (keV); MCA channel");
        linear_hist->SetTitleSize(0.06, "XY");
        linear_hist->SetTitleSize(0.06, "XY");
        linear_hist->SetTitleOffset(0.7, "Y");
        linear_pad->SetLeftMargin(0.15);
        linear_pad->SetRightMargin(0.15);
        linear_pad->SetBottomMargin(0.15);
        linear_pad->SetTopMargin(0.);

        // Plot a graph of true energy against channel number
        TGraphErrors *energy_channel_graph = new TGraphErrors();
        for (int isotope = 0; isotope < roi.true_energy.size(); isotope++) {
            energy_channel_graph->SetPoint(isotope, roi.true_energy[isotope], peak_mean_vect[isotope]);
            energy_channel_graph->SetPointError(isotope, 0., peak_error_vect[isotope]);
        }
        energy_channel_graph->SetTitle("Calibrated energy");
        energy_channel_graph->Draw();
        energy_channel_graph->SetMarkerStyle(20);
        energy_channel_graph->Draw("PSAME");

        // Fit a linear fit to the graph
        energy_channel_graph->Fit("pol1", "", "", 0, 10000);

        // Get the fit parameters
        TF1 *linear_function = energy_channel_graph->GetFunction("pol1");
        double p0 = linear_function->GetParameter(0);
        double p1 = linear_function->GetParameter(1);
        double intercept = -p0 / p1;
        double slope = 1.0 / p1;

        // Setup pad for the residuals graph
        TVirtualPad *residuals_pad = calibration_canvas->cd(2);
        residuals_pad->SetLeftMargin(0.15);
        residuals_pad->SetRightMargin(0.15);
        residuals_pad->SetBottomMargin(0.15);
        residuals_pad->SetTopMargin(0.);

        // Calculate the residual energy
        std::vector<double> residual_energy;
        std::vector<double> residual_error;
        for (int isotope = 0; isotope < roi.true_energy.size(); isotope++) {
            double expected_channel = p0 + p1 * roi.true_energy[isotope];
            residual_energy.push_back((peak_mean_vect[isotope] - expected_channel) / expected_channel);
            residual_error.push_back(peak_error_vect[isotope] / expected_channel);
        }

        // Setup and fill TGraphErrors for residuals
        TGraphErrors *residual_graph = new TGraphErrors();
        double max_res = -9999.;
        double min_res = 9999.;

        for (int isotope = 0; isotope < roi.true_energy.size(); isotope++) {
            residual_graph->SetPoint(isotope, roi.true_energy[isotope], residual_energy[isotope]);
            residual_graph->SetPointError(isotope, 0., residual_error[isotope]);
            if (residual_energy[isotope] + residual_error[isotope] > max_res) {
                max_res = residual_energy[isotope] + residual_error[isotope];
            }
            if (residual_energy[isotope] - residual_error[isotope] < min_res) {
                min_res = residual_energy[isotope] - residual_error[isotope];
            }
        }

        // Setup pad for the residual frame
        TH1F *residuals_hist = residuals_pad->DrawFrame(0.0, 1.1 * min_res, 10000.0, 1.1 * max_res);
        residuals_pad->SetGrid();
        residuals_hist->SetTitleSize(0.06, "XY");
        residuals_hist->SetYTitle("(E_{obs} - E_{fit})/E_{fit}");
        residuals_hist->SetXTitle("True energy (keV)");
        residuals_hist->SetTitleOffset(0.7, "Y");

        residual_graph->SetMarkerStyle(20);
        residual_graph->Draw("PSAME");

        peak_vect.clear();
        peak_mean_vect.clear();
        peak_error_vect.clear();
        residual_energy.clear();
        residual_error.clear();

        // Print out slope and intercept values
        std::cout << "==========================================" << std::endl;
        std::cout << "Slope = " << slope << std::endl;
        std::cout << "Intercept = " << intercept << std::endl;
        std::cout << "==========================================" << std::endl;

        calibration_canvas->Write();
        output_file->Close();
    }

    return true;
}
