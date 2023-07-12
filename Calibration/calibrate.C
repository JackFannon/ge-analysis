#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include "utilities.C"

//====================================================================================================
//======================================== CONFIG OPTIONS ============================================
//====================================================================================================

const std::string DATA_LIST = "newGe_list.txt";
const std::string DATA_LIST_DIR = "../Data/2023/";
const std::string DATA_DIR = DATA_LIST_DIR + "raw/";
const std::string ROI_FILE = "roi.txt";

const std::string OUTPUT_ROOT_FILE = "histograms.root";

const int NUM_OF_ISOTOPES = 4;
const std::string ISOTOPE_SYMBOLS[NUM_OF_ISOTOPES] = {
"Bg(K40_Tl208)",
"Co",
"Cs",
"Ni"
};

bool Ni_sum_flag = false;
const int nhists = 1;
int run_number[nhists] = {
    9999,
};

// Set the number of bins for the detector
const int nbins = 4096;

const bool save_fit = true;

const bool K40norm = false;

//====================================================================================================
//====================================================================================================
//====================================================================================================

void calibrate(){
    // Set gStyle options, see utilities.C for the function
    set_style(132);

    // Read in the list of files from data_list.txt
    std::vector<std::string> ge_data_files = load_data(DATA_LIST, DATA_LIST_DIR);

    // Open the ROOT file to store output histograms in
    TFile* file = new TFile(OUTPUT_ROOT_FILE.c_str() ,"RECREATE");

    /* for(std::string file_name: ge_data_files){ */
    /*     plot_channel_hist(file_name, DATA_DIR)->Write(); */
    /* } */

    std::vector<std::vector<int> > isotopes_in_file;

    // Search for the isotope names in the filenames. If the isotope name is present then the data
    // should contain a peak for that isotope. Otherwise it's just background (0) and should contain peaks
    // from K40 and Tl208.
    // The source types are as follows:
    //     0 -- Background (K40 & Tl208)
    //     1 -- Co60
    //     2 -- Cs137
    //     3 -- NiCf
    for(int file_index = 0; file_index < ge_data_files.size(); file_index++){
        std::vector<int> isotopes;
        isotopes.push_back(0);

        for(int source_index = 0; source_index < NUM_OF_ISOTOPES; source_index++){
            if(ge_data_files[file_index].find(ISOTOPE_SYMBOLS[source_index]) != std::string::npos){
                isotopes.push_back(source_index);
            }
        }
        isotopes_in_file.push_back(isotopes);
    }

    // Print out the isotope type for each file
    int index = 0;
    /* for(std::vector<int> j: isotopes_in_file){ */

    /*     std::cout << "isotope type in : " << ge_data_files[index]; */
    /*     index++; */
    /*     for(int i: j){ */
    /*         std::cout << i << ", "; */
    /*     } */
    /*     std::cout << std::endl; */
    /* } */

    std::ofstream fit_results;
    fit_results.open("calib_fit_results.txt");


    //Define vectors for the:
    //    true_energy_all - the true energy of each isotope peak
    //    ch_mean_all     - the mean channel number for each peak
    //    ch_error_all    - the error on the mean channel number for each peak
    //    res_all         - the residual (E_expected - E_fit)/E_fit for each peak
    //    res_error_all   - the error on the above value
    std::vector<double> true_energy_all;
    std::vector<double> ch_mean_all;
    std::vector<double> ch_error_all;
    std::vector<double> res_all;
    std::vector<double> res_error_all;

//==========================================================================================
//=========================== EXTRACT REGION OF INTEREST DATA ==============================
//==========================================================================================
    // Define vectors for energy and region of interest
    std::vector<double> true_energy[4];
    std::vector<double> roi_low[4];
    std::vector<double> roi_high[4];

    std::string buffer;
    std::ifstream roi_data;

    roi_data.open(ROI_FILE.c_str());

    if(!roi_data.is_open()){
        std::cout << "!!!!! Cannot find " << ROI_FILE.c_str() << "!!!!!" << std::endl;
        return;
    }

    while (!roi_data.eof()) {
        std::getline(roi_data, buffer);
        std::stringstream roi_line(buffer);
        double e_true;
        double ch_low;
        double ch_high;
        int type = -1;
        // Extract the region of interest values
        roi_line >> e_true >> ch_low >> ch_high >> type;
        if(type >= 0){
            true_energy[type].push_back(e_true);
            roi_low[type].push_back(ch_low);
            roi_high[type].push_back(ch_high);
        }
    }
//==========================================================================================
//==========================================================================================
//==========================================================================================

    double mean;
    double error;
    double K40_peak_ref = -1;
    TH1D* hists[ge_data_files.size()];

    TH1D* Ni_sum = new TH1D("Co60+Ni_sum", "Co60+Ni_sum", nbins, 0, nbins);

    // Loop over the data files.
    for(int file_index = 0; file_index < ge_data_files.size(); file_index++){

        hists[file_index] = new TH1D(("h_" + ge_data_files[file_index]).c_str(), (ge_data_files[file_index] + ";Channel;Count").c_str(), nbins, 0, nbins);

        // Read data out of the file ge_data_files[file_index] into the histogram created above
        read_data_into_hist(DATA_DIR + ge_data_files[file_index], hists[file_index]);

        //============================================================================================
        //==================================== FIT THE PEAKS =========================================
        //============================================================================================
        // Loop over isotope peaks that should be present in the data for this file. Using two for loops here as there are three types of data:
        // Background only -- K40 and Tl208 peaks (2 peaks)
        // Co60            -- Two C60 peaks and the background peaks (at least 4 peaks)
        // NiCf            -- Multiple NiCf peaks and the background peaks (at least 15 peaks)
        // Got to loop over the data type and then loop over the peaks that are present in that data type
        for (int isotope_type: isotopes_in_file[file_index]){
            for(int isotope_peak = 0; isotope_peak < roi_high[isotope_type].size(); isotope_peak++){

                // Setup a canvas named after the filename, "file type (BG, Co, Ni)" and the peak number
                std::string canvas_name = ge_data_files[file_index] + ISOTOPE_SYMBOLS[isotope_type] + std::to_string(isotope_peak);

                std::cout << ge_data_files[file_index] << "    " << canvas_name << std::endl;


                TCanvas* my_canvas = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 600, 600);
                my_canvas->SetTitle(ge_data_files[file_index].c_str());

                // Fit the peak with a gaussian using the information from the region of interest file
                fit_peak_ge(hists[file_index], roi_low[isotope_type][isotope_peak], roi_high[isotope_type][isotope_peak], &mean, &error);

                // Draw the histogram
                hists[file_index]->Draw();

                // Store information about the fit (mean and error) and the true energy that the peak should represent
                true_energy_all.push_back(true_energy[isotope_type][isotope_peak]);
                ch_mean_all.push_back(mean);
                ch_error_all.push_back(error);

                // Write the histogram to the root file with a slightly larger axis range than the range of interest.
                if(!save_fit){
                    hists[file_index]->SetAxisRange(roi_low[isotope_type][isotope_peak] - 50, roi_high[isotope_type][isotope_peak] + 50);
                    my_canvas->Write();
                }
                if(save_fit){
                    hists[file_index]->SetAxisRange(roi_low[isotope_type][isotope_peak] - 50, roi_high[isotope_type][isotope_peak] + 50);
                    my_canvas->SaveAs(("../Output/" + canvas_name + ".png").c_str());
                }
            }
        }
        //============================================================================================
        //============================================================================================
        //============================================================================================

        //--------------------------------------------------------------------------------------------

        //============================================================================================
        //============================= PLOT CALIBRATED ENERGY GRAPH =================================
        //============================================================================================
        std::string print_name = "calib" + ge_data_files[file_index];
        TCanvas *c_calibfit = new TCanvas(print_name.c_str() , print_name.c_str(), 600, 600);
        c_calibfit->Divide(1, 2);
        TVirtualPad *pad1 = c_calibfit->cd(1);
        TH1F *h1 = pad1->DrawFrame(0, 0, 10000, nbins, ";True energy (keV); MCA channel");
        pad1->SetLeftMargin(0.15);
        pad1->SetBottomMargin(0.15);
        h1->SetTitleSize(0.06, "XY");
        h1->SetLabelSize(0.06, "XY");

        TGraphErrors *energy_graph = new TGraphErrors();
        for(int point = 0; point < true_energy_all.size(); point++){
            energy_graph->SetPoint(point, true_energy_all[point], ch_mean_all[point]);
            energy_graph->SetPointError(point, 0., ch_error_all[point]);
        }
        energy_graph->SetTitle("Calibrated energy");
        energy_graph->Draw();
        energy_graph->SetMarkerStyle(20);
        energy_graph->Draw("PSAME");

        energy_graph->Fit("pol1", "", "", 0, 10000);
        TF1 *func_calib = energy_graph->GetFunction("pol1");
        Double_t p0 = func_calib->GetParameter(0);
        Double_t p1 = func_calib->GetParameter(1);
        Double_t intercept = -p0 / p1;
        Double_t slope = 1.0 / p1;

        TVirtualPad *pad2 = c_calibfit->cd(2);
        // g_res->SetTitle(";True Energy (keV);(E_{obs} - E_{fit})/E_{fit}");
        pad2->SetLeftMargin(0.15);
        pad2->SetBottomMargin(0.15);

        // Calculate the residual -- the difference between the expected energy and the energy the calibration returns
        for (Int_t i = 0; i < true_energy_all.size(); i++) {
            Double_t ch_exp = p0 + p1 * true_energy_all[i];
            res_all.push_back((ch_mean_all[i] - ch_exp) / ch_exp);
            res_error_all.push_back(ch_error_all[i] / ch_exp);
        }

        TGraphErrors *g_res = new TGraphErrors();

        double max_res =-9999;
        double min_res = 9999;

        for(int point = 0; point < true_energy_all.size(); point++){
            g_res->SetPoint(point, true_energy_all[point], res_all[point]);
            g_res->SetPointError(point, 0., res_error_all[point]);
            if(res_all[point] + res_error_all[point] > max_res){
                max_res = res_all[point] + res_error_all[point];
            }
            if(res_all[point] - res_error_all[point] < min_res){
                min_res = res_all[point] - res_error_all[point];
            }
        }

        TH1F *h2 = pad2->DrawFrame(0.0, min_res + 0.1 * min_res, 10000.0, max_res + 0.1 * max_res);
        pad2->SetGrid();
        h2->SetTitleSize(0.06, "XY");
        h2->SetYTitle("(E_{obs} - E_{fit})/E_{fit}");
        h2->SetXTitle("True energy (keV)");
        h2->SetLabelSize(0.06, "XY");
        h2->SetTitleOffset(1.0, "Y");

        g_res->SetMarkerStyle(20);
        g_res->Draw("PSAME");

        // Write the calibration plot to the ROOT file
        c_calibfit->Write();

        // Clear all of the vectors so they are empty for the next file
        res_all.clear();
        res_error_all.clear();
        ch_mean_all.clear();
        ch_error_all.clear();
        true_energy_all.clear();

        //============================================================================================
        //============================================================================================
        //============================================================================================

        //--------------------------------------------------------------------------------------------

        //============================================================================================
        //============================================================================================
        //============================================================================================
        Double_t e_min = intercept;
        Double_t e_max = intercept + (Double_t)nbins * slope;
        TCanvas *c_calib = new TCanvas("ccalib", "ccalib", 600, 600);

        c_calib->SetTitle(ge_data_files[file_index].c_str());
        std::string hist_name = "hcalib_" + ge_data_files[file_index];

        TH1D* hists_calib = new TH1D(hist_name.c_str(), hist_name.c_str(), nbins, e_min, e_max);
        for (Int_t ibin = 0; ibin < nbins; ibin++) {
            hists_calib->SetBinContent(ibin + 1, hists[file_index]->GetBinContent(ibin + 1));
        }

        hists_calib->Write();

        std::cout << "==========================================" << std::endl;
        std::cout << "Slope = " << slope << std::endl;
        std::cout << "Intercept = " << intercept << std::endl;
        std::cout << "==========================================" << std::endl;
        //============================================================================================
        //============================================================================================
        //============================================================================================
    }



    file->Close();
    return;

}
