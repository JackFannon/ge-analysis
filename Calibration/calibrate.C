#include <forward_list>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
//#include "../all_headers.h"
#include "utilities.C"

//==================================================
//================ LOAD FILE NAMES==================
//==================================================

std::vector<std::string> input_data;

const std::string DATA_DIR = "/Users/jack/Software/GeAnalysis/Data/2023/raw/";

const int nhists = 23;

int data_type = 0;

bool Ni_sum_flag = false;

std::string source_name[3] = {
  "K40",
  "Co60",
  "NiCf"
};

int run_number[nhists] = {
    9999,
};

void calibrate(){

    std::vector<std::string> ge_data_files = load_data("data_list.txt", "/Users/jack/Software/GeAnalysis/Data/2023/");

    TFile* file = new TFile("histograms.root", "RECREATE");
    for(std::string file_name: ge_data_files){
        plot_channel_hist(file_name, DATA_DIR)->Write();
    }
    TH1D* hists[ge_data_files.size()];

    std::vector<std::string> isotope_search_names;
    isotope_search_names.push_back("K");
    isotope_search_names.push_back("Co");
    isotope_search_names.push_back("Cs");
    isotope_search_names.push_back("Ni");


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
        for(int source_label = 0; source_label < isotope_search_names.size(); source_label++){
            if(ge_data_files[file_index].find(isotope_search_names[source_label]) != std::string::npos){
                isotopes.push_back(source_label);
            }
        }
        isotopes_in_file.push_back(isotopes);
    }


    for(std::vector<int> j: isotopes_in_file){
        std::cout << "isotope type: ";
        for(int i: j){
            std::cout << i << ", ";
        }
        std::cout << std::endl;
    }

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

//==========================================================================================
//=========================== EXTRACT REGION OF INTEREST DATA ==============================
//==========================================================================================
    std::string roi_filename;
    roi_filename = "roi.txt";

    std::string buffer;
    std::ifstream roi_data;

    roi_data.open(roi_filename.c_str());

    if(!roi_data.is_open()){
        std::cout << "!!!!! Cannot find " << roi_filename.c_str() << "!!!!!" << std::endl;
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

    if(save_fit){
        TCanvas* c_dummy = new TCanvas();
        c_dummy->Print("fit_results.ps","Portrait");
    }


    double mean;
    double error;
    double K40_peak_ref = -1;
    std::string hist_names[ge_data_files.size()];
    std::string hist_titles[nhists];

    TH1D* Ni_sum = new TH1D("Co60+Ni_sum", "Co60+Ni_sum", nbins, 0, nbins);

    //plot_channel_hist(file_name[0]);
    //return;

    // Loop over the data files.
    for(int file_index = 0; file_index < ge_data_files.size(); file_index++){
        hist_names[file_index] = "h_" + ge_data_files[file_index];
        hist_titles[file_index] = ge_data_files[file_index] + ";Channel;Count";
        hists[file_index] = new TH1D(hist_names[file_index].c_str(), hist_titles[file_index].c_str(), nbins, 0, nbins);

        read_data_into_hist(DATA_DIR + ge_data_files[file_index], hists[file_index]);

        // Loop over isotope peaks that should be present in the data for this file
        for (int isotope_type: isotopes_in_file[file_index]){



            for(int isotope_peak = 0; isotope_peak < roi_high[isotope_type].size(); isotope_peak++){
                std::string canvas_name = std::to_string(isotope_peak) + isotope_search_names[isotope_type] + ge_data_files[file_index] + "fit";
                TCanvas* my_canvas = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 600, 600);

                fit_peak_ge(hists[file_index], roi_low[isotope_type][isotope_peak], roi_high[isotope_type][isotope_peak], &mean, &error);
                my_canvas->SetTitle(ge_data_files[file_index].c_str());
                hists[file_index]->Draw();
                true_energy_all.push_back(true_energy[isotope_type][isotope_peak]);
                ch_mean_all.push_back(mean);
                ch_error_all.push_back(error);

                if(save_fit){
                    hists[file_index]->SetAxisRange(roi_low[isotope_type][isotope_peak] - 50, roi_high[isotope_type][isotope_peak] + 50);
                    my_canvas->Write();
                }
            }
        }
        TGraphErrors *energy_graph = new TGraphErrors();

        for(int point = 0; point < true_energy_all.size(); point++){
            energy_graph->SetPoint(point, true_energy_all[point], ch_mean_all[point]);
            energy_graph->SetPointError(point, 0., ch_error_all[point]);
        }

        energy_graph->SetTitle("");

        energy_graph->Draw();

       std::string print_name = "calib" + ge_data_files[file_index];

        TCanvas *c_calibfit = new TCanvas(print_name.c_str() , print_name.c_str(), 600, 600);
        //  c_calibfit->Divide(1,2,0,0);
        c_calibfit->Divide(1, 2);

        TVirtualPad *pad1 = c_calibfit->cd(1);

        Double_t ch_max = 3000;
        if (data_type == 1){
            ch_max = 32000;
        }else if (data_type == 2){
            ch_max = 16000;
        }

        TH1F *h1 = pad1->DrawFrame(0, 0, 10000, ch_max, ";True energy (keV); MCA channel");
        pad1->SetLeftMargin(0.15);
        pad1->SetBottomMargin(0.15);
        h1->SetTitleSize(0.06, "XY");
        h1->SetLabelSize(0.06, "XY");
        energy_graph->SetMarkerStyle(20);
        energy_graph->Draw("PSAME");
        // energy_graph->Fit("pol1","","",0,2900); // set fit range in order to avoid Ni
        // points... energy_graph->Fit("pol1","","",3000,10000);
        energy_graph->Fit("pol1", "", "", 0, 10000); // set fit range in order to include only Ni points...
        // energy_graph->Fit("pol1");

        TF1 *func_calib = energy_graph->GetFunction("pol1");

        Double_t p0 = func_calib->GetParameter(0);
        Double_t p1 = func_calib->GetParameter(1);
        Double_t intercept = -p0 / p1;
        Double_t slope = 1.0 / p1;

        // calculate residual
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

        // g_res->SetTitle(";True Energy (keV);(E_{obs} - E_{fit})/E_{fit}");
        TVirtualPad *pad2 = c_calibfit->cd(2);
        pad2->SetLeftMargin(0.15);
        pad2->SetBottomMargin(0.15);
    // TH1F * h2 = pad2->DrawFrame(0,-0.03,10000,0.03,";True energy (keV);
    // (E_{obs} - E_{fit})/E_{fit}");
        TH1F *h2 = pad2->DrawFrame(0.0, min_res + 0.1 * min_res, 10000.0, max_res + 0.1 * max_res);
        pad2->SetGrid();
        h2->SetTitleSize(0.06, "XY");
        h2->SetYTitle("(E_{obs} - E_{fit})/E_{fit}");
        h2->SetXTitle("True energy (keV)");
        h2->SetLabelSize(0.06, "XY");
        h2->SetTitleOffset(1.0, "Y");

        g_res->SetMarkerStyle(20);
        g_res->Draw("PSAME");

        c_calibfit->Write();
        res_all.clear();
        res_error_all.clear();
        ch_mean_all.clear();
        ch_error_all.clear();
        true_energy_all.clear();


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
    }
    file->Close();
    return;













//=================================================================================================================

/*

    for (int i = 0; i < ge_data_files.size(); i++) {
        hist_names[i] = "h_"+source_name[i];
        hist_titles[i] = source_name[i]+";Channel;Count";
        hists[i] = new TH1D(hist_names[i].c_str(), hist_titles[i].c_str(), nbins, 0, nbins);

        read_data_into_hist(DATA_DIR + ge_data_files[i], hists[i]);

        if(Ni_sum_flag && source_type[i] == 3){
            Ni_sum->Add(hists[i], 1.0);
            continue;
        }

        my_canvas->SetTitle(source_name[i].c_str());

        hists[i]->Draw();

        if(save_fit){
            hists[i]->SetAxisRange(0, 30000);
            my_canvas->SetLogy();
            my_canvas->Print("fit_results.ps", "Portrait");
            my_canvas->SetLogy(0);
        }


        // Search for the K40, Tl208 and the two C60 peaks
       for (int source = 0; source <= 1; source++){
           for (int isotope = 0; isotope <= roi_high[i].size(); isotope++){
               fit_peak_ge(hists[i], roi_low[source][isotope], roi_high[source][isotope], &mean, &error);

               fit_results << run_number[i] << source << mean << " " << error << std::endl;

               if(save_fit){
                   my_canvas->Print("fit.ps", "Portrait");
               }
           }
       }


        if(save_fit){
            hists[i]->SetAxisRange(roi_low[0][0] - 50, roi_high[0][0] + 50);
            my_canvas->Print("fit_results.ps", "Portrait");
        }

        double K40_peak = mean;
        double K40_peak_error = error;

        fit_peak_ge(hists[i], roi_low[0][1], roi_high[0][1], &mean, &error);
        double Tl208_peak = mean;
        double Tl208_peak_error = error;

        fit_results << run_number[i] << " 0 " << K40_peak << " " << K40_peak_error << std::endl;
        fit_results << run_number[i] << " 1 " << Tl208_peak << " " << Tl208_peak_error << std::endl;

        if(K40_peak_ref < 0){
            K40_peak_ref = mean;
        }

        for (int ip = 0; ip < true_energy[source_type[i]].size(); ip++) {
            std::cout << "Fitting " << i << " " << roi_low[source_type[i]][ip] << " " << roi_high[source_type[i]][ip] << std::endl;
            fit_peak_ge(hists[i], roi_low[source_type[i]][ip],
                        roi_high[source_type[i]][ip], &mean, &error);

            if (save_fit) {
              hists[i]->SetAxisRange(roi_low[source_type[i]][ip] - 50,
                                 roi_high[source_type[i]][ip] + 50);

              my_canvas->Print("fit_results.ps", "Portrait");
            }

            hists[i]->SetAxisRange(0, nbins);

            true_energy_all.push_back(true_energy[source_type[i]][ip]);
            if (K40norm) {
              ch_mean_all.push_back(mean / K40_peak * K40_peak_ref);
              ch_error_all.push_back(error / K40_peak * K40_peak_ref);
            } else {
              ch_mean_all.push_back(mean);
              ch_error_all.push_back(error);
            }
            if (source_type[i] == 1) { // Co60
              fit_results << run_number[i] << " " << ip + 2 << " " << mean
                            << " " << error << std::endl;
            }
        }
    }

    fit_results.close();

    // graph with energy in X-axis (which I think is more natural)
    TGraphErrors *g =
        new TGraphErrors(true_energy_all.size(), &(true_energy_all[0]),
                         &(ch_mean_all[0]), 0, &(ch_error_all[0]));
    energy_graph->SetTitle("");

    // TGraphErrors * g = new TGraphErrors (true_energy_all.size(),
    //				       &(ch_mean_all[0]), &(true_energy_all[0]),
    //				       &(ch_error_all[0]),0);

    // energy_graph->SetTitle("");

    TCanvas *c_calibfit = new TCanvas("c_calibfit", "c_calibfit", 600, 600);
    //  c_calibfit->Divide(1,2,0,0);
    c_calibfit->Divide(1, 2);

    TVirtualPad *pad1 = c_calibfit->cd(1);

    Double_t ch_max = 3000;
    if (data_type == 1)
      ch_max = 32000;
    else if (data_type == 2)
      ch_max = 16000;

    TH1F *h1 =
        pad1->DrawFrame(0, 0, 10000, ch_max, ";True energy (keV); MCA channel");

    pad1->SetLeftMargin(0.15);
    pad1->SetBottomMargin(0.15);
    h1->SetTitleSize(0.06, "XY");
    h1->SetLabelSize(0.06, "XY");
    energy_graph->SetMarkerStyle(20);
    energy_graph->Draw("PSAME");
    //  energy_graph->Fit("pol1","","",0,2900); // set fit range in order to avoid Ni
    //  points... energy_graph->Fit("pol1","","",3000,10000);
    energy_graph->Fit("pol1", "", "", 0,
           10000); // set fit range in order to include only Ni points...
    // energy_graph->Fit("pol1");

    TF1 *func_calib = energy_graph->GetFunction("pol1");

    Double_t p0 = func_calib->GetParameter(0);
    Double_t p1 = func_calib->GetParameter(1);
    Double_t intercept = -p0 / p1;
    Double_t slope = 1.0 / p1;

    // calculate residual
    for (Int_t i = 0; i < true_energy_all.size(); i++) {
      Double_t ch_exp = p0 + p1 * true_energy_all[i];
      res_all.push_back((ch_mean_all[i] - ch_exp) / ch_exp);
      res_error_all.push_back(ch_error_all[i] / ch_exp);
    }

    TGraphErrors *g_res =
        new TGraphErrors(true_energy_all.size(), &(true_energy_all[0]),
                         &(res_all[0]), 0, &(res_error_all[0]));
    // g_res->SetTitle(";True Energy (keV);(E_{obs} - E_{fit})/E_{fit}");
    TVirtualPad *pad2 = c_calibfit->cd(2);
    pad2->SetLeftMargin(0.15);
    pad2->SetBottomMargin(0.15);
    // TH1F * h2 = pad2->DrawFrame(0,-0.03,10000,0.03,";True energy (keV);
    // (E_{obs} - E_{fit})/E_{fit}");
    TH1F *h2 = pad2->DrawFrame(0.0, -0.003, 10000.0, 0.003);
    pad2->SetGrid();
    h2->SetTitleSize(0.06, "XY");
    h2->SetYTitle("(E_{obs} - E_{fit})/E_{fit}");
    h2->SetXTitle("True energy (keV)");
    h2->SetLabelSize(0.06, "XY");
    h2->SetTitleOffset(1.0, "Y");

    g_res->SetMarkerStyle(20);
    g_res->Draw("PSAME");

    Double_t e_min = intercept;
    Double_t e_max = intercept + (Double_t)nbins * slope;

    c_calibfit->Print("linear_1.pdf");

    TCanvas *c_calib = new TCanvas("ccalib", "ccalib", 600, 600);

    for (Int_t i = 0; i < nhists; i++) {
      c_calib->SetTitle(source_name[i].c_str());

      hists_calib[i] = new TH1D(("hcalib_" + source_name[i]).c_str(),
                            source_name[i].c_str(), nbins, e_min, e_max);
      for (Int_t ibin = 0; ibin < nbins; ibin++) {
        hists_calib[i]->SetBinContent(ibin + 1, hists[i]->GetBinContent(ibin + 1));
      }

      hists_calib[i]->Draw();
    }

    if(save_fit){
      my_canvas->Print("fit_results.ps]", "Portrait");
    }

    std::cout << "==========================================" << std::endl;
    std::cout << "Slope = " << slope << std::endl;
    std::cout << "Intercept = " << intercept << std::endl;
    std::cout << "==========================================" << std::endl;

    */
}
