#include "../Calibration/utilities.C"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom.h"
#include "TText.h"
#include "TVirtualPad.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// OPTIONS FOR CROSS CALIBRATION CHECKING
// const double calib_const[2] = {-8.084, 0.1915};

//=================================================================================================================
//======================================== CONFIG OPTIONS =========================================================
//=================================================================================================================
const int nbins = 4096;
const int data_type = 0;

const bool smear_flag = false;

// New detector main data
// const double calib_const[2] = {-0.166597, 0.148259};

const double calib_const[2] = {-2.7842, 0.148073};

const double intercept = -calib_const[0] / calib_const[1];
const double slope = 1.0 / calib_const[1];
const double e_min = 0.001 * intercept;
const double e_max = 0.001 * (intercept + (double)nbins * slope);

const std::vector<float> source_e_true = {1.4608, 2.6145, 1.1732, 1.3325};

float chi2_min_buffer = 99999.;

//=================================================================================================================
//=============================== CALCULATE CHI2 BETWEEN MC AND DATA ==============================================
//=================================================================================================================
float calc_chi2(TH1F *h_data, TH1F *h_mc, double e_min, double e_max, std::string output_filename,
                bool plot_flag = false) {
    // Find bin range where the LINAC peak should be
    int bin_min = h_data->FindBin(e_min);
    int bin_max = h_data->FindBin(e_max);

    // Print out bin and energy values
    std::cout << "================ " << e_min << " " << e_max << " " << bin_min << " " << bin_max << std::endl;

    // Get the integral where the peak should be to normalise
    double norm_data = h_data->Integral(bin_min, bin_max);
    double norm_mc = h_mc->Integral(bin_min, bin_max);

    // Check that axes are consistent
    double xmax_data = h_data->GetXaxis()->GetXmax();
    double xmax_mc = h_mc->GetXaxis()->GetXmax();
    std::cout << "consistency check: " << xmax_data << " " << xmax_mc << std::endl;
    std::cout << "consistency check: " << h_data->GetXaxis()->GetXmin() << " " << h_mc->GetXaxis()->GetXmin()
              << std::endl;

    double hist_x_min = e_min - 0.1;
    double hist_x_max = e_max + 0.05;

    float chi2 = 0;
    double n_data = 0;
    double n_mc = 0;
    double err_data = 0;
    double err_mc = 0;

    int ndf = 0;

    TCanvas *c1 = new TCanvas("", "");

    // Loop over the bins and calculate the chi2 between the data and MC
    for (int i = bin_min; i <= bin_max; i++) {
        n_data = h_data->GetBinContent(i);
        if (n_data == 0) {
            continue;
        } // skip if number of events in data is 0
        err_data = sqrt(n_data);
        n_mc = h_mc->GetBinContent(i) * norm_data / norm_mc;
        err_mc = sqrt(h_mc->GetBinContent(i)) * norm_data / norm_mc;
        chi2 += pow(n_data - n_mc, 2) / (err_data * err_data + err_mc * err_mc);
        ndf++;
    }

    // Plotting
    std::string dummy;
    if (plot_flag) {
        int nbins = h_mc->GetNbinsX();
        for (int i = 0; i < nbins; i++) {
            n_mc = h_mc->GetBinContent(i + 1) * norm_data / norm_mc;
            err_mc = sqrt(h_mc->GetBinContent(i + 1)) * norm_data / norm_mc;
            h_mc->SetBinContent(i + 1, n_mc);
            h_mc->SetBinError(i + 1, err_mc);
        }
        int max_bin_mc = h_mc->GetBinContent(h_mc->GetMaximumBin());
        c1->cd();
        h_data->GetXaxis()->SetRangeUser(hist_x_min, hist_x_max);
        h_data->GetXaxis()->SetTitle("Energy [MeV]");
        h_data->GetYaxis()->SetTitle("Counts");
        h_data->SetMinimum(0);
        // h_mc->GetYaxis()->SetRangeUser(0,1.1 * h_data->GetMaximum());
        h_mc->SetLineColor(2);
        h_mc->SetLineWidth(2);
        h_data->SetLineWidth(2);
        h_data->SetLineColor(1);
        h_data->GetYaxis()->SetRangeUser(0, 1.2 * max_bin_mc);
        h_data->Draw("e");
        h_mc->Draw("hist same");
        h_data->Draw("e same");
        TLine *l1 = new TLine(e_min, 0, e_min, h_mc->GetMaximum());
        TLine *l2 = new TLine(e_max, 0, e_max, h_mc->GetMaximum());
        l1->SetLineColor(4);
        l2->SetLineColor(4);
        l1->SetLineWidth(2);
        l2->SetLineWidth(2);
        l1->Draw();
        l2->Draw();
        c1->Update();
        std::cout << "chi2/NDF = " << chi2 << " / " << ndf - 1 << std::endl;
    }

    std::cout << "MIN =========================== " << chi2_min_buffer << " chi2 ===================== " << chi2
              << std::endl;
    if (chi2_min_buffer == 0 || chi2_min_buffer > chi2) {
        chi2_min_buffer = chi2;
        c1->Print((output_filename + ".pdf").c_str());
    }
    std::cout << std::endl;

    return chi2;
}

//=================================================================================================================
//==================================== APPLY Ni CALIBRATION TO OTHER DATA FILES ===================================
//=================================================================================================================
TH1F *apply_calibration(std::string raw_data_filename) {
    // Histogram for raw data
    TH1F *h_data_raw = new TH1F("h_data", "h_data", nbins, 0, nbins);

    // Histogram for the data once Ni calibration has been applied
    TH1F *h_data_calib = new TH1F("h_calib", "h_calib", nbins, e_min, e_max);

    // Fill a histogram with the channel counts
    read_data_into_hist(raw_data_filename, h_data_raw);

    // Set the bin content of calibrated histogram
    for (int bin = 0; bin < nbins; bin++) {
        h_data_calib->SetBinContent(bin + 1, h_data_raw->GetBinContent(bin + 1));
        h_data_calib->SetBinError(bin + 1, sqrt(h_data_raw->GetBinContent(bin + 1)));
    }

    h_data_calib->GetYaxis()->SetTitle("Counts");
    h_data_calib->GetXaxis()->SetTitle("Energy [MeV]");

    h_data_calib->Draw();

    delete h_data_raw;

    return h_data_calib;
}

//=================================================================================================================
//============================================== SETUP LEGEND =====================================================
//=================================================================================================================
void setup_legend(TLegend *leg, TH1F *input_hist, std::string legend_label, double energy, double momentum,
                  float chi2) {
    leg->SetFillColor(0);
    leg->AddEntry(input_hist, legend_label.c_str(), "1");
    TLatex *energy_label = new TLatex(0.15, 0.8, Form("E_{tot} = %4.3f MeV", energy));
    TText *momentum_label = new TText(0.15, 0.7, Form("P = %4.3f MeV", momentum));
    TLatex *chi2_label = new TLatex(0.15, 0.6, Form("#chi^{2} = %2.1f", chi2));
    energy_label->SetNDC();
    momentum_label->SetNDC();
    chi2_label->SetNDC();
    energy_label->Draw();
    momentum_label->Draw();
    chi2_label->Draw();
}

//=================================================================================================================
//================================================ CALCULATE SMEARING =============================================
//=================================================================================================================
TH1F *smear_mc(TH1F *h_mc, int energy, double source_e_error[source_e_true.size()], double fit_e_min,
               double fit_e_max) {
    // Intialise histogram to store smeared MC
    TH1F *h_mc_smeared = new TH1F("smeared_mc", "smeared_mc", nbins, e_min, e_max);

    // Find the maximum and minimum bin numbers where there will be data
    int bin_min = h_mc->FindBin(fit_e_min);
    int bin_max = h_mc->FindBin(fit_e_max);

    // Create integer to hold MC bin content before smearing
    int mc_entry_bin = 0;

    // Loop over bins and smear with a gaussian
    for (int i = bin_min; i < bin_max; i++) {
        mc_entry_bin = h_mc->GetBinContent(i + 1);
        for (int j = 0; j < 200; j++) {
            double random = gRandom->Gaus(h_mc->GetXaxis()->GetBinCenter(i + 1), 3 * source_e_error[2]);
            double smearing_factor = double(mc_entry_bin) / 200.0;
            h_mc_smeared->Fill(random, smearing_factor);
        }
    }

    return h_mc_smeared;
}

void fit_linac(std::string data_filename, int data_type, int x_min, int x_max, double fit_e_min, double fit_e_max,
               std::string output_filename, int beam_energy, std::string xpos, std::string zpos) {

    std::vector<float> chi2;

    TCanvas *my_canvas = new TCanvas("my_canvas", "my_canvas", 600, 600);

    // Set drawing options using utilities function
    set_style(132);

    // Initalise arrays for the mean and error values of the fits for each radiation source
    double source_e_mean[source_e_true.size()];
    double source_e_error[source_e_true.size()];

    TH1F *h_data_calib = apply_calibration(data_filename);

    // Fit the Co60, K40 and Tl208 peaks from each data file.
    for (int i = 0; i < source_e_true.size(); i++) {
        std::cout << source_e_true[i] << std::endl;
        fit_peak_ge(h_data_calib, 0.991 * source_e_true[i], 1.009 * source_e_true[i], &source_e_mean[i],
                    &source_e_error[i]);
    }

    std::vector<double> chi2_vec;
    std::vector<double> smeared_chi2_vec;
    std::vector<double> x_vec;
    std::vector<double> p_vec;

    std::vector<std::string> filenames;

    double chi2_min[2] = {1e10, 1e10};
    int x_best[2] = {-1, -1};
    double p_best[2] = {0, 0};

    chi2_min_buffer = 99999.;
    // Loop over the energy range between x_min and x_max
    for (int x = x_min; x < x_max; x++) {
        std::string mc_filename;
        // Open MC file
        if (x >= 10000) {
            mc_filename = "../MC/CrossCalibration/OldDetector/" + std::to_string(x) + ".root";
        } else {
            mc_filename = "../MC/CrossCalibration/OldDetector/0" + std::to_string(x) + ".root";
        }

        // TFile for the MC root file
        TFile *mc_file = new TFile(mc_filename.c_str(), "READ");

        // Get the histogram of total energy deposition in the Ge detector from the MC file
        TH1F *h_mc = (TH1F *)mc_file->Get("h21");

        h_mc->Draw();

        // Smear the MC histogram
        TH1F *h_mc_smeared = smear_mc(h_mc, x, source_e_error, fit_e_min, fit_e_max);

        // Energy in MeV
        double etot = 0.001 * (double)x;

        // Momentum conversion
        double p = sqrt(pow(etot, 2) - pow(0.511, 2));

        // Calculate how well the MC compares to the data
        chi2.push_back(calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, output_filename, true));
        if (smear_flag) {
            chi2.push_back(calc_chi2(h_data_calib, h_mc_smeared, fit_e_min, fit_e_max, output_filename, true));
        }

        // Draw the legend and save the canvas
        TLegend *leg = new TLegend(0.15, 0.4, 0.35, 0.55);

        setup_legend(leg, h_mc, "Default MC", etot, p, chi2[0]);
        if (smear_flag) {
            setup_legend(leg, h_mc_smeared, "Smeared MC", etot, p, chi2[1]);
        }

        // Push back the energy and momenta into a vector (?)
        //   not sure why this is done
        x_vec.push_back((double)x);
        p_vec.push_back((double)p);

        // Push back the chi2 value for each case, smeared and not smeared into the
        //   corresponding vectors
        chi2_vec.push_back(chi2[0]);
        smeared_chi2_vec.push_back(chi2[1]);
        // Use the smallest chi2 value to define the best energy/momenta match
        if (chi2[0] < chi2_min[0]) {
            chi2_min[0] = chi2[0];
            x_best[0] = x;
            p_best[0] = p;
        }
        // Do the same for the smeared case
        if (chi2[1] < chi2_min[1]) {
            chi2_min[1] = chi2[1];
            x_best[1] = x;
            p_best[1] = p;
        }

        // Remove all entries from the chi2 vector
        chi2.clear();

        // Delete histograms
        delete h_mc;
        delete h_mc_smeared;
        delete leg;

        // Close the file
        mc_file->Close();
    }

    TCanvas *my_third_canvas = new TCanvas("c3", "c3", 600, 600);
    TGraph *my_graph_p;
    if (smear_flag) {
        my_graph_p = new TGraph(smeared_chi2_vec.size(), &p_vec[0], &smeared_chi2_vec[0]);
    } else {
        my_graph_p = new TGraph(chi2_vec.size(), &p_vec[0], &chi2_vec[0]);
    }

    // Draw graph of momentum against chi^2
    my_graph_p->SetTitle(";Momentum (MeV);#chi^{2}");
    my_graph_p->Draw("ALP");
    my_third_canvas->SaveAs((output_filename + ".root").c_str());

    std::string mc_filename;

    // TFile for the MC root file
    if (x_best[0] >= 10000) {
        mc_filename = "../MC/CrossCalibration/OldDetector/" + std::to_string(x_best[0]) + ".root";
    } else {
        mc_filename = "../MC/CrossCalibration/OldDetector/0" + std::to_string(x_best[0]) + ".root";
    }

    TFile *file_mc = new TFile(mc_filename.c_str(), "READ");

    // Get the histogram of total energy deposition in the Ge detector from the MC file
    TH1F *h_mc = (TH1F *)file_mc->Get("h21");

    TH1F *h_smc = new TH1F("smearing", "smearing", nbins, e_min, e_max);
    if (smear_flag) {
        smear_mc(h_mc, x_best[0], source_e_error, fit_e_min, fit_e_max);
    }

    chi2.push_back(calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, output_filename, true));
    if (smear_flag) {
        chi2.push_back(calc_chi2(h_data_calib, h_smc, fit_e_min, fit_e_max, output_filename, true));
    }

    std::cout << "Best fit total energy = " << x_best[0] << " (keV)" << std::endl;
    std::cout << "Best fit momentum = " << p_best[0] << " (MeV)" << std::endl;
    if (smear_flag) {
        std::cout << "Best fit total energy = " << x_best[1] << " (keV)" << std::endl;
        std::cout << "Best fit momentum = " << p_best[1] << " (MeV)" << std::endl;
    }

    std::ofstream ofs;
    ofs.open("bestfit_momentum_tmp.txt", std::ios::app);

    std::ofstream diff_out;
    diff_out.open("diff_default_smeared_w_new_nicalib.txt", std::ios::app);

    diff_out << xpos << " \t" << zpos << "\t" << beam_energy << "\t" << 0.001 * x_best[0] << "\t"
             << 0.001 * x_best[1] << "\t" << chi2[0] << std::endl;

    ofs << p_best[0];
    if (smear_flag) {
        ofs << "    " << p_best[1];
    }
    ofs << std::endl;

    ofs.close();
    diff_out.close();

    // Delete h_data_calib as it is not used past here
    delete h_data_calib;
    delete h_smc;
    delete h_mc;
    delete my_canvas;
    delete my_third_canvas;

    file_mc->Close();
}

//=================================================================================================================
//=================================== MACRO ENTRANCE -- WRAPPER FOR FIT_LINAC =====================================
//=================================================================================================================
void fitter(std::string data_info, std::string output_filename) {
    // Information required by the fitting function above is as follows:
    // std::string data_filename   - Filename of Ge detector data
    //         int data_type       - Type of data (refers to the type of detector)
    //         int x_min           - Minimum energy value for MC
    //         int x_max           - Maximum energy value for MC
    //      double fit_e_min       - Lower bound for data peak
    //      double fit_e_max       - Upper bound for data peak
    // std::string output_filename - Output name for plots
    //         int beam_energy     - Approximate linac beam energy
    // std::string xpos            - Approximate x position of beam
    // std::string zpos            - Approximate z position of beam

    // Intialise variables to read out the runlist file:
    int run_number;
    int approx_energy;
    int data_type;
    std::string approx_x;
    std::string approx_z;
    std::string filename;
    std::string output_name;
    int min_mc_energy;
    int max_mc_energy;
    float min_data_energy;
    float max_data_energy;

    // Load in the input file list
    std::ifstream input_file(data_info);

    // Check that input file list exists and has opened
    if (!input_file) {
        std::cerr << "Error opening file " << data_info << std::endl;
        return;
    }

    // String to load the file contents into
    std::string line;

    // Loop over the file contents
    while (std::getline(input_file, line)) {
        // If the line is empty of starts with a # then skip
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Buffer to convert between string and the variables
        std::istringstream buffer(line);

        // Load information into local variables
        buffer >> run_number >> approx_energy >> data_type >> approx_x >> approx_z >> filename >> output_name >>
            min_mc_energy >> max_mc_energy >> min_data_energy >> max_data_energy >> calib_p0 >> calib_p1;

        // Store calibration constants in an array
        float calib_const[2] = {calib_p0, calib_p1};

        // Create a string that matches the input data filename, but without the ".csv" on the end
        std::string filename_wo_csv;
        if (filename.substr(filename.length() - 4) == ".csv") {
            filename_wo_csv = filename.substr(0, filename.length() - 4);
        }

        // Call fit_linac to find the best match between MC and data
        fit_linac(filename, data_type, min_mc_energy, max_mc_energy, min_data_energy, max_data_energy, output_name,
                  approx_energy, approx_x, approx_z, calib_const);
    }

    return;
}
