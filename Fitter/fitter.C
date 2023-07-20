#include "../Calibration/utilities.C"
#include <fstream>
#include <string>
#include <vector>

#include "headers.h"

const int nbins = 4096;
const int data_type = 0;

const bool smear_flag = false;

const double calib_const[2] = {-0.166597, 0.148259};
const double intercept = - calib_const[0]/calib_const[1];
const double slope = 1.0/calib_const[1];
const double e_min = 0.001 * intercept;
const double e_max = 0.001 * (intercept + (double)nbins * slope);

const std::vector<float> source_e_true = { 1.4608, 2.6145, 1.1732, 1.3325 };


const std::string mc_filename;

float calc_chi2(TH1D * h_data, TH1F * h_mc, double e_min, double e_max, bool plot_flag = false, TH1D * h_smc = new TH1D("", "", 100, 0, 100), bool flag = false){
    int bin_min = h_data->FindBin(e_min);
    int bin_max = h_data->FindBin(e_max);
    std::cout << "================ " << e_min << " " << e_max << " " << bin_min << " " << bin_max << std::endl;
    double norm_data = h_data->Integral(bin_min,bin_max);
    double norm_mc = h_mc->Integral(bin_min,bin_max);
    double norm_smc = h_smc->Integral(bin_min,bin_max);
    // consistency check:
    double xmax_data = h_data->GetXaxis()->GetXmax();
    double xmax_mc = h_mc->GetXaxis()->GetXmax();
    std::cout << "consistency check: " << xmax_data << " " << xmax_mc << std::endl;
    std::cout << "consistency check: " << h_data->GetXaxis()->GetXmin() << " "
         << h_mc->GetXaxis()->GetXmin() << std::endl;

    double hist_x_min = e_min - 0.1;
    double hist_x_max = e_max + 0.05;


    double chi2 = 0;
    double n_data = 0;
    double n_mc = 0;
    double n_smc = 0;
    double err_data = 0;
    double err_mc = 0;
    double err_smc = 0;

    int ndf = 0;

    TCanvas* c1 = new TCanvas("", "");
    for (int i = bin_min; i <= bin_max; i++){
        n_data = h_data->GetBinContent(i);
        if (n_data == 0){
            continue;
        }// skip if number of events in data is 0
        err_data = sqrt(n_data);
        n_mc = h_mc->GetBinContent(i) * norm_data / norm_mc;
        err_mc =  sqrt(h_mc->GetBinContent(i)) * norm_data / norm_mc;
        n_smc = h_smc->GetBinContent(i) * norm_data / norm_smc;
        err_smc =  sqrt(h_smc->GetBinContent(i)) * norm_data / norm_smc;
        if(flag){
            chi2 += pow(n_data - n_smc,2) / (err_data*err_data + err_smc*err_smc);
        }else{
            chi2 += pow(n_data - n_mc,2) / (err_data*err_data + err_mc*err_mc);
        }
        ndf++;
        //    std::cout << i << " " << n_data << " " << n_mc << std::endl;
    }
    int max_bin_mc = h_mc -> GetMaximumBin();
    double scalefactor;
    //h_smc -> Scale(scalefactor);

    std::string dummy;
    if (plot_flag){
        int nbins = h_mc->GetNbinsX();
        for (int i = 0; i < nbins; i++){
            n_mc = h_mc->GetBinContent(i+1) * norm_data / norm_mc;
            err_mc =  sqrt(h_mc->GetBinContent(i+1)) * norm_data / norm_mc;
            h_mc->SetBinContent(i+1, n_mc);
            h_mc->SetBinError(i+1, err_mc);
            n_smc = h_smc -> GetBinContent(i+1) * norm_data / norm_smc;
            err_smc = sqrt(h_smc->GetBinContent(i+1)) * norm_data / norm_smc;
            h_smc -> SetBinContent(i+1, n_smc);
        }
        c1->cd();
        h_data->GetXaxis()->SetRangeUser(hist_x_min,hist_x_max);
        h_data->GetXaxis()->SetTitle("Energy [MeV]");
        h_data->GetYaxis()->SetTitle("Counts");
        h_data->SetMinimum(0);
        //h_mc->GetYaxis()->SetRangeUser(0,1.1 * h_data->GetMaximum());
        h_mc->SetLineColor(2);
        h_mc->SetLineWidth(2);
        h_smc->SetLineColor(3);
        h_smc->SetLineWidth(2);
        h_data->SetLineWidth(2);
        h_data->SetLineColor(1);
        h_data->Draw("e");
        h_mc->Draw("hist same");
        if(flag){
            h_smc->Draw("hist same");
        }
        h_data->Draw("e same");

        TLine * l1 = new TLine(e_min,0,e_min,h_mc->GetMaximum());
        TLine * l2 = new TLine(e_max,0,e_max,h_mc->GetMaximum());
        l1->SetLineColor(4);
        l2->SetLineColor(4);
        l1->SetLineWidth(2);
        l2->SetLineWidth(2);
        l1->Draw();
        l2->Draw();

        c1->Update();
        std::cout << "chi2/NDF = " << chi2 << " / " << ndf-1 << std::endl;
    }

    return chi2;
}

void fitter(std::string data_filename,
            int data_type,
            int x_min,
            int x_max,
            double fit_e_min,
            double fit_e_max,
            std::string output_filename,
            int beam_energy,
            std::string xpos,
            std::string zpos){

    // Set drawing options using utilities function
    set_style(132);

    // Initiate the histogram for channel counts
    TH1D* h_data_raw = new TH1D("h_data", "h_data", nbins, 0, nbins);

    read_data_into_hist(data_filename, h_data_raw);

    TH1D* h_data_calib = new TH1D("", "", nbins, e_min, e_max);

    for (int bin = 0; bin < nbins; bin++) {
        h_data_calib->SetBinContent(bin + 1, h_data_raw->GetBinContent(bin+1));
        h_data_calib->SetBinError(bin + 1, sqrt(h_data_raw->GetBinContent(bin + 1)));
    }

    TCanvas* my_canvas = new TCanvas("c1", "c1");

    double source_e_mean[source_e_true.size()];
    double source_e_error[source_e_true.size()];

    for (int i = 0 ; i < source_e_true.size(); i++) {
        fit_peak_ge(h_data_calib, 0.975 * source_e_true[i], 1.025 * source_e_true[i], &source_e_mean[i], &source_e_error[i]);
        h_data_calib->GetYaxis()->SetTitle("Counts");
        h_data_calib->GetXaxis()->SetTitle("Energy [MeV]");
        h_data_calib->Draw();
        my_canvas->SaveAs(output_filename.c_str());
    }

    std::vector<double> chi2_vec;
    std::vector<double> smeared_chi2_vec;
    std::vector<double> x_vec;
    std::vector<double> p_vec;

    std::vector <std::string> filenames;

    double chi2_min[2] = {1e10, 1e10};
    int x_best[2] = {-1, -1};
    double p_best[2] = {0, 0};

    for (int x = x_min; x < x_max; x++) {

        TFile *file_mc = new TFile();
        file_mc->Open(mc_filename.c_str());

        TH1F *h_mc = (TH1F*)file_mc->Get("h21");

        int mc_entry_bin = 0;
        TH1D *h_smc = new TH1D("smeared", "smeared", nbins, e_min, e_max);

        int bin_min = h_mc->FindBin(fit_e_min) - 20;
        int bin_max = h_mc->FindBin(fit_e_max) + 10;

        for (int i = bin_min; i < bin_max; i++) {
            mc_entry_bin = h_mc->GetBinContent(i + 1);
            for (int j = 0; j < 200; j++) {
                double random = gRandom->Gaus(h_mc->GetXaxis()->GetBinCenter(i + 1),
                                        source_e_error[2]);
                double sf = double(mc_entry_bin) / 200.0;
                h_smc->Fill(random, sf);
            }
        }

        double etot = 0.001 * (double)x;
        double p = sqrt(pow(etot, 2) - pow(0.511, 2));

        std::vector<float> chi2;
        chi2.push_back(calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, true, h_smc, false));
        chi2.push_back(calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, true, h_smc, true));

        TLegend* leg = new TLegend(0.15, 0.4, 0.35, 0.55);
        leg->SetFillColor(0);
        leg->AddEntry(h_mc, "Default MC" , "1");
        if(smear_flag){
            leg->AddEntry(h_smc, "Smeared MC" , "1");
        }

        if(smear_flag) {
            TLatex* t_etot = new TLatex(0.15, 0.8, Form("E_{tot} = %4.3f MeV", etot));
            TText* t_p = new TText (0.15, 0.7, Form("P = %4.3f MeV", p));
            TLatex* t_chi2 = new TLatex (0.15, 0.6, Form("#chi^{2} = %2.1f", chi2[1]));
            t_etot->SetNDC();
            t_p->SetNDC();
            t_chi2->SetNDC();
            t_etot->Draw();
            t_p->Draw();
            t_chi2->Draw();
        }else{
            TLatex* t_etot = new TLatex (0.15, 0.8, Form("E_{tot} = %4.3f MeV",etot));
            TText* t_p = new TText (0.15, 0.7, Form("P = %4.3f MeV",p));
            TLatex* t_chi2 = new TLatex (0.15, 0.6, Form("#chi^{2} = %2.1f",chi2[0]));
            t_etot->SetNDC();
            t_p->SetNDC();
            t_chi2->SetNDC();
            t_etot->Draw();
            t_p->Draw();
            t_chi2->Draw();
        }
        leg->Draw();

        my_canvas->SaveAs(output_filename.c_str());

        std::cout << x << " " << chi2[0] << " " << chi2[1] << std::endl;

        x_vec.push_back((double)x);
        p_vec.push_back((Double_t)p);


        chi2_vec.push_back(chi2[0]);
        smeared_chi2_vec.push_back(chi2[1]);
        if (chi2[0] < chi2_min[0]){
            chi2_min[0] = chi2[0];
            x_best[0] = x;
            p_best[0] = p;
        }
        if (chi2[1] < chi2_min[1]){
            chi2_min[1] = chi2[1];
            x_best[1] = x;
            p_best[1] = p;
        }
        file_mc->Close();
    }

    TCanvas* my_other_canvas = new TCanvas("c2", "c2", 600, 600);
    TGraph* my_graph = new TGraph(chi2_vec.size(), &x_vec[0], &smeared_chi2_vec[0]);

    if(smear_flag){
        my_graph->GetXaxis()->SetRangeUser(x_best[1] - 15, x_best[1] + 15);
        my_graph->GetYaxis()->SetRangeUser(0, chi2_min[1] + 200);
    }else{
        my_graph->GetXaxis()->SetRangeUser(x_best[0] - 15, x_best[0] + 15);
        my_graph->GetYaxis()->SetRangeUser(0, chi2_min[0] + 200);
    }

    my_graph->SetTitle(";Total energy (keV);#chi^{2}");
    my_graph->Draw("AL");
    my_other_canvas->SaveAs(output_filename.c_str());

    TCanvas* my_third_canvas = new TCanvas("c3", "c3", 600, 600);
    TGraph* my_graph_p;
    if(smear_flag){
        my_graph_p = new TGraph(smeared_chi2_vec.size(), &p_vec[0], &smeared_chi2_vec[0]);
    }else{
        my_graph_p = new TGraph(chi2_vec.size(), &p_vec[0], &chi2_vec[0]);
    }

    my_graph_p->SetTitle(";Momentum (MeV);#chi^{2}");
    my_graph_p->Draw("AL");
    my_third_canvas->SaveAs(output_filename.c_str());

    if(smear_flag){
        my_graph_p->GetXaxis()->SetRangeUser(p_best[1]-0.015,p_best[1]+0.015);
        my_graph_p->GetYaxis()->SetRangeUser(0,chi2_min[1]+200);
    }else{
        my_graph_p->GetXaxis()->SetRangeUser(p_best[0]-0.015,p_best[0]+0.015);
        my_graph_p->GetYaxis()->SetRangeUser(0,chi2_min[0]+200);
    }

    my_graph_p->Draw("AL");
    my_third_canvas->Print(output_filename.c_str());

    TFile* file_mc = new TFile(mc_filename.c_str());
    TH1F* h_mc = (TH1F*)file_mc->Get("h21");
    TH1D* h_smc = new TH1D("smearing", "smearing", nbins, e_min, e_max);
    int bin_min = h_mc->FindBin(fit_e_min) - 20;
    int bin_max = h_mc->FindBin(fit_e_max) + 10;

    for(int i = bin_min; i < bin_max; i++){
        int mc_entry_bin = h_mc->GetBinContent(i+1);
        for(int j = 0; j < 2000; j++){
            double random = gRandom->Gaus(h_mc->GetXaxis()->GetBinCenter(i+1), source_e_error[2]);
            double sf = double(mc_entry_bin)/2000.0;
            h_smc->Fill(random, sf);
        }
    }

    std::vector<float> chi2;

    chi2.push_back(calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, true, h_smc, false));
    chi2.push_back(calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, true, h_smc, true));
    std::cout << "Best fit total energy = " << x_best[0] << " (keV)" <<  std::endl;
    std::cout << "Best fit momentum = " << p_best[0] << " (MeV)" <<  std::endl;
    std::cout << "Best fit total energy = " << x_best[1] << " (keV)" <<  std::endl;
    std::cout << "Best fit momentum = " << p_best[1] << " (MeV)" <<  std::endl;

    if(smear_flag){
        TLatex* t_etot = new TLatex (0.15, 0.8, Form("E_{tot} = %4.3f MeV",0.001*x_best[1]));
        TText* t_p = new TText (0.15, 0.7, Form("P = %4.3f MeV",p_best[1]));
        TLatex* t_chi2 = new TLatex (0.15, 0.6, Form("#chi^{2} = %2.1f",chi2_min[1]));
        t_etot->SetNDC();
        t_printf->SetNDC();
        t_chi2->SetNDC();
        t_etot->Draw();
        t_p->Draw();
        t_chi2->Draw();
    } else{
        TLatex* t_etot = new TLatex (0.15, 0.8, Form("E_{tot} = %4.3f MeV",0.001*x_best[0]));
        TText* t_p = new TText (0.15, 0.7, Form("P = %4.3f MeV",p_best[0]));
        TLatex* t_chi2 = new TLatex (0.15, 0.6, Form("#chi^{2} = %2.1f",chi2_min[0]));
        t_etot->SetNDC();
        t_printf->SetNDC();
        t_chi2->SetNDC();
        t_etot->Draw();
        t_p->Draw();
        t_chi2->Draw();
    }
    TLegend* leg = new TLegend(0.15,0.4,0.35,0.55);
    leg->SetFillColor(0);
    leg->AddEntry(h_mc, "Default MC" , "l");
    if(smear_flag){
        leg->AddEntry(h_smc, "Smeared MC" , "l");
    }
    leg->Draw();

    my_canvas->SaveAs(output_filename.c_str());

    std::ofstream ofs("bestfit_momentum_tmp.txt");
    std::ofstream diff_out("diff_default_smeared_w_new_nicalib.txt", ios::app);
    if(smear_flag){
        ofs << p_best[1];
    }else{
        ofs << p_best[0];
    }
    diff_out << xpos << " \t" << zpos << "\t" << beam_energy << "\t";
    diff_out << 0.001 * x_best[0] << "\t" << 0.001 * x_bestp[1] << "\t" << chi2[0] << std::endl;

    for(int i = 0; i < source_e_true.size(); i++){
        ofs << "\t" << source_e_mean[i] << "\t" << source_e_error[i];
    }
    ofs << endl;
    ofs.close();
    diff_out.close();
}
