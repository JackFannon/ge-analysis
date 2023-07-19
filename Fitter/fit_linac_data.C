#include "../Calibration/utilities.C"
TCanvas * c1;
Int_t stat_smc = 2000;

TRandom3 * gRandom = new TRandom3(0);
gRandom -> SetSeed(0);

Bool_t smear_flag = false;//if smearing is not needed, change to false

void fit_linac_data(string data_filename,
                    Int_t data_type,
                    Int_t x_min, Int_t x_max,
                    Double_t fit_e_min, Double_t fit_e_max,
                    string output_filename,
                    Int_t beam_energy,
                    string xpos,
                    string zpos){

  Int_t fontid = 132;
  gStyle -> SetOptStat(0);
  gStyle -> SetPadBorderSize(0);
  gStyle -> SetLegendBorderSize(0);
  gStyle -> SetFrameFillColor(0);
  gStyle ->  SetCanvasColor(0);
  gStyle->SetStatFont(fontid);
  gStyle->SetLabelFont(fontid, "XYZ");
  gStyle->SetLabelFont(fontid, "");
  gStyle->SetTitleFont(fontid, "XYZ");
  gStyle->SetTitleFont(fontid, "");
  gStyle->SetTitleOffset(1.2, "XYZ");
  gStyle->SetTextFont(fontid);

  //const Int_t nbins = 32768;
  Int_t nbins = 4096;
  if (data_type  == 1 || data_type == 2)
    nbins = 32768;

  /////////////////////////////////////////////////


  Int_t rebin_fac = 8; // rebin factor for new Ge data

  TH1D * h_data_raw =  new TH1D("h_data","h_data",nbins,0,nbins);

  // Parameters for 2019 Ge
  /*
  Double_t calib_const[3][2] = {
    {-6.58976, 0.191351},
      {5.86207, 3.09338},
      {7.39652, 1.58129}
  };
  */
  //parameters for result with Ni calibration
  Double_t calib_const[3][2] = {
    {-0.166597, 0.148259},
      {5.86207, 3.09338},
      {7.39652, 1.58129}
  };

  read_data_into_hist(data_filename, h_data_raw);

  Double_t intercept = - calib_const[data_type][0]/calib_const[data_type][1];
  Double_t slope = 1.0/calib_const[data_type][1];

  Double_t e_min = 0.001* intercept;
  Double_t e_max = 0.001* (intercept + (Double_t)nbins * slope);


  TH1D * h_data_calib =  new TH1D("","",nbins,e_min,e_max);
  for (Int_t ibin = 0; ibin < nbins; ibin++){
    h_data_calib->SetBinContent(ibin+1,h_data_raw->GetBinContent(ibin+1));
    h_data_calib->SetBinError(ibin+1,sqrt(h_data_raw->GetBinContent(ibin+1)));
  }

  c1 = new TCanvas("c1","c1");
  // h_data_calib->Draw();
  // return;
  c1->Print((output_filename+"[").c_str(),"Portrait");

  // Fit K40, Co60 and Th208 peaks before rebinning

  const Int_t n_sources = 3;
  Double_t source_e_true[n_sources]
    = {1.4608,
       // 2.6145,
       1.1732,
       1.3325}; // in MeV
  Double_t source_e_mean[n_sources];
  Double_t source_e_error[n_sources];

  for (Int_t i = 0; i < n_sources; i++){
    fit_peak_ge(h_data_calib,0.975*source_e_true[i] ,1.025*source_e_true[i],
    	     &source_e_mean[i], &source_e_error[i]);
    h_data_calib-> GetYaxis()-> SetTitle("Counts");
    h_data_calib->GetXaxis()->SetTitle("Energy [MeV]");

    h_data_calib->Draw();
    c1->Print(output_filename.c_str(),"Portrait");

  }
  

  if (data_type == 1 || data_type == 2){
    h_data_calib->Rebin(rebin_fac);
  }



  std::vector<double> chi2_vec;
  std::vector<double> smeared_chi2_vec;
  std::vector<double> x_vec;
  std::vector<double> p_vec;

  Double_t chi2_min[2] = {1e10,1e10}; //default->[0] smeared->[1]
  Int_t x_best[2] = {-1,-1};
  Double_t p_best[2] = {0,0};

  Char_t filename[256]; //mc root file

  gStyle->SetOptFit(0);
  for (Int_t x = x_min; x < x_max; x++){
    if (data_type == 0)
      sprintf(filename,"/disk02/usr6/sshima/Ge2021/mc/for_gedeo_oldGe/root/%5.5d.root",x);
    else if (data_type == 1)
      sprintf(filename,"mc/root_newGe/%5.5d.root",x);
    else if (data_type == 2)
      sprintf(filename,"mc/root_newGe_4db/%5.5d.root",x);

    TFile * f_mc = new TFile(filename);
    cout << filename << endl;

    TH1F * h_mc = (TH1F*)f_mc->Get("h21");//default mc
    Int_t mc_entry_bin = 0;
    TH1D * h_smc = new TH1D("smeared","smeared",nbins,e_min,e_max);//smeared mc
    Int_t bin_min = h_mc -> FindBin(fit_e_min) -20;
    Int_t bin_max = h_mc -> FindBin(fit_e_max) +10;

   //smearing mc
    for(i=bin_min;i<bin_max;i++){
      mc_entry_bin = h_mc-> GetBinContent(i+1);
      for(int j = 0; j<200; j++){
	double random = gRandom -> Gaus(h_mc->GetXaxis()->GetBinCenter(i+1),source_e_error[2]);
	double sf = double(mc_entry_bin)/200.0;
	  h_smc -> Fill(random,sf);
      }
    }

    if (data_type == 1 || data_type == 2){
      h_mc->Rebin(rebin_fac);
    }

    Double_t etot = 0.001 * (Double_t)x;
    Double_t p = sqrt(etot*etot-0.511*0.511);

    Double_t chi2[2] ={1e10,1e10};
    chi2[0] = calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, true, h_smc, false);
    chi2[1] = calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, true, h_smc, true);
    TLegend * leg = new TLegend(0.15,0.4,0.35,0.55);
     leg ->SetFillColor(0);
     leg -> AddEntry( h_mc, "Default MC" , "1");
    if(smear_flag){
      leg -> AddEntry( h_smc, "Smeared MC" , "1");
    }


    if(smear_flag){
    TText * t_etot = new TLatex (0.15, 0.8, Form("E_{tot} = %4.3f MeV",etot));
    TText * t_p = new TText (0.15, 0.7, Form("P = %4.3f MeV",p));
    TText * t_chi2 = new TLatex (0.15, 0.6, Form("#chi^{2} = %2.1f",chi2[1]));
    } else {
    TText * t_etot = new TLatex (0.15, 0.8, Form("E_{tot} = %4.3f MeV",etot));
    TText * t_p = new TText (0.15, 0.7, Form("P = %4.3f MeV",p));
    TText * t_chi2 = new TLatex (0.15, 0.6, Form("#chi^{2} = %2.1f",chi2[0]));
    }
    t_etot->SetNDC();
    t_p->SetNDC();
    t_chi2->SetNDC();
    t_etot->Draw();
    t_p->Draw();
    t_chi2->Draw();
    leg -> Draw();

    c1->Print(output_filename.c_str(),"Portrait");

    std::cout << x << " " << chi2[0] <<" "<<chi2[1]<< std::endl;

    x_vec.push_back((Double_t)x);
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
    f_mc->Close();
  }

  TCanvas * c2 = new TCanvas("c2","c2",600,600);
  if(smear_flag){
  TGraph * g = new TGraph (chi2_vec.size(),&x_vec[0],&smeared_chi2_vec[0]);
  g->GetXaxis()->SetRangeUser(x_best[1]-15,x_best[1]+15);
  g->GetYaxis()->SetRangeUser(0,chi2_min[1]+200);
  } else {
  TGraph * g = new TGraph (chi2_vec.size(),&x_vec[0],&smeared_chi2_vec[0]);
  g->GetXaxis()->SetRangeUser(x_best[0]-15,x_best[0]+15);
  g->GetYaxis()->SetRangeUser(0,chi2_min[0]+200);
  }

  g->SetTitle(";Total energy (keV);#chi^{2}");
  g->Draw("AL");
  c2->Print(output_filename.c_str(),"Portrait");

  if(smear_flag){
  g->GetXaxis()->SetRangeUser(x_best[1]-15,x_best[1]+15);
  g->GetYaxis()->SetRangeUser(0,chi2_min[1]+200);
  } else {
  g->GetXaxis()->SetRangeUser(x_best[0]-15,x_best[0]+15);
  g->GetYaxis()->SetRangeUser(0,chi2_min[0]+200);
  }
  g->Draw("AL");
  c2->Print(output_filename.c_str(),"Portrait");
  TCanvas * c3 = new TCanvas("c3","c3",600,600);
  if(smear_flag){
    TGraph * gp = new TGraph (smeared_chi2_vec.size(),&p_vec[0],&smeared_chi2_vec[0]);
  }else {
    TGraph * gp = new TGraph (chi2_vec.size(),&p_vec[0],&chi2_vec[0]);
  }

  gp->SetTitle(";Momentum (MeV);#chi^{2}");
  gp->Draw("AL");
  c3->Print(output_filename.c_str(),"Portrait");
  if(smear_flag){
  gp->GetXaxis()->SetRangeUser(p_best[1]-0.015,p_best[1]+0.015);
  gp->GetYaxis()->SetRangeUser(0,chi2_min[1]+200);
  }else{
    gp->GetXaxis()->SetRangeUser(p_best[0]-0.015,p_best[0]+0.015);
    gp->GetYaxis()->SetRangeUser(0,chi2_min[0]+200);
  }
  gp->Draw("AL");
  c3->Print(output_filename.c_str(),"Portrait");

  if (data_type == 0)
    if(smear_flag){
    sprintf(filename,"/disk02/usr6/sshima/Ge2021/mc/for_gedeo_oldGe/root/%5.5d.root",x_best[0]);
    } else {
    sprintf(filename,"/disk02/usr6/sshima/Ge2021/mc/for_gedeo_oldGe/root/%5.5d.root",x_best[1]);
    }
  else if (data_type == 1)
    sprintf(filename,"mc/root_newGe/%5.5d.root",x_best[0]);
  else if (data_type == 2)
    sprintf(filename,"mc/root_newGe_4db/%5.5d.root",x_best[0]);

  TFile * f_mc = new TFile(filename);
  TH1F * h_mc = (TH1F*)f_mc->Get("h21");
  TH1D * h_smc = new TH1D("smearing","smearing",nbins,e_min,e_max);
  Int_t bin_min = h_mc -> FindBin(fit_e_min) -20;
  Int_t bin_max = h_mc -> FindBin(fit_e_max) +10;

  for(i=bin_min;i<bin_max;i++){
    mc_entry_bin = h_mc-> GetBinContent(i+1);
    for(int j = 0; j<2000; j++){
	double random = gRandom -> Gaus(h_mc->GetXaxis()->GetBinCenter(i+1),source_e_error[2]);
	double sf = double(mc_entry_bin)/2000.0;
	h_smc -> Fill(random,sf);
      }
  }

  if (data_type == 1 || data_type == 2){
    h_mc->Rebin(rebin_fac);
  }
  Double_t chi2[2] = {1e10,1e10};
  chi2[0] = calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, true, h_smc, false);
  chi2[1] = calc_chi2(h_data_calib, h_mc, fit_e_min, fit_e_max, true, h_smc, true);
  cout << "Best fit total energy = " << x_best[0] << " (keV)" <<  endl;
  cout << "Best fit momentum = " << p_best[0] << " (MeV)" <<  endl;
  cout << "Best fit total energy = " << x_best[1] << " (keV)" <<  endl;
  cout << "Best fit momentum = " << p_best[1] << " (MeV)" <<  endl;

  if(smear_flag){
    TText * t_etot = new TLatex (0.15, 0.8, Form("E_{tot} = %4.3f MeV",0.001*x_best[1]));
    TText * t_p = new TText (0.15, 0.7, Form("P = %4.3f MeV",p_best[1]));
    TText * t_chi2 = new TLatex (0.15, 0.6, Form("#chi^{2} = %2.1f",chi2_min[1]));
  } else{
    TText * t_etot = new TLatex (0.15, 0.8, Form("E_{tot} = %4.3f MeV",0.001*x_best[0]));
    TText * t_p = new TText (0.15, 0.7, Form("P = %4.3f MeV",p_best[0]));
    TText * t_chi2 = new TLatex (0.15, 0.6, Form("#chi^{2} = %2.1f",chi2_min[0]));
  }
  TLegend * leg = new TLegend(0.15,0.4,0.35,0.55);
  leg -> SetFillColor(0);
  leg -> AddEntry( h_mc, "Default MC" , "l");
  if(smear_flag){
  leg -> AddEntry( h_smc, "Smeared MC" , "l");
  }
  t_etot->SetNDC();
  t_p->SetNDC();
  t_chi2->SetNDC();
  t_etot->Draw();
  t_p->Draw();
  t_chi2->Draw();
  leg -> Draw();

  c1->Print(output_filename.c_str(),"Portrait");

  c1->Print((output_filename+"]").c_str(),"Portrait");

  ofstream ofs ("bestfit_momentum_tmp.txt");
  ofstream diff_out ("diff_default_smeared_w_new_nicalib.txt",ios::app);
  if(smear_flag){
    ofs << p_best[1];
  }else {
    ofs << p_best[0];
  }
  diff_out << xpos << " \t" << zpos <<"\t"<< beam_energy << "\t";
  diff_out << 0.001*x_best[0] << "\t" << 0.001*x_best[1] << "\t" << chi2[0] << "\t" << chi2[1] << std::endl;

  for (Int_t i = 0; i < n_sources; i++){
    ofs << "\t" << source_e_mean[i] << "\t" << source_e_error[i];
  }
  ofs << endl;
  ofs.close();
  diff_out.close();

  // h_mc->Draw();

}

Double_t calc_chi2(TH1D * h_data, TH1F * h_mc, Double_t e_min, Double_t e_max, Bool_t plot_flag = false,TH1D * h_smc, Bool_t flag = false){
  Int_t bin_min = h_data->FindBin(e_min);
  Int_t bin_max = h_data->FindBin(e_max);
  cout << "================ " << e_min << " " << e_max << " " << bin_min << " " << bin_max << endl;
  Double_t norm_data = h_data->Integral(bin_min,bin_max);
  Double_t norm_mc = h_mc->Integral(bin_min,bin_max);
  Double_t norm_smc = h_smc->Integral(bin_min,bin_max);
  // consistency check:
  Double_t xmax_data = h_data->GetXaxis()->GetXmax();
  Double_t xmax_mc = h_mc->GetXaxis()->GetXmax();
  cout << "consistency check: " << xmax_data << " " << xmax_mc << endl;
  cout << "consistency check: " << h_data->GetXaxis()->GetXmin()
       << " " << h_mc->GetXaxis()->GetXmin() << endl;


  Double_t hist_x_min = e_min - 0.1;
  Double_t hist_x_max = e_max + 0.05;


  Double_t chi2 = 0;
  Double_t n_data = 0;
  Double_t n_mc = 0;
  Double_t n_smc = 0;
  Double_t err_data = 0;
  Double_t err_mc = 0;
  Double_t err_smc = 0;

  Int_t ndf = 0;

  for (Int_t i = bin_min; i <= bin_max; i++){
    n_data = h_data->GetBinContent(i);
    if (n_data == 0) continue; // skip if number of events in data is 0
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
    //    cout << i << " " << n_data << " " << n_mc << endl;
  }
  Int_t max_bin_mc = h_mc -> GetMaximumBin();
  Double_t scalefactor;
  //h_smc -> Scale(scalefactor);

  string dummy;
  if (plot_flag){
    Int_t nbins = h_mc->GetNbinsX();
    for (Int_t i = 0; i < nbins; i++){
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
    cout << "chi2/NDF = " << chi2 << " / " << ndf-1 << endl;
  }

  return chi2;
}
