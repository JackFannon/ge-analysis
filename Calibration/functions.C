
void read_old_ge_raw_data(std::string filename = "", TH1D* h = 0){
    ifstream ifs(filename.c_str());
    if(!ifs.is_open()){
        cout << "!!!!! Cannot find " << filename.c_str() << "!!!!!" << endl;
        return;
    }
    string buf;
    int ibin = 0;
    int nbins = h->GetNbinsX();

    // Need to skip the first 9 lines of the data-file. None of these lines contain data
    int linesToSkip = 9;

    for(int i = 0; i < linesToSkip; i++){
        getline(ifs, buf);
    }

    string dummy1;

    while (!ifs.eof()){
        getline(ifs, buf);
        std::stringstream ss(buf);
        ss >> dummy1;

        for(int i = 0; i < 16; i++){
            int count = 0;
            ss >> count;
            ibin++;
            if(ibin <= nbins){
                h->SetBinContent(ibin, count);
            }
        }
    }

}

void fit_peak_ge(TH1D* h_in, double min, double max, double *mean, double *error){
    TF1* fit = new TF1("gausexpo", "gaus(0)", min - 20, max + 20);

    cout << "min = " << min << " max = " << max << endl;
}
